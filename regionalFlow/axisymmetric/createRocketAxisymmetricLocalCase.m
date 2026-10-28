function item=createRocketAxisymmetricLocalCase(options)
%CREATEROCKETAXISYMMETRICLOCALCASE 局部Jacobian新链的完整箭体工况及初猜身份。
% 模块一：继承v73真实完整几何、HLLE/远场/Tw/SA；数值配置显式保留。
% 模块二：均匀初猜或内部冻结Q23可选；后者完整校验SHA/几何/气体/边界，不作为给定入口。
% 模块三：初猜历史和本链接受步分开；未知Te/Tr/求解Tw及物理解资格不改变。
if nargin<1,options=struct;end
assert(isstruct(options)&&isscalar(options),'RocketLocalCase:Options','工况须为单一结构。');
base=options;guess='v73-uniform23';relative=2e-6;initialCFL=[];maximumCFL=1e8;
keys={'initialGuess','solverDifferenceStep','solverInitialPseudoCFL','solverMaximumPseudoCFL'};
for k=1:numel(keys)
    key=keys{k};if ~isfield(base,key),continue;end;value=base.(key);base=rmfield(base,key);
    switch key
        case 'initialGuess',guess=value;
        case 'solverDifferenceStep',relative=value;
        case 'solverInitialPseudoCFL',initialCFL=value;
        case 'solverMaximumPseudoCFL',maximumCFL=value;
    end
end
assert(((ischar(guess)&&isrow(guess))||(isstring(guess)&&isscalar(guess)&&~ismissing(guess)))&& ...
    ismember(char(guess),{'uniform','v73-uniform23'}),'RocketLocalCase:Guess','只支持均匀初猜或显式冻结Q23。');guess=char(guess);
if isempty(initialCFL),if strcmp(guess,'uniform'),initialCFL=.5;else,initialCFL=1e8;end,end
for value={relative,initialCFL,maximumCFL},validateattributes(value{1},{'numeric'},{'real','finite','positive','scalar'});end
assert(initialCFL<=maximumCFL,'RocketLocalCase:CFL','初始CFL不能大于上限。');
item=createRocketAxisymmetricFluxCase(base);
seed=struct('kind','uniform','sourceEvidenceSHA256','','priorAcceptedSteps',0, ...
    'initialGuessIsBoundaryCondition',false,'externalCFDImported',false,'physicalSolution',false);
if strcmp(guess,'v73-uniform23')
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    candidates={fullfile(root,'validation','v73','frozen_uniform23_evidence.mat'), ...
        fullfile(root,'outputs','v73_validation','frozen_uniform23_evidence.mat')};path='';
    for k=1:numel(candidates),if exist(candidates{k},'file'),path=candidates{k};break;end,end
    assert(~isempty(path),'RocketLocalCase:Evidence','缺少内部冻结Q23，不能退回均匀初猜。');
    expected='915de3818d4f8a06056356f71fd2f6114ad843e511df106c4e7565707f40060b';
    assert(strcmp(fileSHA(path),expected),'RocketLocalCase:Evidence','冻结初猜证据字节不匹配。');
    old=load(path,'caseDefinition','implicit');c=old.caseDefinition;
    assert(old.implicit.acceptedSteps==23&&~old.implicit.physicalSolution&& ...
        isequaln(item.mesh,c.mesh)&&isequaln(item.gas,c.gas)&&isequaln(item.boundary,c.boundary)&& ...
        isequal(item.residualScale,c.residualScale)&&isequal(item.Q,old.implicit.initialQ), ...
        'RocketLocalCase:Inputs','冻结初猜与当前几何、气体、边界或尺度不一致。');
    item.Q=old.implicit.Q;seed.kind=guess;seed.sourceEvidenceSHA256=expected;seed.priorAcceptedSteps=23;
    item.initialization.kind='internal-frozen-complete-domain-numerical-initial-guess';
    item.initialization.allFiveComponentsUniform=false;
    item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
end
item.initialization.seed=seed;item.initialization.timeStepsTaken=0;
item.solverConfiguration=struct('maxSteps',1,'initialPseudoCFL',initialCFL,'maximumPseudoCFL',maximumCFL,'pseudoCFLGrowth',2, ...
    'maxNewtonIterations',24,'jacobianRefreshEvery',1,'differenceScheme','central','differenceStep',relative, ...
    'residualScale',item.residualScale,'steadyResidualTolerance',2e-6,'steadyRequiredSteps',3, ...
    'minimumSteadySteps',3,'conservationTolerance',5e-10,'nonlinearAbsoluteTolerance',1e-11,'nonlinearRelativeTolerance',1e-9);
item.caseFactory='createRocketAxisymmetricLocalCase';item.options=options;
item.inputContract.jacobianDifferencePolicy='local-conservative-acoustic-SA-roundoff-floor-v1';
item.inputContract.jacobianRelativeStep=relative;item.inputContract.internalInitialGuess=seed;
item.inputContract.newChainIncludesPriorAcceptedSteps=false;item.inputContract.changedFromV73= ...
    'same full spatial RANS/SA and physical boundaries; explicit local-conservative numerical Jacobian and internally bound initial guess';
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
