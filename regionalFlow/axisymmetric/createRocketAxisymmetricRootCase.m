function item=createRocketAxisymmetricRootCase(options)
%CREATEROCKETAXISYMMETRICROOTCASE 从内部直接定常解建立独立BE确认链工况。
% 模块一：校验真实定常胶囊、计算源码、完整几何/气体/边界；不是导入外部CFD。
% 模块二：Newton准备和本链BE步单列，原Q6/此前Q23不计入新链历史。
% 模块三：给定Tw、未知Te/Tr/求解Tw及物理/面积资格保持；普通算术为默认。
if nargin<1,options=struct;end
assert(isstruct(options)&&isscalar(options),'RocketRootCase:Options','工况须为单一结构。');base=options;arithmetic='double';
if isfield(base,'rootArithmetic'),arithmetic=base.rootArithmetic;base=rmfield(base,'rootArithmetic');end
assert(((ischar(arithmetic)&&isrow(arithmetic))||(isstring(arithmetic)&&isscalar(arithmetic)&&~ismissing(arithmetic)))&& ...
    ismember(char(arithmetic),{'double','two-double-euler'}),'RocketRootCase:Arithmetic','未知根求解算术。');arithmetic=char(arithmetic);
if isfield(base,'initialGuess')
    assert(strcmp(base.initialGuess,'v75-spatial-root'),'RocketRootCase:Guess','此工况须明确使用内部定常根。');base=rmfield(base,'initialGuess');
end
base.initialGuess='uniform';
if ~isfield(base,'solverInitialPseudoCFL'),base.solverInitialPseudoCFL=5;end
if ~isfield(base,'solverMaximumPseudoCFL'),base.solverMaximumPseudoCFL=1e6;end
item=createRocketAxisymmetricLocalCase(base);root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
tag=strrep(arithmetic,'-','_');relative=['spatial_root_' tag '/spatial_root_evidence.mat'];
if strcmp(arithmetic,'double'),expected='44100288ee800df1c2de3f07cb68bfe62979a8dcc00cf870eaf9adf5c599b2b5';
else,expected='651b17529cec5c2ec32afa588fcfe9d262c394255e98e2a2879e5485cc6dab3a';end
folder=fullfile(root,'validation','v75');if ~exist(fullfile(folder,relative),'file'),folder=fullfile(root,'outputs','v75_validation');end
path=fullfile(folder,relative);assert(exist(path,'file')&&strcmp(fileSHA(path),expected),'RocketRootCase:Evidence','内部定常胶囊缺失或字节不匹配。');
d=load(path,'solution','confirmation','row','c');actual=d.c.boundary;bare=actual;
if isfield(bare,'interiorFluxArithmetic'),bare=rmfield(bare,'interiorFluxArithmetic');end
assert(isequaln(item.mesh,d.c.mesh)&&isequaln(item.gas,d.c.gas)&&isequaln(item.boundary,bare)&& ...
    isequal(item.residualScale,d.c.residualScale)&&d.solution.spatialRootConverged&& ...
    d.solution.options.steadyResidualTolerance==2e-6&&d.solution.options.conservationTolerance==5e-10&& ...
    d.confirmation.acceptedSteps==3&&d.confirmation.steadyConverged, ...
    'RocketRootCase:Inputs','胶囊与实际物理输入不同或未经过原标准求解/确认。');
% 只核对实际求根依赖，旧的无关导出/测试代码另存来源，不重绑旧链。
process=jsondecode(fileread(fullfile(folder,'spatial_root01_process.json')));
assert(strcmp(process.status,'terminal')&&process.exitCode==0,'RocketRootCase:Process','内部求根没有正常终态。');
needed={'solveAxisymmetricSteadyNewton','axisymmetricFVJacobianLocal','axisymmetricJacobianLocalReference', ...
    'axisymmetricRANSResidual','axisymmetricNavierStokesResidual','axisymmetricFVState','axisymmetricFVGas', ...
    'axisymmetricHLLEFlux','axisymmetricHLLEFluxAccurate','axisymmetricDoubleDouble'};
for k=1:numel(needed)
    name=['regionalFlow/axisymmetric/' needed{k} '.m'];key=matlab.lang.makeValidName(name);
    assert(isfield(process.sourcesBefore,key)&&strcmp(process.sourcesBefore.(key),fileSHA(fullfile(root,name))), ...
        'RocketRootCase:Source','生成定常根的实际源码已变化。');
end
item.Q=d.solution.Q;item.boundary=actual;item.caseFactory='createRocketAxisymmetricRootCase';item.options=options;
seed=struct('kind','v75-spatial-root','sourceEvidenceSHA256',expected,'priorAcceptedSteps',0, ...
    'priorNewtonIterations',d.solution.acceptedNewtonIterations,'parentSeedKind','v74-local6','parentAcceptedSteps',6, ...
    'parentSeedPriorAcceptedSteps',23,'initialGuessIsBoundaryCondition',false,'externalCFDImported',false, ...
    'spatialRootResidualNorm',d.solution.residualNorm,'rootSolver','solveAxisymmetricSteadyNewton','rootArithmetic',arithmetic, ...
    'physicalSolution',false);
item.initialization.kind='internal-direct-spatial-Newton-root-not-physical-certification';
item.initialization.allFiveComponentsUniform=false;item.initialization.timeStepsTaken=0;item.initialization.seed=seed;
item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
item.initialization.spatialRootHistory=d.solution.history;
item.inputContract.internalInitialGuess=seed;item.inputContract.interiorFluxArithmetic=arithmetic;
item.inputContract.initializationSolver='solveAxisymmetricSteadyNewton';item.inputContract.newChainIncludesPriorAcceptedSteps=false;
item.inputContract.changedFromV74='same spatial RANS/SA; actual direct spatial Newton root used as internal initial guess; separate original BE confirmation';
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
