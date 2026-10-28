function item=createRocketAxisymmetricContactCase(options)
%CREATEROCKETAXISYMMETRICCONTACTCASE 同网格/同物理输入的已有HLLC受控工况。
% 模块一：显式载入v77的真实HLLE有效根，只作为内部数值初猜，不导入外部CFD。
% 模块二：仅改内部Euler通量为原HLLC，黏性/墙/远场/源/SA合同和所有护栏保持。
% 模块三：HLLC新根、父HLLE历史与新BE确认分别记录，不重绑或覆盖旧链。
if nargin<1,options=struct;end
defaults=struct('kind','collar72','rootEvidenceSHA256','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'RocketContact:Options','只接受明确网格类型和新根SHA。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
kind=char(options.kind);expected=char(options.rootEvidenceSHA256);
assert(ismember(kind,{'collar72','nested144'}),'RocketContact:Kind','只比较已有72/144层实际箭体。');
assert(isempty(expected)||~isempty(regexp(expected,'^[a-f0-9]{64}$','once')),'RocketContact:Evidence','根SHA须为64位小写十六进制。');
hash72='fc47b5dad82a379cc9b8c57eaa921029c5c300bf67dd179b815036e03cc5a1e3';
hash144='d33a95de11609e08b2e919d4c4aea26d8d83c28314871515d3c22ff7cf9d591d';
parentOptions=struct('kind',kind,'parent72SHA256','','rootEvidenceSHA256',hash72);parentHash=hash72;
if strcmp(kind,'nested144'),parentOptions.parent72SHA256=hash72;parentOptions.rootEvidenceSHA256=hash144;parentHash=hash144;end
parent=createRocketAxisymmetricWallResolutionCase(parentOptions);item=parent;
item.boundary.interiorConvectiveFlux='hllc';physical=item.boundary;physical.interiorConvectiveFlux=parent.boundary.interiorConvectiveFlux;
assert(isequal(physical,parent.boundary)&&strcmp(parent.boundary.interiorConvectiveFlux,'hlle-einfeldt'), ...
    'RocketContact:Inputs','受控试验不得更改物理边界或偷换父通量。');
item.initialization.kind='same-grid-internal-HLLE-spatial-root-for-explicit-HLLC-numerical-initial-guess';
if isfield(item.initialization,'spatialRootHistory'),item.initialization=rmfield(item.initialization,'spatialRootHistory');end
if isfield(item.initialization,'mapping'),item.initialization=rmfield(item.initialization,'mapping');end
seed=struct('kind','v77-HLLE-steady-root-for-HLLC','sourceEvidenceSHA256',parentHash, ...
    'priorAcceptedSteps',0,'parentHLLEPreparationNewtonIterations',parent.initialization.seed.priorNewtonIterations, ...
    'initialGuessIsBoundaryCondition',false,'externalCFDImported',false,'physicalSolution',false);
item.initialization.seed=seed;item.initialization.timeStepsTaken=0;
item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
study=struct('kind',kind,'parentHLLERootSHA256',parentHash,'parentFactory','createRocketAxisymmetricWallResolutionCase', ...
    'sameGridAndConservativeInitialQ',true,'physicalBoundaryChanged',false,'stateClipped',false,'fallbackUsed',false, ...
    'interiorConvectiveFlux','hllc','parentInteriorConvectiveFlux','hlle-einfeldt', ...
    'physicalSolution',false,'experimentalValidation',false,'gridIndependentClaimed',false,'newProductionSolvedAreaM2',0);
item.contactStudy=study;item.inputContract.contactStudy=study;
item.inputContract.meshProvenance=parent.wallResolution;
item.inputContract=rmfield(item.inputContract,'wallResolution');item=rmfield(item,'wallResolution');
item.inputContract.internalInitialGuess=seed;item.inputContract.interiorConvectiveFluxPolicy='hllc';
item.inputContract.initializationSolver='none-for-HLLC-initial-guess-only';
item.inputContract.changedFromV77='same real viscous rocket grid and physical inputs; original HLLC internal flux only; explicit old HLLE steady Q as internal initial guess';
item.caseFactory='createRocketAxisymmetricContactCase';item.options=options;
if ~isempty(expected)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    path=fullfile(root,'validation','v78',kind,'spatial_root_evidence.mat');
    if ~exist(path,'file'),path=fullfile(root,'outputs','v78_validation',kind,'spatial_root_evidence.mat');end
    assert(exist(path,'file')&&strcmp(fileSHA(path),expected),'RocketContact:Evidence','新HLLC根胶囊缺失或改字节。');
    d=load(path,'solution','c','sourceBindings');original=item;original.options.rootEvidenceSHA256='';
    assert(isequaln(d.c,original)&&isequaln(d.sourceBindings,rocketContactSourceBindings)&& ...
        d.solution.spatialRootConverged&&isequal(d.solution.initialQ,item.Q)&& ...
        d.solution.options.steadyResidualTolerance==2e-6&&d.solution.options.conservationTolerance==5e-10, ...
        'RocketContact:Root','新HLLC根的真实工况、初猜、源码或原门槛不符。');
    item.Q=d.solution.Q;item.initialization.kind='internal-HLLC-spatial-root-for-independent-original-BE-confirmation';
    item.initialization.seed.kind='v78-HLLC-spatial-root';item.initialization.seed.sourceEvidenceSHA256=expected;
    item.initialization.seed.parentHLLERootSHA256=parentHash;
    item.initialization.seed.priorNewtonIterations=d.solution.acceptedNewtonIterations;
    item.initialization.spatialRootHistory=d.solution.history;
    item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
    item.inputContract.internalInitialGuess=item.initialization.seed;item.inputContract.initializationSolver='solveAxisymmetricSteadyNewton';
end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
