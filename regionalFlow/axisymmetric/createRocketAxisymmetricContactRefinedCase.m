function item=createRocketAxisymmetricContactRefinedCase(options)
%CREATEROCKETAXISYMMETRICCONTACTREFINEDCASE 从已收敛HLLC72根建立同方程144层初猜。
% 模块一：保留原144层几何与物理边界；明确区分同格HLLE初猜的失败试验。
% 模块二：五量凸插值只作初猜，不裁温度/SA，不放宽任何星态或稳态门槛。
% 模块三：新根和新BE链独立绑定，不续写或改判原失败试验。
if nargin<1,options=struct;end
defaults=struct('rootEvidenceSHA256','');assert(isstruct(options)&&isscalar(options)&& ...
    all(ismember(fieldnames(options),fieldnames(defaults))),'RocketContactRefined:Options','只允许新根SHA。');
if ~isfield(options,'rootEvidenceSHA256'),options.rootEvidenceSHA256='';end
expected=char(options.rootEvidenceSHA256);assert(isempty(expected)||~isempty(regexp(expected,'^[a-f0-9]{64}$','once')), ...
    'RocketContactRefined:Evidence','新根SHA格式不符。');
parentHash='8e9bc2894d43a85a618b2fbe3e01d9d912a15ef9c91e0654e4fa6b2ba729b372';
parent=createRocketAxisymmetricContactCase(struct('kind','collar72','rootEvidenceSHA256',parentHash));
item=createRocketAxisymmetricContactCase(struct('kind','nested144','rootEvidenceSHA256',''));
[item.Q,mapping]=mapRocketAxisymmetricInitialGuess(parent,parent.Q,item);
item.initialization.kind='internal-converged-HLLC72-root-convex-interpolation-for-same-equation-fine-grid';
item.initialization.mapping=mapping;item.initialization.seed=struct('kind','v78-HLLC72-root-for-fine-grid', ...
    'sourceEvidenceSHA256',parentHash,'parentHLLCPreparationNewtonIterations',parent.initialization.seed.priorNewtonIterations, ...
    'priorAcceptedSteps',0,'initialGuessIsBoundaryCondition',false,'externalCFDImported',false,'physicalSolution',false);
item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
item.inputContract.internalInitialGuess=item.initialization.seed;item.contactStudy.kind='nested144_from_hllc72';
item.contactStudy.sameGridAndConservativeInitialQ=false;item.contactStudy.samePhysicalGridAndBoundaryAsFailedDirect144=true;
item.contactStudy.initializationIsConservativeRemapping=false;item.contactStudy.coarseHLLCRootSHA256=parentHash;
item.inputContract.contactStudy=item.contactStudy;item.caseFactory='createRocketAxisymmetricContactRefinedCase';item.options=options;
if ~isempty(expected)
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));path=fullfile(root,'validation','v78','nested144_from_hllc72','spatial_root_evidence.mat');
    if ~exist(path,'file'),path=fullfile(root,'outputs','v78_validation','nested144_from_hllc72','spatial_root_evidence.mat');end
    assert(exist(path,'file')&&strcmp(fileSHA(path),expected),'RocketContactRefined:Evidence','新根胶囊缺失或改字节。');
    d=load(path,'c','solution','sourceBindings');original=item;original.options.rootEvidenceSHA256='';
    assert(isequaln(d.c,original)&&isequaln(d.sourceBindings,rocketContactRefinedSourceBindings)&& ...
        d.solution.spatialRootConverged&&isequal(d.solution.initialQ,item.Q)&& ...
        d.solution.options.steadyResidualTolerance==2e-6&&d.solution.options.conservationTolerance==5e-10, ...
        'RocketContactRefined:Root','新根的实际工况、初猜、源码或原门槛不一致。');
    item.Q=d.solution.Q;item.initialization.kind='internal-fine-grid-HLLC-spatial-root-for-independent-original-BE-confirmation';
    item.initialization.seed.kind='v78-HLLC-refined-spatial-root';item.initialization.seed.sourceEvidenceSHA256=expected;
    item.initialization.seed.coarseHLLCRootSHA256=parentHash;item.initialization.seed.priorNewtonIterations=d.solution.acceptedNewtonIterations;
    item.initialization.spatialRootHistory=d.solution.history;
    item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
    item.inputContract.internalInitialGuess=item.initialization.seed;item.inputContract.initializationSolver='solveAxisymmetricSteadyNewton';
end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
