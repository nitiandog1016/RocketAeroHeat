function item=createRocketAxisymmetricWallResolutionCase(options)
%CREATEROCKETAXISYMMETRICWALLRESOLUTIONCASE 有界的近壁/整径向分辨率对照。
% 模块一：72层中20个法向层，2μm首节点偏移、1.25增长；保持同一真实墙/外域。
% 模块二：144层为72层逐格二分，可显式从内部72层根作数值初猜。
% 模块三：插值、直接Newton、三步BE确认分开；未取得物理/实验资格。
if nargin<1,options=struct;end
defaults=struct('kind','collar72','parent72SHA256','','rootEvidenceSHA256','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'RocketWallResolution:Options','只接受明确网格类型与内部胶囊SHA。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
kind=char(options.kind);parentHash=char(options.parent72SHA256);rootHash=char(options.rootEvidenceSHA256);
assert(ismember(kind,{'collar72','nested144'}),'RocketWallResolution:Kind','仅支持72层及其逐格二分144层。');
for hash={parentHash,rootHash},assert(isempty(hash{1})||~isempty(regexp(hash{1},'^[a-f0-9]{64}$','once')), ...
        'RocketWallResolution:Evidence','胶囊SHA必须为64位小写十六进制。');end
assert(isempty(parentHash)||strcmp(kind,'nested144'),'RocketWallResolution:Parent','72层工况不能给自身作为父根。');
referenceHash='ed4faf40c2b01f35560ca3d5033d96bc9267ebf6cc68d98be722538f54cfbbcd';
reference=createRocketAxisymmetricSensitivityCase(struct('kind','radial48','rootEvidenceSHA256',referenceHash));
domain=reference.domain.options;domain.radialLayers=72;
collar=struct('normalLayers',20,'firstNormalSpacing',2e-6,'normalGrowth',1.25);
base=struct('initialGuess','uniform','domain',domain,'collar',collar,'flight',reference.inputContract.flight, ...
    'solverInitialPseudoCFL',5,'solverMaximumPseudoCFL',1e6);
item=createRocketAxisymmetricLocalCase(base);parent=reference;parentQ=reference.Q;seedHash=referenceHash;
refinement=struct('kind','none','parentCellCount',item.mesh.cellCount,'childCellCount',item.mesh.cellCount);
if strcmp(kind,'nested144')
    parent72=item;
    [item.mesh,item.domain,refinement]=refineRocketAxisymmetricRadial(parent72.mesh,parent72.domain);
    item.Q=repmat(item.Q(1,:),item.mesh.cellCount,1);
    if ~isempty(parentHash)
        d=loadCapsule('collar72',parentHash);assert(d.solution.spatialRootConverged&& ...
            isequal(d.c.mesh,parent72.mesh)&&isequal(d.c.boundary,parent72.boundary)&& ...
            isequaln(d.sourceBindings,rocketWallResolutionSourceBindings), ...
            'RocketWallResolution:Parent','内部72层父根与真实网格、方程或实际源码不一致。');
        parent=d.c;parentQ=d.solution.Q;seedHash=parentHash;
    end
end
assert(isequaln(item.surfaceMesh,reference.surfaceMesh)&&isequal(item.boundary,reference.boundary)&& ...
    isequal(item.gas,reference.gas)&&isequal(item.domain.wallPolylineXR,reference.domain.wallPolylineXR), ...
    'RocketWallResolution:Inputs','不能改变原生几何、气体或物理边界。');
[item.Q,mapping]=mapRocketAxisymmetricInitialGuess(parent,parentQ,item);
study=struct('kind',kind,'reference48RootSHA256',referenceHash,'initialGuessSourceSHA256',seedHash, ...
    'cellCount',item.mesh.cellCount,'effectiveDomainOptions',item.domain.options, ...
    'effectiveCollarOptions',item.domain.collar.options,'radialRefinement',refinement, ...
    'sameWallAndSourceSurface',true,'streamwiseRefinementPerformed',false,'gridIndependentClaimed',false, ...
    'experimentalValidation',false,'newProductionSolvedAreaM2',0);
item.initialization.kind='internal-root-on-new-wall-resolution-grid-numerical-initial-guess';
item.initialization.allFiveComponentsUniform=false;item.initialization.mapping=mapping;
item.initialization.seed=struct('kind','v77-wall-resolution-interpolated','sourceEvidenceSHA256',seedHash, ...
    'priorAcceptedSteps',0,'initialGuessIsBoundaryCondition',false,'externalCFDImported',false,'physicalSolution',false);
item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
item.caseFactory='createRocketAxisymmetricWallResolutionCase';item.options=options;item.wallResolution=study;
item.inputContract.wallResolution=study;item.inputContract.internalInitialGuess=item.initialization.seed;
item.inputContract.effectiveDomainOptions=item.domain.options;item.inputContract.effectiveCollarOptions=item.domain.collar.options;
item.inputContract.interiorFluxArithmetic='double';
if ~isempty(rootHash)
    d=loadCapsule(kind,rootHash);original=item;original.options.rootEvidenceSHA256='';
    assert(isequaln(d.c,original)&&isequaln(d.sourceBindings,rocketWallResolutionSourceBindings)&& ...
        d.solution.spatialRootConverged&&isequal(d.solution.initialQ,item.Q)&& ...
        d.solution.options.steadyResidualTolerance==2e-6&&d.solution.options.conservationTolerance==5e-10, ...
        'RocketWallResolution:Root','新根的初猜、真实工况、依赖或原门槛不成立。');
    item.Q=d.solution.Q;item.initialization.kind='internal-wall-resolution-spatial-root-for-separate-BE-confirmation';
    item.initialization.seed.kind='v77-wall-resolution-root';item.initialization.seed.sourceEvidenceSHA256=rootHash;
    item.initialization.seed.priorNewtonIterations=d.solution.acceptedNewtonIterations;
    item.initialization.spatialRootHistory=d.solution.history;
    item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
    item.inputContract.internalInitialGuess=item.initialization.seed;item.inputContract.initializationSolver='solveAxisymmetricSteadyNewton';
end
end
function d=loadCapsule(kind,expected)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
path=fullfile(root,'validation','v77',kind,'spatial_root_evidence.mat');
if ~exist(path,'file'),path=fullfile(root,'outputs','v77_validation',kind,'spatial_root_evidence.mat');end
fid=fopen(path,'rb');assert(fid>=0,'RocketWallResolution:Evidence','缺少明确指定的内部根胶囊。');
guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
assert(strcmp(hash,expected),'RocketWallResolution:Evidence','内部根胶囊改字节。');d=load(path);
end
