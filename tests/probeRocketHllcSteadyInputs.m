function report=probeRocketHllcSteadyInputs
%PROBEROCKETHLLCSTEADYINPUTS 在相同有效定常Q上检查已有HLLC及完整五量耦合。
% 模块一：严格核对v77真实根SHA/依赖，不能重复使用v72早期失败状态。
% 模块二：只改内部Euler离散；物理边界、黏性、几何源与SA状态不改。
% 模块三：原星态护栏拒绝时保留真实错误及输入，不裁值/不回退，也不算求解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
folder=fullfile(root,'outputs','v78_validation','initial_hllc');assert(~exist(folder,'dir'),'RocketHllcProbe:Overwrite','禁止覆盖真实初态检查。');
mkdir(folder);rows=struct([]);
kinds={'collar72','nested144'};expected={'fc47b5dad82a379cc9b8c57eaa921029c5c300bf67dd179b815036e03cc5a1e3','d33a95de11609e08b2e919d4c4aea26d8d83c28314871515d3c22ff7cf9d591d'};
for k=1:2
    path=fullfile(root,'outputs','v77_validation',kinds{k},'spatial_root_evidence.mat');
    if ~exist(path,'file'),path=fullfile(root,'validation','v77',kinds{k},'spatial_root_evidence.mat');end
    assert(strcmp(fileSHA(path),expected{k}),'RocketHllcProbe:Evidence','原有效定常根改字节。');
    d=load(path);assert(d.solution.spatialRootConverged&&isequaln(d.sourceBindings,rocketWallResolutionSourceBindings));
    c=d.c;Q=d.solution.Q;before=c.boundary;c.boundary.interiorConvectiveFlux='hllc';
    physical=c.boundary;physical.interiorConvectiveFlux=before.interiorConvectiveFlux;assert(isequal(physical,before));
    succeeded=false;identifier='';message='';norm=NaN;rate=[];info=[];
    try
        [rate,info]=axisymmetricRANSResidual(Q,c.mesh,c.gas,c.boundary);succeeded=true;
        norm=max(abs(rate)./c.residualScale,[],'all');
        old=d.solution.finalInfo;
        assert(isequal(info.faceViscousFluxPerRadian,old.faceViscousFluxPerRadian)&& ...
            isequal(info.cellSourceIntegral,old.cellSourceIntegral)&&isequaln(info.wall,old.wall)&& ...
            isequal(info.faceConvectiveFluxPerRadian(c.mesh.boundaryFaces,:),old.faceConvectiveFluxPerRadian(c.mesh.boundaryFaces,:))&& ...
            isequal(info.sa.faceDiffusiveFluxPerRadian,old.sa.faceDiffusiveFluxPerRadian), ...
            'RocketHllcProbe:Controlled','同Q黏性、墙、源或实际开放边界发生了变化。');
    catch ex
        identifier=ex.identifier;message=ex.message;
        if ~strcmp(identifier,'AxisFV:HLLC'),rethrow(ex);end
    end
    row=struct('kind',kinds{k},'sourceEvidenceSHA256',expected{k},'originalHLLESteadyResidual',d.solution.residualNorm, ...
        'HLLCAtSameQAcceptedByOriginalGuards',succeeded,'HLLCOriginalResidualAtSameQ',norm,'errorIdentifier',identifier, ...
        'errorMessage',message,'physicalBoundaryChanged',false,'stateClipped',false,'fallbackUsed',false, ...
        'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
    save(fullfile(folder,[kinds{k} '_initial_hllc.mat']),'c','Q','rate','info','row','-v7');
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    fprintf('%s\n',jsonencode(row));
end
report=struct('status','completed-same-valid-steady-Q-original-HLLC-guard-and-complete-RANS-coupling-check', ...
    'cases',rows,'physicalSolution',false,'newProductionSolvedAreaM2',0);
fid=fopen(fullfile(folder,'initial_hllc_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
