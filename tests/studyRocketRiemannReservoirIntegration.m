function report=studyRocketRiemannReservoirIntegration(folder)
%STUDYROCKETRIEMANNRESERVOIRINTEGRATION 在真实Q19比较旧特征与显式精确储库的完整三维NS。
% 模块一：旧完整rate/info必须逐字段复现；旧历史不并入改变边界后的新方程。
% 模块二：只替换开放边界类型，原19752面/217272CV、壁温、物性及内部精确通量保持。
% 模块三：重新计算一个旧被拒候选，任何新错误均保留；合法候选不提升为接受或物理解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','krylov19','wall_velocity_krylov_evidence.mat');s=load(sourceFile,'result','gas','boundary');
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;
Q=s.result.Q;gas=s.gas;oldBoundary=s.boundary;scale=s.result.residualScale;
[legacyRate,legacyInfo]=cutCellNavierStokesResidual(Q,cut,gas,oldBoundary);
assert(isequaln(legacyRate,s.result.rate)&&isequaln(legacyInfo,s.result.info),'NSRiemannIntegration:Legacy','新显式选项破坏了旧完整NS路径。');
boundary=oldBoundary;
for key={'xmin','xmax','ymin','ymax','zmin','zmax'}
    b=oldBoundary.(key{1});assert(strcmp(b.type,'characteristic')&&strcmp(b.mode,'riemannFarfield'));
    boundary.(key{1})=struct('type','riemannReservoir','primitive',b.reservoirPrimitive);
end
timer=tic;[rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);elapsed=toc(timer);mapped=mapCutCellWallFields(cut,info,false);
for key={'primitive','wallPressure','wallHeatFluxIntoWall','wallTemperatureInput','wallFlux'}
    assert(isequaln(info.(key{1}),legacyInfo.(key{1})),'NSRiemannIntegration:Scope','同Q给定静止墙的原面结果不应改变：%s',key{1});
end
assert(~isequaln(rate,legacyRate)&&numel(info.riemannReservoirBoundaryAudit.faceID)==19752&&~mapped.eligibleForSurfaceCoverage);
e=s.result.linearLog{1}.steadyDirectionEvidence;candidateAlpha=1/8192;
candidateQ=s.result.initialQ+candidateAlpha*e.directionU.*scale;
a=s.result.lineLog([s.result.lineLog.alpha]==candidateAlpha);assert(isscalar(a)&&~a.accepted&&strcmp(a.errorIdentifier,'RocketAeroHeat:NSDraft:CharacteristicRegime'));
candidateRate=[];candidateInfo=struct;candidateError='';candidatePassed=false;
try
    [candidateRate,candidateInfo]=cutCellNavierStokesResidual(candidateQ,cut,gas,boundary);candidatePassed=true;
catch ex,candidateError=ex.identifier;end
f=cut.openFaces;outside=f.cellPlus==0;w=cut.walls;
denom=max(1,sum(abs(info.faceFlux(outside,:)).*f.area(outside),1)+sum(abs(info.wallFlux).*w.area,1));
balance=max(abs(info.balanceResidual)./denom);assert(balance<1e-9);
report=struct('status','tested-real-20deg-explicit-Riemann-reservoir-same-Q-complete-NS-and-formerly-rejected-candidate', ...
    'sourceFile',sourceFile,'sourceSHA256',rocketStreamingFileSHA256(sourceFile),'oldCompleteNSReplayed',true, ...
    'cellCount',size(Q,1),'sourceFaces',size(cut.inputGeometry.F,1),'riemannReservoirFaces',numel(info.riemannReservoirBoundaryAudit.faceID), ...
    'sourceMaximumResidual',max(abs(legacyRate)./scale,[],'all'),'newBoundaryMaximumResidual',max(abs(rate)./scale,[],'all'), ...
    'sourceRmsMerit',sqrt(mean((legacyRate./scale).^2,'all')),'newBoundaryRmsMerit',sqrt(mean((rate./scale).^2,'all')), ...
    'sameQWallPressureAndHeatKept',true,'spatialConservationRelative',balance,'newBaselineElapsedSeconds',elapsed, ...
    'candidateAlpha',candidateAlpha,'candidateOriginalError',a.errorIdentifier,'candidateNewFullNSPassed',candidatePassed,'candidateNewError',candidateError, ...
    'sameGoverningNSAndInteriorMethod',true,'openBoundaryDiscretizationChanged',true,'sourceHistoryAppended',false, ...
    'separateLegacyReplayNSEvaluations',1,'newBoundaryActualNSEvaluations',2,'newFlowAcceptances',0, ...
    'physicalTimeAdvanced',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'riemann_reservoir_real_evidence.mat'),'Q','gas','boundary','oldBoundary','scale','rate','info','legacyRate','legacyInfo', ...
    'mapped','candidateQ','candidateRate','candidateInfo','report','-v7');
fid=fopen(fullfile(folder,'riemann_reservoir_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
