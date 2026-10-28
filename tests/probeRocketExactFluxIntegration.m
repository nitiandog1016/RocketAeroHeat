function report=probeRocketExactFluxIntegration(folder)
%PROBEROCKETEXACTFLUXINTEGRATION 在同一真实六次端点测试新显式全内面通量及默认完整复现。
% 模块一：旧默认/显式hllc完整rate/info精确保持，不重绑v81旧源码。
% 模块二：批量849336原内部面精确Riemann，保持黏性/墙热/外边界与全空间五量账。
% 模块三：失败保留真实标识，不将截断未知区/填补/降低温域视为通过。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CutExactIntegration:Overwrite','不得覆盖真实试验。');mkdir(folder);
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
s=load(fullfile(root,'outputs','v81_validation','arithmetic02','steady_krylov_evidence.mat'),'result','gas','boundary');
Q=s.result.Q;gas=s.gas;boundary=s.boundary;cut=d.domain;
[r0,i0]=cutCellNavierStokesResidual(Q,cut,gas,boundary);assert(isequaln(r0,s.result.rate)&&isequaln(i0,s.result.info));
boundary.interiorEulerFlux='hllc';[r1,i1]=cutCellNavierStokesResidual(Q,cut,gas,boundary);assert(isequaln(r0,r1)&&isequaln(i0,i1));clear r1 i1;
fprintf('EXACT_REAL_COMPLETE_LEGACY_INFO_PASSED\n');
bad={struct('interiorEulerFlux','wrong'),struct('interiorEulerFlux',17),struct('solidWall',struct('interiorEulerFlux','exact-godunov'))};
for j=1:numel(bad),rejected=false;try,cutCellInteriorEulerPolicy(bad{j});catch ex,rejected=strcmp(ex.identifier,'CutNS:InteriorFluxPolicy');end;assert(rejected);end
boundary.interiorEulerFlux='exact-godunov';timer=tic;
try
    [rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);
catch ex
    report=struct('status','rejected-real-six-correction-state-explicit-exact-flux','identifier',ex.identifier,'message',ex.message, ...
        'defaultCompleteInfoExactlyReplayed',true,'elapsedSeconds',toc(timer),'physicalSolution',false,'newProductionSolvedAreaM2',0);
    write(fullfile(folder,'integration_failure.json'),jsonencode(report));rethrow(ex);
end
elapsed=toc(timer);assert(isequaln(info.primitive,i0.primitive)&&isequaln(info.gradient,i0.gradient)&& ...
    isequaln(info.viscousFaceFlux,i0.viscousFaceFlux)&&isequaln(info.wallFlux,i0.wallFlux)&& ...
    isequaln(info.wallHeatFluxIntoWall,i0.wallHeatFluxIntoWall)&&isequaln(info.wallThermalGradientAudit,i0.wallThermalGradientAudit)&& ...
    isequaln(info.netBoundaryFlux,i0.netBoundaryFlux)&&isequaln(info.viscousWorkAudit,i0.viscousWorkAudit));
f=cut.openFaces;outside=f.cellMinus==0|f.cellPlus==0;
assert(isequaln(info.convectiveFaceFlux(outside,:),i0.convectiveFaceFlux(outside,:)));
scale=s.result.residualScale;residual=max(abs(rate)./scale,[],'all');rms=sqrt(mean((rate./scale).^2,'all'));
denom=max(1,sum(abs(info.faceFlux(outside,:)).*f.area(outside),1)+sum(abs(info.wallFlux).*cut.walls.area,1));
balance=max(abs(info.balanceResidual)./denom);assert(balance<1e-9);
a=info.exactEulerFluxAudit;key=find(a.faceID==70995);assert(numel(key)==1&&a.pressureStarPa(key)>100);
report=struct('status','passed-real-six-correction-state-full-explicit-exact-Euler-integration-and-legacy-info', ...
    'cellCount',numel(cut.cells.volume),'sourceFaces',19752,'interiorFaces',numel(a.faceID), ...
    'maximumResidual',residual,'rmsMerit',rms,'spatialConservationRelative',balance,'elapsedSeconds',elapsed, ...
    'actualFacePressureStarPa',a.pressureStarPa(key),'actualFaceSampleTemperatureK',a.samplePrimitive(key,5), ...
    'maxPressureRootIterations',max(a.pressureRootIterations),'defaultCompleteInfoExactlyReplayed',true, ...
    'viscousWallAndBoundaryExactlyKept',true,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'exact_integration_evidence.mat'),'report','Q','rate','info','gas','boundary','-v7');
write(fullfile(folder,'integration_summary.json'),jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
