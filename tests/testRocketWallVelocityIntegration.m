function report=testRocketWallVelocityIntegration(folder)
%TESTROCKETWALLVELOCITYINTEGRATION 在真实20度原墙同Q检查静止壁速度约束的完整NS作用域。
% 模块一：缺省/显式legacy完整rate/info精确复现已交付v82末态。
% 模块二：只改变solidWall的Newton应力，质量/能量/墙压/Fourier和全部内部通量保持。
% 模块三：逐墙闭式牵引、机械耗散/静止墙熵配对及完整空间五量账，不作物理资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'WallVelocityIntegration:Overwrite','不得覆盖真实同Q记录。');mkdir(folder);
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;
s=load(fullfile(root,'outputs','v82_validation','exact_krylov01','steady_krylov_evidence.mat'),'result','gas','boundary');
Q=s.result.Q;gas=s.gas;boundary=s.boundary;scale=s.result.residualScale;
[legacyRate,legacyInfo]=cutCellNavierStokesResidual(Q,cut,gas,boundary);
assert(isequaln(legacyRate,s.result.rate)&&isequaln(legacyInfo,s.result.info), ...
    'WallVelocityIntegration:Default','缺省完整空间算子没有精确复现v82源状态。');clear s;
boundary.solidWall.velocityGradient='legacy';[r,i]=cutCellNavierStokesResidual(Q,cut,gas,boundary);
assert(isequaln(r,legacyRate)&&isequaln(i,legacyInfo));clear r i;
fprintf('WALL_VELOCITY_DEFAULT_COMPLETE_V82_REPLAY_PASSED\n');
bad=boundary;bad.velocityGradient='dirichlet-tangent';rejected=false;
try,cutCellNavierStokesResidual(Q,cut,gas,bad);catch ex,rejected=strcmp(ex.identifier,'CutNS:WallVelocityBoundary');end;assert(rejected);clear bad;
boundary.solidWall.velocityGradient='dirichlet-tangent';timer=tic;
[rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);elapsed=toc(timer);
for name={'primitive','gradient','faceFlux','convectiveFaceFlux','viscousFaceFlux','wallPressure', ...
        'wallHeatFluxIntoWall','wallTemperatureInput','wallMechanicalPowerIntoFluid','wallThermalGradientAudit', ...
        'openBoundaryFlux','outerWallBoundaryFlux','viscousWorkAudit','exactEulerFluxAudit'}
    key=name{1};assert(isequaln(info.(key),legacyInfo.(key)),'WallVelocityIntegration:Scope','不该改变的量发生变化：%s',key);
end
assert(isequaln(rate(:,[1 5]),legacyRate(:,[1 5]))&&isequaln(info.wallFlux(:,[1 5]),legacyInfo.wallFlux(:,[1 5]))&& ...
    isequaln(info.netBoundaryFlux([1 5]),legacyInfo.netBoundaryFlux([1 5]))&&any(info.wallFlux(:,2:4)~=legacyInfo.wallFlux(:,2:4),'all'));
w=cut.walls;owner=w.ownerCell;n=w.normal;P=info.primitive;u=P(owner,2:4);T=P(owner,5);
dn=sum((w.centroid-cut.cells.centroid(owner,:)).*n,2);traction=info.wallPressure.*n-info.wallFlux(:,2:4);
oldTraction=legacyInfo.wallPressure.*n-legacyInfo.wallFlux(:,2:4);
expected=-gas.mu./dn.*(u+sum(u.*n,2).*n/3);
tractionDifference=max(abs(traction-expected),[],'all')/max(1,max(abs(expected),[],'all'));assert(tractionDifference<3e-12);
mechanical=-sum(u.*traction,2);oldMechanical=-sum(u.*oldTraction,2);diffusionEnergy=-info.wallFlux(:,5);
entropy=w.area.*(mechanical./T+(1./T-1./info.wallTemperatureInput).*diffusionEnergy);
oldEntropy=w.area.*(oldMechanical./T+(1./T-1./info.wallTemperatureInput).*diffusionEnergy);
assert(all(mechanical>=-128*eps*max(1,max(abs(mechanical)))))
assert(all(entropy>=-128*eps*max(1,max(abs(entropy)))))
f=cut.openFaces;outside=f.cellPlus==0;denom=max(1,sum(abs(info.faceFlux(outside,:)).*f.area(outside),1)+sum(abs(info.wallFlux).*w.area,1));
spatial=max(abs(info.balanceResidual)./denom);assert(spatial<1e-9);
cv=gas.R/(gas.gamma-1);cellID=14281;ur=P(cellID,2:4);rho=P(cellID,1);tc=P(cellID,5);
slope=@(r)(r(5)-sum(ur.*r(2:4))+(.5*sum(ur.^2)-cv*tc)*r(1))/(rho*cv);
report=struct('status','passed-real-20deg-static-wall-velocity-complete-NS-scope-and-legacy-identity', ...
    'cellCount',size(Q,1),'wallCount',numel(w.area),'defaultCompleteInfoExactlyReplayed',true, ...
    'massEnergyAllInteriorPressureAndFourierKept',true,'tractionRelativeDifference',tractionDifference, ...
    'oldMechanicalCoolingWalls',nnz(oldMechanical<0),'newMechanicalCoolingWalls',nnz(mechanical<0), ...
    'oldNegativeWallEntropyFaces',nnz(oldEntropy<0),'newNegativeWallEntropyFaces',nnz(entropy<0), ...
    'initialNewMaximumResidual',max(abs(rate)./scale,[],'all'),'initialNewRmsMerit',sqrt(mean((rate./scale).^2,'all')), ...
    'oldColdCellTemperatureRateKs',slope(legacyRate(cellID,:)),'newColdCellTemperatureRateKs',slope(rate(cellID,:)), ...
    'spatialConservationRelative',spatial,'elapsedSeconds',elapsed,'wholeEntropyBalanceProved',false, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'wall_velocity_integration_evidence.mat'),'report','Q','rate','info','gas','boundary','scale','legacyRate','legacyInfo','-v7');
fid=fopen(fullfile(folder,'wall_velocity_integration_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
