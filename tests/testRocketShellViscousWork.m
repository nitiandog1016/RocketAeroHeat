function report=testRocketShellViscousWork(folder)
%TESTROCKETSHELLVISCOUSWORK 真实20度原状态的完整旧路径及只变内面功的受控校核。
% 模块一：缺省与显式legacy完整rate/info重放v80，不重绑旧源码。
% 模块二：新策略只改共享内部能量通量；牵引、热梯度、壁热和全部边界保持。
% 模块三：独立由旧WLS及真实几何重构旧功速度和半功恒等式；不证明一般正性。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'ShellWorkTest:Overwrite','不得覆盖已有检查。');mkdir(folder);
base=fullfile(root,'outputs','v80_validation');if ~exist(fullfile(base,'advance01','advance_evidence.mat'),'file'),base=fullfile(root,'validation','v80');end
d=load(fullfile(base,'domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;
s=load(fullfile(base,'advance01','advance_evidence.mat'),'result','gas','boundary');Q=s.result.Q;gas=s.gas;boundary=s.boundary;
[r0,i0]=cutCellNavierStokesResidual(Q,cut,gas,boundary);assert(isequaln(r0,s.result.residual)&&isequaln(i0,s.result.finalSpatialInfo));
boundary.viscousWorkVelocity='legacy-wls';[r1,i1]=cutCellNavierStokesResidual(Q,cut,gas,boundary);assert(isequaln(r0,r1)&&isequaln(i0,i1));clear r1 i1;
bad={struct('viscousWorkVelocity','bad'),struct('viscousWorkVelocity',17),struct('viscousWorkVelocity',[]), ...
    struct('solidWall',struct('viscousWorkVelocity','arithmetic')),struct('ymax',struct('viscousWorkVelocity','arithmetic'))};
for j=1:numel(bad),rejected=false;try,cutCellViscousWorkPolicy(bad{j});catch ex,rejected=strcmp(ex.identifier,'CutNS:ViscousWorkBoundary');end;assert(rejected);end
boundary.viscousWorkVelocity='arithmetic';[rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);a=info.viscousWorkAudit;
assert(isequaln(rate(:,1:4),r0(:,1:4))&&isequaln(info.primitive,i0.primitive)&&isequaln(info.gradient,i0.gradient)&& ...
    isequaln(info.viscousFaceFlux(:,1:4),i0.viscousFaceFlux(:,1:4))&&isequaln(info.convectiveFaceFlux,i0.convectiveFaceFlux)&& ...
    isequaln(info.wallFlux,i0.wallFlux)&&isequaln(info.wallHeatFluxIntoWall,i0.wallHeatFluxIntoWall)&& ...
    isequaln(info.wallThermalGradientAudit,i0.wallThermalGradientAudit)&&isequaln(info.netBoundaryFlux,i0.netBoundaryFlux));
f=cut.openFaces;ids=a.faceID;m=f.cellMinus(ids);p=f.cellPlus(ids);n=size(Q,1);P=info.primitive;
old=zeros(numel(ids),3);
for j=1:numel(ids)
    face=ids(j);left=P(m(j),2:4)+(i0.gradient(1:3,:,m(j))*(f.centroid(face,:)-cut.cells.centroid(m(j),:))')';
    right=P(p(j),2:4)+(i0.gradient(1:3,:,p(j))*(f.centroid(face,:)-cut.cells.centroid(p(j),:))')';old(j,:)=.5*(left+right);
end
assert(isequaln(old,a.oldVelocityXYZ)&&isequaln(.5*(P(m,2:4)+P(p,2:4)),a.velocityXYZ));
delta=sum((a.velocityXYZ-old).*a.tractionXYZ,2);actual=info.viscousFaceFlux(ids,5)-i0.viscousFaceFlux(ids,5);
workError=max(abs(actual-delta))/max(1,max(abs(delta)));assert(workError<1e-12);
expected=(accumarray(m,f.area(ids).*delta,[n 1])-accumarray(p,f.area(ids).*delta,[n 1]))./cut.cells.volume;
rateError=max(abs((rate(:,5)-r0(:,5))-expected))/max(1,max(abs(expected)));assert(rateError<1e-11);
pairError=max(abs([a.ownerHalfPairDifference;a.neighbourHalfPairDifference]))/max(1,max(abs(a.pairMomentumToInternalEnergyPerArea)));assert(pairError<1e-12);
boundaryFaces=setdiff((1:numel(f.area))',ids);assert(isequaln(info.viscousFaceFlux(boundaryFaces,:),i0.viscousFaceFlux(boundaryFaces,:)));
cv=gas.R/(gas.gamma-1);u=P(:,2:4);rho=P(:,1);
oldTemperatureRate=(r0(:,5)-sum(u.*r0(:,2:4),2)+(.5*sum(u.^2,2)-cv*P(:,5)).*r0(:,1))./(rho*cv);
newTemperatureRate=(rate(:,5)-sum(u.*rate(:,2:4),2)+(.5*sum(u.^2,2)-cv*P(:,5)).*rate(:,1))./(rho*cv);
report=struct('status','passed-actual-20deg-complete-legacy-replay-internal-work-only-and-pair-identity', ...
    'cellCount',n,'sourceFaces',19752,'internalFaces',numel(ids),'legacyCompleteInfoExactlyReplayed',true, ...
    'unchangedMassMomentumStressHeatWallAndBoundary',true,'illegalPolicyRejections',numel(bad), ...
    'workFluxRelativeDifference',workError,'energyResidualChangeRelativeDifference',rateError,'halfPairRelativeDifference',pairError, ...
    'maximumOldWorkSpeed',max(vecnorm(a.oldVelocityXYZ,2,2)),'maximumNewWorkSpeed',max(vecnorm(a.velocityXYZ,2,2)), ...
    'limitingCV',16340,'oldLimitingTemperatureRateKs',oldTemperatureRate(16340),'newLimitingTemperatureRateKs',newTemperatureRate(16340), ...
    'negativePairCount',a.negativePairCount,'generalSkewGridPositivityProved',false,'physicalSolution',false, ...
    'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'work_velocity_evidence.mat'),'report','rate','info','gas','boundary','-v7');
fid=fopen(fullfile(folder,'work_velocity_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
