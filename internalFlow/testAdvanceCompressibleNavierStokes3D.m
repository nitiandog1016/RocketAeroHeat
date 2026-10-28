function report=testAdvanceCompressibleNavierStokes3D(outputFolder)
%TESTADVANCECOMPRESSIBLENAVIERSTOKES3D 小网格真实推进/稳态/全局守恒与回退校核。
% 模块作用：用解析平壁导热比较最终墙热，不把步数上限当收敛；
% 人工强能量移除算子仅检查整步回退，不是NS物理算例或基底场。
sourceFolder=fileparts(mfilename('fullpath'));addpath(sourceFolder,'-begin');
if nargin<1,folder=fullfile(fileparts(sourceFolder),'outputs','v57_validation','solver');
else,folder=char(outputFolder);end
if ~exist(folder,'dir'),mkdir(folder);end
kernelPath=which('compressibleNavierStokesResidual3D');kernelHash=sha256File(kernelPath);
gas=struct('gamma',1.4,'R',287.05,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]);
periodic=struct('type','periodic');wall=@(T)struct('type','wall','temperature',T,'velocity',[0 0 0]);
bc=struct('xmin',periodic,'xmax',periodic,'ymin',periodic,'ymax',periodic,'zmin',periodic,'zmax',periodic);

%% 1. 均匀超声速流在真实RK2时间推进中保持；仅到时不能声称稳态收敛
shape=[3 3 2];grid=struct('spacing',[.05 .02 .05]);P=zeros([shape 5]);
P(:,:,:,1)=1.2;P(:,:,:,2)=1100;P(:,:,:,3)=27;P(:,:,:,4)=-19;P(:,:,:,5)=300;
Q=pack(P,gas);freeBC=bc;freeBC.xmin=struct('type','supersonicInflow','primitive',[1.2 1100 27 -19 300]);
freeBC.xmax=struct('type','supersonicOutflow');
uniform=advanceCompressibleNavierStokes3D(Q,grid,false(shape),gas,freeBC,struct('targetTime',2.5e-5));
assert(isequal(uniform.Q,Q) && uniform.completedRequestedTime && ~uniform.steadyConverged);
assert(strcmp(uniform.stopReason,'target-time-reached'));
assert(uniform.time==2.5e-5 && uniform.acceptedSteps>0);
assert(all(uniform.history.residualNorm==0) && uniform.globalBalanceRelative==0);

%% 2. 固体内部始终NaN；当前墙温是给定，最后墙热由最后态派生
shapeSolid=[4 4 3];solid=false(shapeSolid);solid(2:3,2:3,2)=true;
P=zeros([shapeSolid 5]);P(:,:,:,1)=1;P(:,:,:,5)=300;Q=pack(P,gas);
for component=1:5,v=Q(:,:,:,component);v(solid)=NaN;Q(:,:,:,component)=v;end
solidBC=bc;solidBC.solidWall=wall(300);
static=advanceCompressibleNavierStokes3D(Q,grid,solid,gas,solidBC, ...
    struct('steadyResidualTolerance',1e-8,'steadyRequiredSteps',3,'minimumSteadySteps',3,'maxSteps',20));
assert(static.steadyConverged && strcmp(static.stopReason,'steady-residual-converged'));
f=reshape(static.Q,[],5);assert(all(isnan(f(solid(:),:)),'all'));
assert(all([static.wallFaces.heatFluxIntoWall]==0));
assert(strcmp(static.wallTemperatureValueKind,'prescribed') && ~static.eligibleForSurfaceCoverage);

%% 3. 非均匀闭周期三维真实NS：五量总改变等于RK2净边界积分
shape=[4 3 2];P=zeros([shape 5]);
for i=1:shape(1)
    for j=1:shape(2)
        for z=1:shape(3)
            P(i,j,z,:)=[1+.03*sin(2*pi*i/shape(1)),10+.3*cos(2*pi*j/shape(2)), ...
                .2*sin(2*pi*i/shape(1)),.1*cos(pi*z),300+2*cos(2*pi*j/shape(2))];
        end
    end
end
Q=pack(P,gas);closed=advanceCompressibleNavierStokes3D(Q,grid,false(shape),gas,bc, ...
    struct('targetTime',5e-5,'maxSteps',100));
assert(closed.completedRequestedTime && ~closed.steadyConverged);
assert(closed.globalBalanceRelative<1e-12 && max(abs(closed.accumulatedBoundaryIntegral))==0);
assert(max(abs(closed.finalConservedIntegral-closed.initialConservedIntegral)./ ...
    max(1,abs(closed.initialConservedIntegral)))<1e-12);
assert(~isequal(closed.Q,Q),'非均匀场必须真实推进，不得原样返回假历史。');

%% 4. 两墙导热从扰动状态真正收敛：常物性模型的解析T(y)、墙热和质量
thermalGas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
thermalGrid=struct('spacing',[1 .01 1]);thermalShape=[1 6 1];L=.06;Tlo=300;Thi=330;
P=zeros([thermalShape 5]);initialPressure=315;
for j=1:thermalShape(2)
    y=(j-.5)*.01;T=Tlo+(Thi-Tlo)*y/L+8*sin(pi*y/L);
    P(:,j,:,1)=initialPressure/(thermalGas.R*T);P(:,j,:,5)=T;
end
thermalBC=bc;thermalBC.ymin=wall(Tlo);thermalBC.ymax=wall(Thi);
options=struct('steadyResidualTolerance',2e-6,'steadyRequiredSteps',8,'minimumSteadySteps',8, ...
    'maxSteps',4000,'convectiveCFL',.45,'diffusiveCFL',.45);
timer=tic;thermal=advanceCompressibleNavierStokes3D(pack(P,thermalGas),thermalGrid, ...
    false(thermalShape),thermalGas,thermalBC,options);thermalSeconds=toc(timer);
assert(thermal.steadyConverged,'两墙导热没有达到真实残差停止条件。');
assert(thermal.residualNorm<=options.steadyResidualTolerance);
assert(all(thermal.history.residualNorm(end-7:end)<=options.steadyResidualTolerance));
T=unpackTemperature(thermal.Q,thermalGas);analytical=Tlo+(Thi-Tlo)*((1:6)-.5)/6;
temperatureError=max(abs(T(:)'-analytical));
assert(temperatureError<1e-4,'稳态温度未匹配线性解析解。');
cp=thermalGas.gamma*thermalGas.R/(thermalGas.gamma-1);k=thermalGas.mu*cp/thermalGas.Pr;
heatExact=k*(Thi-Tlo)/L;bottom=arrayfun(@(w)w.normal(2)<0,thermal.wallFaces);
top=arrayfun(@(w)w.normal(2)>0,thermal.wallFaces);
heatError=max([abs([thermal.wallFaces(bottom).heatFluxIntoWall]/heatExact-1), ...
    abs([thermal.wallFaces(top).heatFluxIntoWall]/heatExact+1)]);
assert(heatError<1e-5,'最终墙热未匹配独立Fourier解析量。');
assert(thermal.globalBalanceRelative<1e-10);
assert(all([thermal.wallFaces(bottom).temperatureInput]==Tlo) && ...
    all([thermal.wallFaces(top).temperatureInput]==Thi));
assert(~thermal.physicalSolution && ~thermal.eligibleForSurfaceCoverage);

%% 5. Couette从非稳态迭代：墙机械功必须变为两墙总热，而非只算初始残差
% 这是人工常gamma气体（R=1 J/kg/K、mu=.5 Pa*s），不是实际空气工况。
% 对该有限体积半单元壁梯度，二次温度的中心值有明确O(dy^2)偏移。
% 单独核对解析点值网格误差与修正后停止误差，不放宽残差门限。
couetteGas=thermalGas;couetteGrid=struct('spacing',[1 .0075 1]);couetteShape=[1 8 1];
couetteL=.06;wallT=300;upperU=3;P=zeros([couetteShape 5]);
P(:,:,:,1)=1;P(:,:,:,5)=wallT; % 初始u=0、T=Tw，与解析u/T均不同。
couetteBC=bc;couetteBC.ymin=wall(wallT);
couetteBC.ymax=struct('type','wall','temperature',wallT,'velocity',[upperU 0 0]);
timer=tic;couette=advanceCompressibleNavierStokes3D(pack(P,couetteGas),couetteGrid, ...
    false(couetteShape),couetteGas,couetteBC,options);couetteSeconds=toc(timer);
assert(couette.steadyConverged,'Couette黏性发热没有达到真实残差停止条件。');
assert(all(couette.history.residualNorm(end-7:end)<=options.steadyResidualTolerance));
eta=((1:8)-.5)/8;analyticU=upperU*eta;
raw=couette.Q;rho=raw(:,:,:,1);velocity=raw(:,:,:,2)./rho;
couetteVelocityError=max(abs(velocity(:)'-analyticU));
assert(couetteVelocityError<2e-6,'最终Couette速度未匹配独立解析线性剖面。');
cp=couetteGas.gamma*couetteGas.R/(couetteGas.gamma-1);couetteK=couetteGas.mu*cp/couetteGas.Pr;
coefficient=couetteGas.mu*upperU^2/(2*couetteK);
analyticT=wallT+coefficient*eta.*(1-eta);
discreteOffset=coefficient/(4*couetteShape(2)^2);
couetteT=unpackTemperature(couette.Q,couetteGas);
couettePointwiseTemperatureError=max(abs(couetteT(:)'-analyticT));
couetteTemperatureStoppingError=max(abs(couetteT(:)'-(analyticT+discreteOffset)));
assert(couetteTemperatureStoppingError<3e-5, ...
    '温度未匹配二阶网格误差模型；不能用粗温度门限隐藏发热符号错误。');
assert(abs(couettePointwiseTemperatureError-discreteOffset)<3e-5);
bottom=arrayfun(@(w)w.normal(2)<0,couette.wallFaces);top=arrayfun(@(w)w.normal(2)>0,couette.wallFaces);
exactShear=couetteGas.mu*upperU/couetteL;exactWallHeat=couetteGas.mu*upperU^2/(2*couetteL);
couetteHeatError=max(abs([couette.wallFaces.heatFluxIntoWall]/exactWallHeat-1));
assert(couetteHeatError<1e-5,'两墙热流不符合黏性机械功的独立解析值。');
assert(max(abs(arrayfun(@(w)w.traction(1)+exactShear,couette.wallFaces(top))))/exactShear<1e-5);
assert(max(abs(arrayfun(@(w)w.traction(1)-exactShear,couette.wallFaces(bottom))))/exactShear<1e-5);
totalWallHeat=sum([couette.wallFaces.heatFluxIntoWall].*[couette.wallFaces.area]);
movingWallPower=-sum(arrayfun(@(w)dot(w.velocityInput,w.traction)*w.area,couette.wallFaces(top)));
couetteSteadyPowerBalance=abs(totalWallHeat-movingWallPower)/movingWallPower;
assert(couetteSteadyPowerBalance<1e-5);
assert(abs(movingWallPower/(couetteGas.mu*upperU^2/couetteL)-1)<1e-5);
assert(couette.globalBalanceRelative<1e-10);
assert(all([couette.wallFaces.temperatureInput]==wallT));
assert(strcmp(couette.wallTemperatureValueKind,'prescribed') && ~couette.eligibleForSurfaceCoverage);

%% 6. 可核对的整步回退：人工守恒能量移除算子，不充当任何物理模型
shape=[2 2 1];P=zeros([shape 5]);P(:,:,:,1)=1;P(:,:,:,5)=300;Q=pack(P,gas);
a=sqrt(gas.gamma*gas.R*300);dtReference=1/(sum(a./grid.spacing)/.45 + ...
    2*max((4/3)*gas.mu,gas.mu*(gas.gamma*gas.R/(gas.gamma-1))/gas.Pr/(gas.R/(gas.gamma-1))) ...
    *sum(1./grid.spacing.^2)/.45);
sink=gas.R/(gas.gamma-1)*300/(.1*dtReference);
controlled=@(q,g,s,air,b)energySinkFixture(q,g,s,air,b,sink);
retry=advanceCompressibleNavierStokes3D(Q,grid,false(shape),gas,bc, ...
    struct('targetTime',2*dtReference,'maxSteps',1,'maxRetries',12,'residualFunction',controlled));
assert(retry.rejectedAttempts>0 && retry.acceptedSteps==1 && ~retry.converged);
assert(strcmp(retry.stopReason,'maximum-steps-reached') && ~retry.completedRequestedTime);
assert(retry.history.acceptedDt(end)<retry.history.proposedDt(end));
assert(abs(retry.history.acceptedDt(end)/retry.history.proposedDt(end)-2^-retry.history.retries(end))<1e-12);
assert(retry.globalBalanceRelative<1e-12);
assert(all(unpackTemperature(retry.Q,gas)>=80,'all') && retry.positivityClippingCount==0);
reject(@()advanceCompressibleNavierStokes3D(Q,grid,false(shape),gas,bc, ...
    struct('targetTime',1e-5,'convectiveCFL',1.1)),'RocketAeroHeat:NSAdvance:CFL');
reject(@()advanceCompressibleNavierStokes3D(Q,grid,false(shape),gas,bc, ...
    struct('targetTime',2*dtReference,'maxSteps',1,'maxRetries',0,'residualFunction',controlled)), ...
    'RocketAeroHeat:NSAdvance:RetryLimit');
assert(all(retry.globalBalanceRelative<=retry.options.conservationTolerance));
assert(strcmp(kernelHash,sha256File(kernelPath)),'测试期间空间内核源文件改变，证据不能归于同一版本。');
report=struct('date','2026-10-14','tested',true,'experimentalValidation',false, ...
    'eligibleForSurfaceCoverage',false,'uniformSteps',uniform.acceptedSteps, ...
    'closedPeriodicSteps',closed.acceptedSteps,'closedPeriodicGlobalBalance',closed.globalBalanceRelative, ...
    'thermalSteadySteps',thermal.acceptedSteps,'thermalSteadyTime',thermal.time, ...
    'thermalSeconds',thermalSeconds,'thermalResidual',thermal.residualNorm, ...
    'thermalTemperatureErrorK',temperatureError,'thermalWallHeatRelativeError',heatError, ...
    'thermalWallHeatExactWm2',heatExact,'thermalGlobalBalance',thermal.globalBalanceRelative, ...
    'couetteTestGas',couetteGas,'couetteSteadySteps',couette.acceptedSteps,'couetteSteadyTime',couette.time, ...
    'couetteSeconds',couetteSeconds,'couetteResidual',couette.residualNorm, ...
    'couetteVelocityErrorMps',couetteVelocityError,'couettePointwiseTemperatureGridErrorK',couettePointwiseTemperatureError, ...
    'couetteAnalyticDiscreteOffsetK',discreteOffset,'couetteTemperatureStoppingErrorK',couetteTemperatureStoppingError, ...
    'couetteWallHeatRelativeError',couetteHeatError,'couetteWallHeatExactWm2',exactWallHeat, ...
    'couetteTotalWallHeatW',totalWallHeat,'couetteMovingWallPowerW',movingWallPower, ...
    'couetteSteadyPowerRelativeBalance',couetteSteadyPowerBalance,'couetteGlobalBalance',couette.globalBalanceRelative, ...
    'controlledRetries',retry.rejectedAttempts,'controlledAcceptedDt',retry.history.acceptedDt(end), ...
    'controlledProposedDt',retry.history.proposedDt(end),'controlledGlobalBalance',retry.globalBalanceRelative, ...
    'residualKernelPath',kernelPath,'residualKernelSha256',kernelHash, ...
    'scope','小网格RK2、平壁解析导热及Couette真实稳态；人工常gamma气体非空气工况；人工回退fixture不属物理验证；无真实底区/箭体/实验/全网格边界收敛');
save(fullfile(folder,'advance_test_evidence.mat'),'report','thermal','couette','uniform','closed','retry','static');
fid=fopen(fullfile(folder,'advance_test_summary.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
fprintf('NS RK2 advance passed: uniform and periodic conservation, solid NaN, true two-wall thermal steady (%d steps), viscous Couette steady (%d steps, %.9g residual), controlled whole-step retry=%d; no full-area eligibility.\n', ...
    thermal.acceptedSteps,couette.acceptedSteps,couette.residualNorm,retry.rejectedAttempts);disp(report);
end
function Q=pack(P,gas)
rho=P(:,:,:,1);u=P(:,:,:,2);v=P(:,:,:,3);w=P(:,:,:,4);T=P(:,:,:,5);Q=P;
Q(:,:,:,2)=rho.*u;Q(:,:,:,3)=rho.*v;Q(:,:,:,4)=rho.*w;
Q(:,:,:,5)=rho.*(gas.R/(gas.gamma-1)*T+.5*(u.^2+v.^2+w.^2));
end
function T=unpackTemperature(Q,gas)
rho=Q(:,:,:,1);kinetic=.5*(Q(:,:,:,2).^2+Q(:,:,:,3).^2+Q(:,:,:,4).^2)./rho;
T=(Q(:,:,:,5)-kinetic)./(rho*gas.R/(gas.gamma-1));
end
function [R,info]=energySinkFixture(Q,grid,solid,gas,boundary,sink)
% 模块作用：人工一致体积残差/净边界积分，只测试RK拒绝与守恒累计。
[R,info]=compressibleNavierStokesResidual3D(Q,grid,solid,gas,boundary);
v=R(:,:,:,5);v(~solid)=v(~solid)-sink;R(:,:,:,5)=v;
info.netBoundaryFlux(5)=info.netBoundaryFlux(5)+sink*nnz(~solid)*prod(grid.spacing);
info.method='controlled energy-sink test fixture; not a physical NS boundary';
end
function reject(action,identifier)
caught=false;try,action();catch exception,assert(strcmp(exception.identifier,identifier),exception.message);caught=true;end
assert(caught,['Expected failure: ' identifier]);
end
function value=sha256File(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(typecast(bytes,'int8'));
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
