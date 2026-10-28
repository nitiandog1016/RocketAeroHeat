function report=testCompressibleNavierStokesResidual3D(outputFolder)
%TESTCOMPRESSIBLENAVIERSTOKESRESIDUAL3D 正式独立库空间算子的解析及守恒检查。
% 模块作用：不用迭代图伪造解，比较自由流、实际黏性Couette功、Fourier
% 导热、接触面/激波通量及独立流体边界积分。这里不验证真实底区精度。
folder=fileparts(mfilename('fullpath'));addpath(folder,'-begin');
% 模块作用：证据单独写入输出目录，避免MAT/JSON混入可替换源码。
if nargin<1 || isempty(outputFolder)
    outputFolder=fullfile(fileparts(folder),'outputs','internal_flow_kernel_test');
end
if ~isfolder(outputFolder),mkdir(outputFolder);end
gas=struct('gamma',1.4,'R',287.05,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]);
gas.k=gas.mu*(gas.gamma*gas.R/(gas.gamma-1))/gas.Pr;
grid=struct('spacing',[.025 .005 .02]);shape=[4 8 3];solid=false(shape);
periodic=struct('type','periodic');wall=@(T,u)struct('type','wall','temperature',T,'velocity',u);
bc=struct('xmin',periodic,'xmax',periodic,'ymin',periodic,'ymax',periodic,'zmin',periodic,'zmax',periodic);

%% 1. 向内/向外超声速开放边界：均匀流应不变且真实质量/能量穿流相等
P=zeros([shape 5]);P(:,:,:,1)=1.2;P(:,:,:,2)=1100;P(:,:,:,3)=27;P(:,:,:,4)=-19;P(:,:,:,5)=300;
Q=pack(P,gas);freeBC=bc;freeBC.xmin=struct('type','supersonicInflow','primitive',[1.2 1100 27 -19 300]);
freeBC.xmax=struct('type','supersonicOutflow');
[r,f]=compressibleNavierStokesResidual3D(Q,grid,solid,gas,freeBC);
freeError=max(abs(r(:)));assert(freeError<1e-7);
assert(norm(f.balanceResidual,inf)<1e-7);
expectedMass=1.2*1100*(shape(2)*grid.spacing(2))*(shape(3)*grid.spacing(3));
F=f.faceFlux{1};fluxIn=squeeze(F(1,1,1,:))';fluxOut=squeeze(F(end,1,1,:))';
assert(relative(fluxIn,fluxOut)<1e-13);
assert(abs(fluxIn(1)*grid.spacing(2)*grid.spacing(3)*shape(2)*shape(3)/expectedMass-1)<1e-13);

%% 2. 内部solid的静止等温壁；固体保守量与输出必须NaN
s=false([6 6 4]);s(3:4,3:4,2:3)=true;P=zeros([size(s) 5]);P(:,:,:,1)=1;P(:,:,:,5)=300;
Q=pack(P,gas);for q=1:5,t=Q(:,:,:,q);t(s)=NaN;Q(:,:,:,q)=t;end
staticBC=bc;staticBC.solidWall=wall(300,[0 0 0]);
[r,w]=compressibleNavierStokesResidual3D(Q,grid,s,gas,staticBC);
rf=reshape(r,[],5);assert(all(isnan(rf(s(:),:)),'all'));
assert(max(abs(rf(~s(:),:)),[],'all')<1e-7);
assert(~isempty(w.wallFaces) && all([w.wallFaces.heatFluxIntoWall]==0));
assert(norm(w.volumeIntegratedResidual,inf)<1e-9);
assert(norm(w.netBoundaryFlux,inf)<1e-9);

%% 3. Couette：壁剪力与上墙机械功为独立解析量，不把等温初场当稳态热解
P=zeros([shape 5]);P(:,:,:,1)=1;P(:,:,:,5)=300;
Utop=40;L=shape(2)*grid.spacing(2);A=Utop/L;
for j=1:shape(2),P(:,j,:,2)=A*(j-.5)*grid.spacing(2);end
couetteBC=bc;couetteBC.ymin=wall(300,[0 0 0]);couetteBC.ymax=wall(300,[Utop 0 0]);
[r,c]=compressibleNavierStokesResidual3D(pack(P,gas),grid,solid,gas,couetteBC);
traction=gas.mu*A;wallArea=(shape(1)*grid.spacing(1))*(shape(3)*grid.spacing(3));
top=arrayfun(@(w)w.normal(2)>0,c.wallFaces);bottom=arrayfun(@(w)w.normal(2)<0,c.wallFaces);
assert(max(abs(arrayfun(@(w)w.traction(1)+traction,c.wallFaces(top))))<1e-11);
assert(max(abs(arrayfun(@(w)w.traction(1)-traction,c.wallFaces(bottom))))<1e-11);
assert(max(abs([c.wallFaces.heatFluxIntoWall]))<1e-10);
energy=r(:,:,:,5);expectedHeating=gas.mu*A^2;
couetteError=max(abs(energy(:)-expectedHeating))/expectedHeating;
assert(couetteError<1e-8);
assert(abs(c.volumeIntegratedResidual(5)/(traction*Utop*wallArea)-1)<1e-8);
assert(norm(c.balanceResidual,inf)<1e-8);

%% 4. 线性温度：常压静止导热的局部壁热流及全域零净热率
P=zeros([shape 5]);Tlo=300;Thi=330;B=(Thi-Tlo)/L;pressure=60000;
for j=1:shape(2)
    T=Tlo+B*(j-.5)*grid.spacing(2);P(:,j,:,5)=T;P(:,j,:,1)=pressure/(gas.R*T);
end
thermalBC=bc;thermalBC.ymin=wall(Tlo,[0 0 0]);thermalBC.ymax=wall(Thi,[0 0 0]);
[r,t]=compressibleNavierStokesResidual3D(pack(P,gas),grid,solid,gas,thermalBC);
top=arrayfun(@(w)w.normal(2)>0,t.wallFaces);bottom=arrayfun(@(w)w.normal(2)<0,t.wallFaces);
expectedFlux=gas.k*B;
thermalError=max(abs([t.wallFaces(bottom).heatFluxIntoWall]/expectedFlux-1));
assert(thermalError<1e-12);
assert(max(abs([t.wallFaces(top).heatFluxIntoWall]/expectedFlux+1))<1e-12);
assert(abs(t.netBoundaryFlux(5))<1e-9);
assert(max(abs(r(:,:,:,5)),[],'all')<1e-5);

%% 5. 静止接触界面：HLLC不能制造质量穿流；非零热流来自真实Fourier项
P=zeros([shape 5]);P(:,:,:,1)=1;P(:,:,:,5)=300;P(3:4,:,:,1)=2;P(3:4,:,:,5)=150;
[~,contact]=compressibleNavierStokesResidual3D(pack(P,gas),grid,solid,gas,bc);
F=squeeze(contact.faceFlux{1}(3,2,2,:))';
assert(abs(F(1))<1e-12 && abs(F(2)/(gas.R*300)-1)<1e-13);
assert(abs(F(5)/(gas.k*150/grid.spacing(1))-1)<1e-12);

%% 6. 正激波独立Rankine-Hugoniot：远离跃迁面的相同质量/动量/焓通量
shapeShock=[12 3 2];P=zeros([shapeShock 5]);M=2;T1=300;rho1=1;g=gas.gamma;
u1=M*sqrt(g*gas.R*T1);densityRatio=(g+1)*M^2/((g-1)*M^2+2);
pressureRatio=(2*g*M^2-(g-1))/(g+1);T2=T1*pressureRatio/densityRatio;
P(:,:,:,1)=rho1;P(:,:,:,2)=u1;P(:,:,:,5)=T1;
P(7:end,:,:,1)=rho1*densityRatio;P(7:end,:,:,2)=u1/densityRatio;P(7:end,:,:,5)=T2;
[~,shock]=compressibleNavierStokesResidual3D(pack(P,gas),grid,false(shapeShock),gas,bc);
before=squeeze(shock.faceFlux{1}(3,2,1,:))';after=squeeze(shock.faceFlux{1}(10,2,1,:))';
shockConservationError=relative(before,after);assert(shockConservationError<1e-12);
assert(all(isfinite(shock.faceFlux{1}(:))));
% 没有时间迭代，不能称“激波位置/厚度或NS激波剖面已收敛”。

%% 7. 任意内部流场：独立逐流体单元边界枚举与体积分核对五量守恒
shapeBox=[5 6 4];s=false(shapeBox);s(3,3:4,2:3)=true;
P=zeros([shapeBox 5]);
for i=1:shapeBox(1)
    for j=1:shapeBox(2)
        for z=1:shapeBox(3)
            P(i,j,z,:)=[1+.02*i,.9*j,-.7*i,.3*z,300+2*i+3*j+z];
        end
    end
end
Q=pack(P,gas);for q=1:5,v=Q(:,:,:,q);v(s)=NaN;Q(:,:,:,q)=v;end
boxBC=struct('xmin',wall(300,[0 0 0]),'xmax',wall(300,[0 0 0]), ...
    'ymin',wall(300,[0 0 0]),'ymax',wall(300,[0 0 0]), ...
    'zmin',wall(300,[0 0 0]),'zmax',wall(300,[0 0 0]),'solidWall',wall(300,[0 0 0]));
[r,b]=compressibleNavierStokesResidual3D(Q,grid,s,gas,boxBC);
expectedBoundary=independentWallIntegral(P,s,grid.spacing,gas,300);
rf=reshape(r,[],5);integral=sum(rf(~s(:),:),1)*prod(grid.spacing);
closedConservationError=relative(integral,-expectedBoundary);
assert(closedConservationError<1e-10);
assert(relative(b.netBoundaryFlux,expectedBoundary)<1e-12);
assert(norm(b.balanceResidual,inf)<1e-8);

%% 8. SI尺度协变：固定速度温度，所有长度翻倍使Couette加热体密度降四倍
P=zeros([shape 5]);P(:,:,:,1)=1;P(:,:,:,5)=300;
for j=1:shape(2),P(:,j,:,2)=Utop*(j-.5)/shape(2);end
largeGrid=grid;largeGrid.spacing=2*grid.spacing;
[r2,c2]=compressibleNavierStokesResidual3D(pack(P,gas),largeGrid,solid,gas,couetteBC);
assert(max(abs(r2(:,:,:,5)*4-expectedHeating),[],'all')/expectedHeating<1e-8);
assert(abs(c2.volumeIntegratedResidual(5)/c.volumeIntegratedResidual(5)-2)<1e-8);
assert(isequal(c.residualUnits,{'kg/(m^3*s)','N/m^3','N/m^3','N/m^3','W/m^3'}));

%% 9. 明确拒绝非法能量、热边界、开放回流；没有裁温或密度/能量填值
bad=pack(P,gas);bad(1,1,1,5)=-1;
reject(@()compressibleNavierStokesResidual3D(bad,grid,solid,gas,couetteBC),'RocketAeroHeat:NSDraft:InternalEnergy');
badBC=couetteBC;badBC.ymin.temperature=1;
reject(@()compressibleNavierStokesResidual3D(pack(P,gas),grid,solid,gas,badBC),'RocketAeroHeat:NSDraft:Temperature');
badBC=bc;badBC.xmin=wall(300,[0 0 0]);badBC.xmax=struct('type','supersonicOutflow');
reject(@()compressibleNavierStokesResidual3D(pack(P,gas),grid,solid,gas,badBC),'RocketAeroHeat:NSDraft:OutflowMach');
% 独立解析强膨胀：rho=1,T=300,uL=-a,uR=a，Davis波速为±2a。
% 此时SM=0，p*=p-rho*a^2=(1-gamma)*p=-.4p；所选rho*=rho/2、
% T*=264K仍为正。只查Q*内能会漏过负法向动量通量，必须拒绝。
expansion=zeros(2,1,1,5);expansion(:,:,:,1)=1;expansion(:,:,:,5)=300;
sound=sqrt(gas.gamma*gas.R*300);expansion(1,1,1,2)=-sound;expansion(2,1,1,2)=sound;
pressure=gas.R*300;counterPressure=pressure-sound^2;
counterTemperature=300*(1+(.5*gas.gamma-1)*(gas.gamma-1));
assert(abs(counterPressure/pressure+.4)<1e-12 && abs(counterTemperature-264)<1e-12);
reject(@()compressibleNavierStokesResidual3D(pack(expansion,gas),grid,false(2,1,1),gas,bc), ...
    'RocketAeroHeat:NSDraft:HLLCStarPressure');
assert(~c.physicalSolution && ~c.eligibleForSurfaceCoverage && c.positivityClippingCount==0);
report=struct('status','tested experimental spatial library; not a solved regional flow field', ...
    'eligibleForSurfaceCoverage',false,'physicalSolution',false,'experimentalValidation',false, ...
    'uniformSupersonicResidualMaximum',freeError,'couetteEnergyRelativeError',couetteError, ...
    'linearFourierWallRelativeError',thermalError,'shockFarFluxRelativeError',shockConservationError, ...
    'closedFiveQuantityRelativeError',closedConservationError,'couetteWallPowerW',c.volumeIntegratedResidual(5), ...
    'expectedCouetteWallPowerW',traction*Utop*wallArea,'linearWallHeatFluxWm2',expectedFlux, ...
    'strongExpansionStarPressureRejected',true,'counterexampleStarPressureRatio',counterPressure/pressure, ...
    'counterexampleSelectedStarTemperatureK',counterTemperature, ...
    'testsPassed',9,'scope','仅空间算子独立校核；推进由独立测试验证；不验证实际区域或全面积', ...
    'timeAdvancementTestedHere',false,'meshConvergenceStudy',false);
save(fullfile(outputFolder,'kernel_test_evidence.mat'),'report','c','t','b');
fid=fopen(fullfile(outputFolder,'kernel_test_summary.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
fprintf('3D viscous NS spatial library passed 9 independent groups: free flow, solid NaN/isothermal wall, Couette wall power, linear Fourier heat, contact, shock conservation, closed five-field balance, SI scaling, rejection including negative HLLC star pressure. No solved regional flow or coverage.\n');
disp(report);
end
function Q=pack(P,gas)
rho=P(:,:,:,1);u=P(:,:,:,2);v=P(:,:,:,3);w=P(:,:,:,4);T=P(:,:,:,5);
Q=P;Q(:,:,:,2)=rho.*u;Q(:,:,:,3)=rho.*v;Q(:,:,:,4)=rho.*w;
Q(:,:,:,5)=rho.*(gas.R/(gas.gamma-1)*T+.5*(u.^2+v.^2+w.^2));
end
function r=relative(a,b)
r=max(abs(a-b)./max(1,abs(b)));
end
function expected=independentWallIntegral(P,solid,h,gas,Tw)
% 模块作用：不用被测faceFlux/wallFaces，按流体单元六方向直接积分解析墙通量。
shape=[size(solid,1),size(solid,2),size(solid,3)];expected=zeros(1,5);
for id=find(~solid(:))'
    [i,j,k]=ind2sub(shape,id);sub=[i,j,k];state=reshape(P(i,j,k,:),1,5);pressure=state(1)*gas.R*state(5);
    for axis=1:3
        for sign=[-1 1]
            next=sub;next(axis)=next(axis)+sign;outside=next(axis)<1 || next(axis)>shape(axis);
            if ~outside && ~solid(next(1),next(2),next(3)),continue;end
            area=prod(h)/h(axis);distance=h(axis)/2;n=zeros(1,3);n(axis)=sign;
            drag=gas.mu*state(2:4)/distance;drag(axis)=(4/3)*drag(axis);
            heat=gas.k*(state(5)-Tw)/distance;
            expected=expected+area*[0,pressure*n+drag,heat];
        end
    end
end
end
function reject(action,id)
caught=false;try,action();catch exception,assert(strcmp(exception.identifier,id),exception.message);caught=true;end
assert(caught,['Expected rejection: ' id]);
end
