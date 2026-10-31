function result = solveAxisymmetricEuler(grid,freestream,options)
%SOLVEAXISYMMETRICEULER 内部求解零攻角轴对称、无旋流 Euler 守恒方程。
% 默认一阶有限体积，可选受限线性重构；HLL 或上风总焓能量试验通量。
% 反射滑移壁面 + 显式局部伪时间推进，支持同工况同网格重启。
% 定 gamma 完全气体；局部时间只用于求稳态，不能解释为物理瞬态。
% 输出始终带收敛/边界诊断；未通过者不能当作边界层入口。
if nargin<3,options=struct;end
for name={'gamma','R','pressure','temperature','density','speed','Mach'}
    validateattributes(freestream.(name{1}),{'numeric'},{'real','finite','scalar','positive'});
end
if freestream.gamma<=1 || freestream.Mach<=1
    error('RocketAeroHeat:Euler:Domain','当前入口要求超声速来流和 gamma>1。');
end
for name={'alpha','beta'}
    if isfield(freestream,name{1})
        validateattributes(freestream.(name{1}),{'numeric'},{'real','finite','scalar'});
        if abs(freestream.(name{1}))>1e-10
            error('RocketAeroHeat:Euler:Incidence','轴对称 Euler 仅适用于零攻角、零侧滑角。');
        end
    end
end
g=freestream.gamma; a0=sqrt(g*freestream.R*freestream.temperature);
rho0=freestream.density; L=grid.length;
if abs(freestream.pressure/(rho0*freestream.R*freestream.temperature)-1)>1e-8 || ...
        abs(freestream.speed/a0/freestream.Mach-1)>1e-8
    error('RocketAeroHeat:Euler:State','来流状态不满足定比热理想气体关系。');
end
order=option(options,'spatialOrder',1);
variables=validatestring(option(options,'reconstructionVariables','primitive'),{'primitive','enthalpy'});
if order==1,variables='primitive';end % 常量重构不使用额外的高阶热力变量选项。
validateattributes(order,{'numeric'},{'scalar','integer','>=',1,'<=',2});
maxIter=option(options,'maxIterations',50000); cfl=option(options,'CFL',.5/order);
tol=option(options,'residualTolerance',2e-5);
scheme=validatestring(option(options,'fluxScheme','HLLenthalpy'),{'HLL','HLLenthalpy'});
validateattributes(maxIter,{'numeric'},{'scalar','integer','>=',1});
validateattributes(cfl,{'numeric'},{'real','finite','scalar','>',0,'<=',0.9});
validateattributes(tol,{'numeric'},{'real','finite','scalar','>',0,'<=',1e-2});

%% 1. 无量纲状态与几何；使用自由来流初始化，不引入工程压力作答案
nx=grid.nx; nr=grid.nr;
reconstruction=[];limiterDiagnostics=struct('minimumLimiter',1,'limitedCellCount',0);
if order==2
    reconstruction=prepareEulerReconstruction(grid,options);geo=reconstruction.geometry;
    % 守恒量代表旋转体积平均值，线性重构必须使用相应积分中心和源项矩。
    grid.volume=geo.volume;grid.meridionalArea=geo.area;
    grid.centerX=geo.centerX;grid.centerR=geo.centerR;
end
Uinf=reshape([1,-freestream.Mach,0,1/(g*(g-1))+freestream.Mach^2/2],1,1,4);
U=repmat(Uinf,nx,nr,1); V=grid.volume/L^3;
% 可从同工况、同网格的既有守恒场继续迭代，用于稳健初始化及加长迭代。
% 这里只恢复状态，不强制修正总焓，也不把旧解的收敛标志沿用到新解。
initialization='均匀自由来流';
restart=struct('used',false,'sourceFluxScheme','','sourceConverged',false,'sourceIterations',0);
restart.sourceSpatialOrder=NaN;
if isfield(options,'initialState')
    previous=options.initialState;
    if ~isequal(previous.grid.X,grid.X) || ~isequal(previous.grid.R,grid.R) || ...
            ~strcmp(option(previous.grid,'boundaryTopology','upstreamInlet'), ...
            option(grid,'boundaryTopology','upstreamInlet'))
        error('RocketAeroHeat:Euler:RestartGrid','重启场与当前网格不一致。');
    end
    for name={'gamma','R','pressure','temperature','density','speed','Mach'}
        if ~isequal(previous.freestream.(name{1}),freestream.(name{1}))
            error('RocketAeroHeat:Euler:RestartFlight','重启场与当前来流不一致。');
        end
    end
    for name={'density','pressure','velocityX','velocityR'}
        validateattributes(previous.(name{1}),{'numeric'},{'real','finite','size',[nx,nr]});
    end
    for name={'alpha','beta'}
        if option(previous.freestream,name{1},0)~=option(freestream,name{1},0)
            error('RocketAeroHeat:Euler:RestartFlight','重启场的攻角/侧滑角不一致。');
        end
    end
    if any(previous.density<=0,'all') || any(previous.pressure<=0,'all')
        error('RocketAeroHeat:Euler:RestartState','重启状态的密度/压力必须为正。');
    end
    rho=previous.density/rho0;u=previous.velocityX/a0;v=previous.velocityR/a0;
    pressure=previous.pressure/(rho0*a0^2);
    U=cat(3,rho,rho.*u,rho.*v,pressure/(g-1)+.5*rho.*(u.^2+v.^2));
    initialization='同网格同工况的既有 SI 流场';
    restart.used=true;restart.sourceFluxScheme=option(previous,'fluxScheme','未记录');
    restart.sourceConverged=option(previous,'converged',false);
    restart.sourceIterations=option(previous,'iterations',NaN);
    restart.sourceSpatialOrder=option(previous,'spatialOrder',1);
    oldBoundary='oneSided';if restart.sourceSpatialOrder==1,oldBoundary='firstOrderReflection';end
    restart.sourceBoundaryReconstruction=option(previous,'boundaryReconstruction',oldBoundary);
    restart.sourceReconstructionVariables=option(previous,'reconstructionVariables','primitive');
end
Ax=grid.axialArea/L^2; Bx=grid.radialAreaX/L^2; Br=grid.radialAreaR/L^2;
Ar=option(grid,'axialAreaR',zeros(size(Ax)))/L^2;
boundaryTopology=validatestring(option(grid,'boundaryTopology','upstreamInlet'),{'upstreamInlet','axisAtNose'});
axisAtNose=strcmp(boundaryTopology,'axisAtNose');
if axisAtNose && (any(grid.R(end,:)~=0) || any(Ax(end,:)~=0) || any(Ar(end,:)~=0))
    error('RocketAeroHeat:Euler:AxisBoundary','鼻前轴线必须满足 R=0 且旋转面积为零。');
end
S=grid.meridionalArea/L^2;
scale=reshape([1,max(freestream.Mach,1),max(freestream.Mach,1),Uinf(4)],1,1,4);
history=nan(maxIter,2); timer=tic; rejected=0; converged=false;
fallbackFaces=0;

%% 2. 稳态迭代：密度/压力不为正时缩短步长重试，绝不裁剪守恒量
for iter=1:maxIter
    [res,rate]=residual(U);
    normR=max(abs(res./V./scale),[],'all');
    history(iter,:)=[iter,normR];
    if isfield(options,'progressCallback') && (iter==1 || mod(iter,100)==0)
        options.progressCallback(iter,normR);
    end
    if normR<tol,converged=true;break;end
    dt=cfl*V./rate;
    accepted=false;
    for trial=1:20
        candidate=U-dt./V.*res;
        [rho,~,~,p]=primitive(candidate,g);
        if order==2 && all(isfinite(candidate),'all') && all(rho>1e-10,'all') && all(p>1e-12,'all')
            % 受限线性空间离散采用 SSP-RK2 局部伪时间步；同一 dt 用于
            % 两阶段。单用前向 Euler 可能放大平滑波，不能把限幅器当时间稳定器。
            [secondResidual,~]=residual(candidate);
            candidate=.5*U+.5*(candidate-dt./V.*secondResidual);
            [rho,~,~,p]=primitive(candidate,g);
        end
        if all(isfinite(candidate),'all') && all(rho>1e-10,'all') && all(p>1e-12,'all')
            accepted=true;break;
        end
        dt=0.5*dt; rejected=rejected+1;
    end
    if ~accepted
        error('RocketAeroHeat:Euler:Positivity', ...
            '第 %d 步减小步长后仍出现非正状态（rho_min=%.3g，p_min=%.3g），停止计算。', ...
            iter,min(rho,[],'all'),min(p,[],'all'));
    end
    U=candidate;
end
% 最后一次更新后的残差也必须计算，历史尾行对应实际返回场。
[res,~,Fx,Fb,faceStates,pressureSource]=residual(U);
normR=max(abs(res./V./scale),[],'all'); converged=normR<tol;
history=history(1:iter,:); history(end,2)=normR;

%% 3. 恢复 SI 场及壁面数值通量；实际捕捉的激波由压力梯度显示
[rho,u,v,p,a]=primitive(U,g);
result=struct('grid',grid,'density',rho*rho0,'velocityX',u*a0, ...
    'velocityR',v*a0,'pressure',p*rho0*a0^2, ...
    'temperature',p./rho*a0^2/freestream.R,'Mach',hypot(u,v)./a);
wallArea=hypot(Bx(:,1),Br(:,1));
wallPressure=nan(nx,1); active=grid.bodyWall & wallArea>0;
wallPressure(active)=(Fb(active,1,2).*Bx(active,1)+Fb(active,1,3).*Br(active,1)) ...
    ./wallArea(active).^2*rho0*a0^2;
result.wall=struct('x',grid.wallX,'radius',grid.wallR,'body',active, ...
    'pressure',wallPressure,'cellPressure',result.pressure(:,1));
if order==2
    [rw,uw,vw,pw]=primitive(faceStates.wall,g);
    result.wall.x=geo.faceBX(:,1);result.wall.radius=geo.faceBR(:,1);
    result.wall.trace=struct('density',rw*rho0,'velocityX',uw*a0,'velocityR',vw*a0, ...
        'pressure',pw*rho0*a0^2,'temperature',pw./rw*a0^2/freestream.R, ...
        'x',geo.faceBX(:,1),'radius',geo.faceBR(:,1),'fractionFromTail',geo.faceBFraction(:,1), ...
        'method','受限线性单边壁面迹；r 加权面位置；滑移反射前的内侧状态');
end
result.freestream=freestream; result.converged=converged;
result.iterations=iter; result.residual=normR; result.residualHistory=history;
result.elapsedSeconds=toc(timer); result.positivityRetries=rejected;
result.method=sprintf('轴对称 Euler PDE；一阶有限体积 %s 通量；滑移壁面；局部伪时间',scheme);
if strcmp(scheme,'HLLenthalpy')
    result.method='轴对称 Euler PDE；HLL 质量/动量 + 上风总焓能量通量（试验）；一阶有限体积、滑移壁、局部伪时间';
end
result.fluxScheme=scheme;
result.spatialOrder=order;
result.reconstructionVariables=variables;
result.boundaryReconstruction='firstOrderReflection';
if order==2
    result.boundaryReconstruction=reconstruction.boundaryMethod;
    result.method=[strrep(result.method,'一阶有限体积','有限体积') ...
        '；物理坐标受限线性重构（试验，边界/极值可退化一阶）；SSP-RK2 局部伪时间；Venkatakrishnan 型平滑限幅；压力源项积分矩修正'];
    if strcmp(reconstruction.boundaryMethod,'slipGhost')
        result.method=[result.method '；壁面/轴线镜像状态参与梯度与限幅（平面滑移近似）'];
    else
        result.method=[result.method '；边界梯度仅使用单边域内邻居'];
    end
    if strcmp(variables,'enthalpy'),result.method=[result.method '；rho/u/v/H 重构并恢复 p（试验）'];end
end
result.gridMethod=grid.method;
result.initialization=initialization;
result.restart=restart;
result.scope='零攻角/侧滑、单值回转体、定比热理想气体；不含黏性、底部尾迹、高温化学';
result.thermodynamicModel='caloricallyPerfect';
% 有限体积总残差应等于边界通量减去径向几何源；量纲归一后存档。
sumResidual=reshape(sum(res,[1 2]),1,4);
boundary=reshape(sum(Fx(end,:,:)-Fx(1,:,:),2)+sum(Fb(:,end,:)-Fb(:,1,:),1),1,4);
source=[0,0,sum(pressureSource.*S,'all'),0];
result.conservation.identityError=max(abs(sumResidual-boundary+source));
result.conservation.netBoundaryMinusSource=boundary-source;
result.conservation.relativeImbalance=abs(boundary-source)./reshape(scale,1,4);
result.conservation.wallMassLeak=max(abs(Fb(:,1,1)));
result.conservation.wallEnergyLeak=max(abs(Fb(:,1,4)));
% 左端外推仅适用于全超声速轴向流出；外边界自由来流必须远离扰动。
exitArea=hypot(Ax(1,:),Ar(1,:));
exitNormalVelocity=-(u(1,:).*Ax(1,:)+v(1,:).*Ar(1,:))./exitArea;
result.boundary.supersonicOutflow=all(exitNormalVelocity>a(1,:));
if order==2
    [~,ue,ve,~,ae]=primitive(faceStates.left,g);
    result.boundary.supersonicOutflow=result.boundary.supersonicOutflow && ...
        all(-(ue.*Ax(1,:)+ve.*Ar(1,:))./exitArea>ae);
end
result.boundary.outerPressureDisturbance=max(abs(p(:,end)*g-1));
result.boundary.upstreamPressureDisturbance=max(abs(p(end,:)*g-1));
% 压力平静不能排除熵波/速度扰动，因此同时检查全部原始状态。
disturbance=max(cat(3,abs(rho-1),abs(p*g-1), ...
    abs(u+freestream.Mach)/freestream.Mach,abs(v)/freestream.Mach),[],3);
result.boundary.outerStateDisturbance=max(disturbance(:,end));
result.boundary.upstreamStateDisturbance=max(disturbance(end,:));
if order==2
    result.boundary.outerStateDisturbance=max(result.boundary.outerStateDisturbance,stateDisturbance(faceStates.outer,g,freestream.Mach));
    result.boundary.upstreamStateDisturbance=max(result.boundary.upstreamStateDisturbance,stateDisturbance(faceStates.right,g,freestream.Mach));
end
result.boundary.topology=boundaryTopology;
if axisAtNose
    % 鼻前轴线上的流动本来就含激波和滞止减速，不应要求其等于来流。
    % 鼻前自由来流包含在外层曲边边界检查中；NaN 明确表示没有独立入口边。
    result.boundary.upstreamPressureDisturbance=NaN;
    result.boundary.upstreamStateDisturbance=NaN;
end
result.boundary.clear= result.boundary.outerStateDisturbance<0.01 && ...
    (axisAtNose || result.boundary.upstreamStateDisturbance<0.01) && result.boundary.supersonicOutflow;
result.available=converged && result.boundary.clear;
result.edgeStateAvailable=false; % 尚未完成壁面外缘速度/熵的网格精度验证
% 绝热无外功的稳态流线应保持总焓；一阶数值耗散仍可能引起局部偏差。
% 单独存储此精度诊断，不能把残差变小等同于热力学精度足够。
H=g/(g-1)*p./rho+0.5*(u.^2+v.^2);
Hinf=1/(g-1)+0.5*freestream.Mach^2;
result.diagnostics.totalEnthalpyRelativeError=H/Hinf-1;
result.diagnostics.maximumTotalEnthalpyDeviation=max(abs(H/Hinf-1),[],'all');
result.diagnostics.gridAccuracyVerified=false;
result.diagnostics.fluxFallbackFaceCount=fallbackFaces;
result.diagnostics.reconstruction=limiterDiagnostics;
result.status=sprintf('收敛=%d；边界合格=%d；残差=%.3g；外边界扰动=%.3g', ...
    converged,result.boundary.clear,normR,result.boundary.outerPressureDisturbance);

    function [res,rate,Fx,Fb,trace,sourceP]=residual(Q)
        % 轴向左端流出、右端超声速流入；径向外端固定自由来流。
        [rhoQ,uQ,vQ,pressure]=primitive(Q,g);
        if order==2
            if strcmp(variables,'enthalpy')
                [reconstructed,sourceP,limiterDiagnostics]=reconstructEulerEnthalpy(cat(3,rhoQ,uQ,vQ,pressure),reconstruction,g);
            else
                [reconstructed,sourceP,limiterDiagnostics]=reconstructEulerPrimitive(cat(3,rhoQ,uQ,vQ,pressure),reconstruction);
            end
            qm=conservative(reconstructed{1},g);qp=conservative(reconstructed{2},g);
            wm=conservative(reconstructed{3},g);wp=conservative(reconstructed{4},g);
            right=repmat(Uinf,1,nr,1);
            if axisAtNose,right=qp(end,:,:);right(:,:,3)=-right(:,:,3);end
            [Fx,rx,fxFallback]=axisymmetricEulerFlux(cat(1,qm(1,:,:),qp),cat(1,qm,right),Ax,Ar,g,scheme);
            ghost=reflect(wm(:,1,:),Bx(:,1),Br(:,1));
            [Fb,rb,fbFallback]=axisymmetricEulerFlux(cat(2,ghost,wp),cat(2,wm,repmat(Uinf,nx,1,1)),Bx,Br,g,scheme);
            trace=struct('left',qm(1,:,:),'right',qp(end,:,:),'wall',wm(:,1,:),'outer',wp(:,end,:));
        else
        right=repmat(Uinf,1,nr,1);
        if axisAtNose
            right=Q(end,:,:);right(:,:,3)=-right(:,:,3);
        end
        Qx=cat(1,Q(1,:,:),Q,right);
        [Fx,rx,fxFallback]=axisymmetricEulerFlux(Qx(1:end-1,:,:),Qx(2:end,:,:),Ax,Ar,g,scheme);
        ghost=Q(:,1,:); den=hypot(Bx(:,1),Br(:,1));
        % 轴线面加权面积为零；给任意有限单位法向，通量自动为零。
        nnx=zeros(nx,1); nnr=ones(nx,1); nonzero=den>0;
        nnx(nonzero)=Bx(nonzero,1)./den(nonzero); nnr(nonzero)=Br(nonzero,1)./den(nonzero);
        mn=ghost(:,:,2).*nnx+ghost(:,:,3).*nnr;
        ghost(:,:,2)=ghost(:,:,2)-2*mn.*nnx;
        ghost(:,:,3)=ghost(:,:,3)-2*mn.*nnr;
        Qb=cat(2,ghost,Q,repmat(Uinf,nx,1,1));
        [Fb,rb,fbFallback]=axisymmetricEulerFlux(Qb(:,1:end-1,:),Qb(:,2:end,:),Bx,Br,g,scheme);
        sourceP=pressure;
        trace=struct('left',Q(1,:,:),'right',Q(end,:,:),'wall',Q(:,1,:),'outer',Q(:,end,:));
        end
        fallbackFaces=fxFallback+fbFallback;
        res=diff(Fx,1,1)+diff(Fb,1,2);
        res(:,:,3)=res(:,:,3)-sourceP.*S;
        rate=rx(1:end-1,:)+rx(2:end,:)+rb(:,1:end-1)+rb(:,2:end);
    end
end
function Q=conservative(P,g)
rho=P(:,:,1);u=P(:,:,2);v=P(:,:,3);p=P(:,:,4);
Q=cat(3,rho,rho.*u,rho.*v,p/(g-1)+.5*rho.*(u.^2+v.^2));
end
function ghost=reflect(Q,sx,sr)
den=hypot(sx,sr);nx=zeros(size(sx));nr=ones(size(sr));active=den>0;
nx(active)=sx(active)./den(active);nr(active)=sr(active)./den(active);
ghost=Q;mn=Q(:,:,2).*nx+Q(:,:,3).*nr;
ghost(:,:,2)=Q(:,:,2)-2*mn.*nx;ghost(:,:,3)=Q(:,:,3)-2*mn.*nr;
end
function value=stateDisturbance(Q,g,M)
[rho,u,v,p]=primitive(Q,g);
value=max(cat(3,abs(rho-1),abs(p*g-1),abs(u+M)/M,abs(v)/M),[],'all');
end
function [rho,u,v,p,a]=primitive(U,g)
rho=U(:,:,1); u=U(:,:,2)./rho; v=U(:,:,3)./rho;
p=(g-1)*(U(:,:,4)-0.5*rho.*(u.^2+v.^2));
if nargout>4,a=sqrt(g*p./rho);end
end
function value=option(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
