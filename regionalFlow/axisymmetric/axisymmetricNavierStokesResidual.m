function [rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary,saState)
%AXISYMMETRICNAVIERSTOKESRESIDUAL 贴体无旋轴对称四守恒量NS空间算子。
% 模块一：定gamma/常Pr，显式常数或Sutherland输运；真实四守恒量，不裁状态。
% 模块二：共享面HLLC一阶对流、WLS+位移约束黏性梯度、Fourier导热。
% 模块三：源项[0,0,integral(p-tau_theta_theta)dA,0]使用平面面积；
% 速度散度保留ur/r，轴边界用正则极限，不能用eps代替半径。
% 模块四：Tw及储库始终为输入；输出当前墙q、机械功和完整边界/源项账。
% 边界iMin/iMax/jMin/jMax全部显式：wall(T,velocity=[ux ur])、axis、
% periodic(仅几何显式配对的i)、supersonicInflow(primitive=[rho ux ur T])、
% supersonicOutflow、characteristic(复用三维特征边界，储库为五原始量)。
% riemannReservoir：明确四原始量储库，开放面统一使用HLLE两波数值通量；
% 黏性入流给储库迹，按实际质量通量判定的出流为零法向梯度，零质量按出流。
% prescribedTrace仅供有完整解析边界的数值校核，不能冒充通用亚声速入口。
% 可显式指定wall.thermalGradient='dirichlet-tangent'；函数Tw须给沿t=[-nr,nx]的导数。
% 可显式指定静止wall.velocityGradient='dirichlet-tangent'约束壁速度切向导数0。
% 可显式指定全局boundary.viscousWorkVelocity='arithmetic'，只改内部面黏性功速度。
% 可显式指定全局boundary.pressureSourceQuadrature='cell-constant'，用分片常数p积分几何源。
% 默认仍为原linear-planar；显式选择才附带积分审计，不按压力正负自动切换。
% 可显式指定boundary.interiorConvectiveFlux为hllc/hlle-davis/hlle-einfeldt。
% 仅改内部Euler通量；默认HLLC和开放面HLLE合同保持，SA使用实际新质量通量。
% 可显式指定boundary.interiorFluxArithmetic='two-double-euler'，仅提高内部HLLE形成精度。
% 该选项不提高黏性/边界/几何源的精度，不是新的物理模型；省略保持原完整数值路径。
% axis.symmetryGradient='mirror-parity'显式启用偏斜轴镜像WLS，真实面和虚样本分开；
% 零面积轴面本构状态用owner偶量投影代理，原始切向外推仅诊断而非已解轴迹。
% 当前提供层流及显式SA耦合的数值模块，未通过三维横流/真实箭体物理验收。
% 显式第五参量仅由五量RANS残差使用；省略时保留原四量表达式和完整info。
turbulent=nargin>=5;
if turbulent
    assert(isstruct(saState)&&isscalar(saState)&&isequal(sort(fieldnames(saState)),sort({'nuTilde';'wallDistance'})), ...
        'AxisSA:State','SA内部耦合需要明确N与真实壁距。');
end
gas=axisymmetricFVGas(gas);P=axisymmetricFVState(Q,gas);nc=mesh.cellCount;
variableTransport=isfield(gas,'transport');
assert(~turbulent||variableTransport,'AxisSA:Gas','SA需要显式局部输运合同。');
assert(strcmp(mesh.kind,'axisymmetric-body-fitted-fv-v1')&&size(Q,1)==nc, ...
    'AxisFV:Geometry','状态须匹配本贴体控制体。');
f=mesh.faces;owner=f.owner;other=f.neighbour;nf=numel(owner);inside=other>0;
V=mesh.volumePerRadian;A=mesh.planarArea;C=mesh.incidence;
assert(all(isfinite(V)&V>0)&&all(isfinite(A)&A>0)&&all(isfinite(f.areaPerRadian)&f.areaPerRadian>=0)&& ...
    all(isfinite(f.normal),'all')&&max(abs(sum(f.normal.^2,2)-1))<1e-10&& ...
    max(abs(f.areaVectorPerRadian-f.normal.*f.areaPerRadian),[],'all')<=1e-12*max(1,max(f.areaPerRadian))&& ...
    max(abs(C*f.areaVectorPerRadian-[zeros(nc,1),A])./A,[],'all')<1e-9, ...
    'AxisFV:Geometry','几何度量改变或无效。');
names={'iMin','iMax','jMin','jMax'};
assert(isstruct(boundary)&&all(isfield(boundary,names)),'AxisFV:Boundary','须明确四条边的边界条件。');
if mesh.periodicI
    assert(strcmp(boundary.iMin.type,'periodic')&&strcmp(boundary.iMax.type,'periodic'), ...
        'AxisFV:Boundary','几何配对周期须与边界一致。');
else
    assert(~strcmp(boundary.iMin.type,'periodic')&&~strcmp(boundary.iMax.type,'periodic'), ...
        'AxisFV:Boundary','周期条件需要先配对几何。');
end

axisSymmetry=axisymmetricAxisSymmetry(mesh,boundary);

%% 一、明确边界迹及梯度约束；不可由表面外缘标量虚构整个入口剖面
phi=P(:,[2 3 4 5]);bf=mesh.boundaryFaces;nb=numel(bf);trace=P(owner(bf),:);
kind=zeros(nf,1);wall=false(nf,1);axisFace=false(nf,1);normalZero=false(nf,1);
wallInput=nan(nf,3);incoming=false(nf,1);
riemannFace=false(nf,1);riemannState=nan(nf,6);riemannFlux=zeros(nf,4);
constantReservoir=nan(4,6);reservoirCached=false(4,1);
for b=1:nb
    face=bf(b);id=owner(face);n=f.normal(face,:);point=f.centroid(face,:);bc=boundary.(names{f.boundarySide(face)});
    assert(isstruct(bc)&&isscalar(bc)&&isfield(bc,'type'),'AxisFV:Boundary','边界合同无效。');
    switch bc.type
        case 'wall'
            assert(all(isfield(bc,{'temperature','velocity'})),'AxisFV:Boundary','墙温和墙速度须明确给定。');
            value=bc.temperature;if isa(value,'function_handle'),value=value(point);end
            velocity=bc.velocity;if isa(velocity,'function_handle'),velocity=velocity(point);end
            validateattributes(value,{'numeric'},{'real','finite','scalar','>=',gas.temperatureRange(1),'<=',gas.temperatureRange(2)});
            validateattributes(velocity,{'numeric'},{'real','finite','numel',2});velocity=reshape(velocity,1,2);
            assert(abs(dot(velocity,n))<1e-12*max(1,norm(velocity)),'AxisFV:WallMotion','固定网格墙只允许切向运动。');
            trace(b,2:4)=[velocity,value];trace(b,1)=trace(b,5)/(gas.R*value);
            wallInput(face,:)=[velocity,value];wall(face)=true;kind(face)=1;incoming(face)=true;
        case 'axis'
            assert(f.areaPerRadian(face)==0&&point(2)==0&&abs(n(1))<1e-12, ...
                'AxisFV:Axis','轴边界必须位于真实r=0，不能指定在有限半径。');
            if ~axisSymmetry.active||~axisSymmetry.faceMask(face)
                assert(abs(point(1)-mesh.volumeCentroid(id,1))<=1e-10*abs(point(2)-mesh.volumeCentroid(id,2)), ...
                    'AxisFV:AxisSkew','缺省轴边界要求首层正交；偏斜轴须显式选择mirror-parity并验证。');
            end
            trace(b,3)=0;axisFace(face)=true;kind(face)=2;incoming(face)=true;
        case {'supersonicInflow','prescribedTrace'}
            assert(isfield(bc,'primitive'),'AxisFV:Boundary','入口/解析迹缺少完整原始量。');
            value=bc.primitive;if isa(value,'function_handle'),value=value(point);end
            state=axisymmetricFVState(axisymmetricFVConservative(reshape(value,1,[]),gas),gas);
            if strcmp(bc.type,'supersonicInflow')
                assert(dot(state(2:3),n)<-state(6),'AxisFV:InflowMach','输入不是全特征向内的超声速。');kind(face)=3;
            else,kind(face)=4;end
            trace(b,:)=state;incoming(face)=true;
        case 'supersonicOutflow'
            assert(dot(trace(b,2:3),n)>trace(b,6),'AxisFV:OutflowMach','出口存在入射特征/回流，须使用明确特征边界。');
            normalZero(face)=true;kind(face)=5;
        case 'characteristic'
            assert(exist('resolveNSCharacteristicBoundary3D','file')==2,'AxisFV:Dependency','请将internalFlow加入路径。');
            if isfield(bc,'reservoirPrimitive')&&isa(bc.reservoirPrimitive,'function_handle'),bc.reservoirPrimitive=bc.reservoirPrimitive(point);end
            assert(isfield(bc,'reservoirPrimitive')&&numel(bc.reservoirPrimitive)==5&&bc.reservoirPrimitive(4)==0, ...
                'AxisFV:Swirl','当前无旋轴对称接口不接受非零周向储库速度。');
            [state,characteristic]=resolveNSCharacteristicBoundary3D([trace(b,1:3),0,trace(b,4)],bc,[n,0],gas);
            trace(b,:)=axisymmetricFVState(axisymmetricFVConservative(state([1 2 3 5]),gas),gas);
            incoming(face)=characteristic.isInflow;normalZero(face)=~characteristic.isInflow;kind(face)=6;
        case 'riemannReservoir'
            assert(isfield(bc,'primitive'),'AxisFV:RiemannReservoir','Riemann储库须明确四原始量[rho ux ur T]。');
            sideID=f.boundarySide(face);constant=~isa(bc.primitive,'function_handle');
            if constant&&reservoirCached(sideID)
                exterior=constantReservoir(sideID,:);
            else
                value=bc.primitive;if ~constant,value=value(point);end
                assert(isnumeric(value)&&isreal(value)&&isvector(value)&&numel(value)==4&&all(isfinite(value(:))), ...
                    'AxisFV:RiemannReservoir','无旋储库必须是四个有限实数，不能省略或加入周向速度。');
                exteriorQ=axisymmetricFVConservative(reshape(value,1,4),gas);
                exterior=axisymmetricFVState(exteriorQ,gas);
                if constant,constantReservoir(sideID,:)=exterior;reservoirCached(sideID)=true;end
            end
            riemannState(face,:)=exterior;riemannFace(face)=true;kind(face)=7;
            trace(b,:)=exterior;
        otherwise
            error('AxisFV:Boundary','未知边界类型%s；不自动补入流或出口数据。',bc.type);
    end
end
% 批量开放面Riemann通量；先确定真实质量方向，再建立黏性边界迹。
reservoirFaces=find(riemannFace);
if ~isempty(reservoirFaces)
    exteriorQ=axisymmetricFVConservative(riemannState(reservoirFaces,1:4),gas);
    riemannFlux(reservoirFaces,:)=hlleBoundary(Q(owner(reservoirFaces),:),exteriorQ, ...
        P(owner(reservoirFaces),:),riemannState(reservoirFaces,:),f.normal(reservoirFaces,:));
    incoming(reservoirFaces)=riemannFlux(reservoirFaces,1)<0;
    normalZero(reservoirFaces)=~incoming(reservoirFaces);
    outFaces=reservoirFaces(~incoming(reservoirFaces));
    boundaryIndex=zeros(nf,1);boundaryIndex(bf)=1:nb;
    trace(boundaryIndex(outFaces),:)=P(owner(outFaces),:);
end
trace=axisymmetricFVState(axisymmetricFVConservative(trace(:,1:4),gas),gas);
if axisSymmetry.active
    % 虚样本是完整镜像状态，偶rho/ux/T/p相同、奇ur反号；EOS在镜像点成立。
    ai=axisSymmetry.boundaryIndex;ao=axisSymmetry.ownerCell;dx=axisSymmetry.physicalDisplacementXR(:,1);
    ghostPhi=[P(ao,2),-P(ao,3),P(ao,4:5)];
    gradientBoundary=trace(:,[2 3 4 5]);gradientBoundary(ai,:)=ghostPhi;
    sample=[phi;gradientBoundary];gx=axisSymmetry.gradientX*sample;gr=axisSymmetry.gradientR*sample;
    densityBoundary=trace(:,1);densityBoundary(ai)=P(ao,1);
    densitySamples=[P(:,1);densityBoundary];densityGX=axisSymmetry.gradientX*densitySamples;
    densityGR=axisSymmetry.gradientR*densitySamples;
    % 原始切向外推仅用于审计：鼻/尾端偏斜薄单元可能给出无界迹，不能当物理解。
    rawAxisPrimitive=[P(ao,1)+densityGX(ao).*dx,P(ao,2)+gx(ao,1).*dx, ...
        zeros(numel(ao),1),P(ao,4)+gx(ao,3).*dx];
    % 实际零面积轴面采用一阶owner投影代理：偶rho/ux/T同owner，奇ur=0。
    % 这是显式离散合同，不是raw越界后的限幅/回退；任何raw都不进入通量或局部温界。
    axisPrimitive=[P(ao,1:2),zeros(numel(ao),1),P(ao,4)];
    assert(all(isfinite(axisPrimitive),'all')&&all(axisPrimitive(:,1)>0)&& ...
        all(axisPrimitive(:,4)>=gas.temperatureRange(1)&axisPrimitive(:,4)<=gas.temperatureRange(2)), ...
        'AxisFV:AxisTraceState','零面积轴面owner投影状态非有限、rho非正或温度越域，未裁值。');
    trace(ai,:)=[axisPrimitive,axisPrimitive(:,1).*gas.R.*axisPrimitive(:,4), ...
        sqrt(gas.gamma*gas.R*axisPrimitive(:,4))];
else
    sample=[phi;trace(:,[2 3 4 5])];gx=mesh.gradientX*sample;gr=mesh.gradientR*sample;
end
assert(all(isfinite(gx),'all')&&all(isfinite(gr),'all'),'AxisFV:Gradient','梯度非有限。');

%% 二、每个共享面仅一次Riemann/黏性通量，使用真实体积与面权重
left=P(owner,:);right=left;right(inside,:)=P(other(inside),:);right(bf,:)=trace;
QL=Q(owner,:);QR=QL;QR(inside,:)=Q(other(inside),:);QR(bf,:)=axisymmetricFVConservative(trace(:,1:4),gas);
interiorFluxRequested=isfield(boundary,'interiorConvectiveFlux');interiorPolicy='hllc';
if interiorFluxRequested,interiorPolicy=boundary.interiorConvectiveFlux;end
assert(((ischar(interiorPolicy)&&isrow(interiorPolicy))||(isstring(interiorPolicy)&&isscalar(interiorPolicy)&&~ismissing(interiorPolicy)))&& ...
    ismember(char(interiorPolicy),{'hllc','hlle-davis','hlle-einfeldt'}), ...
    'AxisFV:InteriorFluxPolicy','内部对流通量只支持hllc/hlle-davis/hlle-einfeldt。');interiorPolicy=char(interiorPolicy);
for side=1:4
    assert(~isfield(boundary.(names{side}),'interiorConvectiveFlux'), ...
        'AxisFV:InteriorFluxPolicy','内部通量是全域离散合同，不可放在单条边界内而被忽略。');
    assert(~isfield(boundary.(names{side}),'interiorFluxArithmetic'), ...
        'AxisFV:InteriorArithmetic','内部通量算术必须是全域选项，不能放在单条边界。');
end
arithmeticRequested=isfield(boundary,'interiorFluxArithmetic');arithmetic='double';
if arithmeticRequested
    arithmetic=boundary.interiorFluxArithmetic;
    assert(((ischar(arithmetic)&&isrow(arithmetic))||(isstring(arithmetic)&&isscalar(arithmetic)&&~ismissing(arithmetic)))&& ...
        ismember(char(arithmetic),{'double','two-double-euler'}),'AxisFV:InteriorArithmetic','未知内部通量算术。');arithmetic=char(arithmetic);
    assert(~strcmp(arithmetic,'two-double-euler')||ismember(interiorPolicy,{'hlle-davis','hlle-einfeldt'}), ...
        'AxisFV:InteriorArithmetic','双double形成仅支持显式HLLE，不自动替换HLLC。');
end
flux=physical(QL,left,f.normal);
if strcmp(interiorPolicy,'hllc')
    flux(inside,:)=hllc(QL(inside,:),QR(inside,:),left(inside,:),right(inside,:),f.normal(inside,:));
    if interiorFluxRequested
        interiorUnL=sum(left(inside,2:3).*f.normal(inside,:),2);interiorUnR=sum(right(inside,2:3).*f.normal(inside,:),2);
        interiorConvectiveFlux=struct('kind','explicit-original-interior-HLLC-flux-v1', ...
            'waveSpeedsMPerS',[min(interiorUnL-left(inside,6),interiorUnR-right(inside,6)),max(interiorUnL+left(inside,6),interiorUnR+right(inside,6))], ...
            'fluxPerArea',flux(inside,:),'stateClipped',false,'adaptiveFluxFallback',false);
    end
else
    wavePolicy=extractAfter(interiorPolicy,'hlle-');
    if strcmp(arithmetic,'two-double-euler')
        [flux(inside,:),interiorConvectiveFlux]=axisymmetricHLLEFluxAccurate(QL(inside,:),QR(inside,:),f.normal(inside,:),gas,wavePolicy);
    else
        [flux(inside,:),interiorConvectiveFlux]=axisymmetricHLLEFlux(QL(inside,:),QR(inside,:),f.normal(inside,:),gas,wavePolicy);
    end
end
if interiorFluxRequested
    interiorConvectiveFlux.policy=interiorPolicy;interiorConvectiveFlux.faceID=find(inside);
    interiorConvectiveFlux.scope='interior Euler flux only; original open-face solver, viscous flux and geometric sources unchanged; SA uses actual mass flux';
end
open=bf(~wall(bf)&~axisFace(bf));flux(open,:)=physical(QR(open,:),right(open,:),f.normal(open,:));
flux(riemannFace,:)=riemannFlux(riemannFace,:);
flux(wall|axisFace,:)=0;flux(wall,2:3)=left(wall,5).*f.normal(wall,:);
if any(riemannFace),convectiveFlux=flux;end
Gx=gx(owner,1:3);Gr=gr(owner,1:3);facePhi=left(:,2:4);d=f.delta;
Gx(inside,:)=.5*(gx(owner(inside),1:3)+gx(other(inside),1:3));
Gr(inside,:)=.5*(gr(owner(inside),1:3)+gr(other(inside),1:3));
% 用两侧线性迹平均求黏性功速度，线性场在偏斜和r加权质心上仍精确。
dl=f.centroid(inside,:)-mesh.volumeCentroid(owner(inside),:);
dr=f.centroid(inside,:)-f.neighbourShift(inside,:)-mesh.volumeCentroid(other(inside),:);
facePhi(inside,:)=.5*(phi(owner(inside),1:3)+gx(owner(inside),1:3).*dl(:,1)+gr(owner(inside),1:3).*dl(:,2)+ ...
    phi(other(inside),1:3)+gx(other(inside),1:3).*dr(:,1)+gr(other(inside),1:3).*dr(:,2));
facePhi(bf,:)=trace(:,2:4);
diff=right(:,2:4)-left(:,2:4);correction=(diff-Gx.*d(:,1)-Gr.*d(:,2))./sum(d.^2,2);
Gx=Gx+correction.*d(:,1);Gr=Gr+correction.*d(:,2);
z=find(normalZero);normalGradient=Gx(z,:).*f.normal(z,1)+Gr(z,:).*f.normal(z,2);
Gx(z,:)=Gx(z,:)-normalGradient.*f.normal(z,1);Gr(z,:)=Gr(z,:)-normalGradient.*f.normal(z,2);
if axisSymmetry.active
    af=axisSymmetry.faceID;ao=axisSymmetry.ownerCell;
    % 偶ux/T径向导数0；奇ur轴值0且轴向导数0，gn=urOwner/rOwner为一阶正则代理。
    Gx(af,:)=[gx(ao,1),zeros(numel(ao),1),gx(ao,3)];
    Gr(af,:)=[zeros(numel(ao),1),P(ao,3)./mesh.volumeCentroid(ao,2),zeros(numel(ao),1)];
end
wallVelocityRequested=false;
for side=1:4
    bc=boundary.(names{side});wallVelocityRequested=wallVelocityRequested|| ...
        isfield(bc,'velocityGradient')||isfield(bc,'velocityTangentialGradient');
end
if wallVelocityRequested
    % 仅显式静止墙的Newton应力梯度受约束；cell WLS、SA源和Fourier梯度保持。
    wids=find(wall);
    [Gx(:,1:2),Gr(:,1:2),wallVelocityGradient]=axisymmetricWallVelocityGradient(mesh,P(:,2:3),boundary, ...
        wids,wallInput(wids,1:2),Gx(:,1:2),Gr(:,1:2));
end
wallThermalRequested=false;
for side=1:4
    bc=boundary.(names{side});wallThermalRequested=wallThermalRequested|| ...
        isfield(bc,'thermalGradient')||isfield(bc,'temperatureTangentialGradient');
end
if wallThermalRequested
    % 改动进入同一次能量残差/BE/Newton，而非对保存的q后处理换色或取绝对值。
    wids=find(wall);
    [Gx(:,3),Gr(:,3),wallThermalGradient]=axisymmetricWallThermalGradient(mesh,P(:,4),boundary, ...
        wids,wallInput(wids,3),Gx(:,3),Gr(:,3));
end
hoop=zeros(nf,1);positive=f.centroid(:,2)>0;
hoop(positive)=facePhi(positive,2)./f.centroid(positive,2);hoop(~positive)=Gr(~positive,2);
div=Gx(:,1)+Gr(:,2)+hoop;
if variableTransport
    % 局部面迹与Newton/Fourier梯度同点：墙面严格使用输入Tw。
    % 边界原始量经守恒量往返可能差一个ulp；输入Tw不能改成该舍入副本。
    facePhi(wall,3)=wallInput(wall,3);
    % 拉伸真实网格的无限制迹可为负；只限制本构采样重构，不裁单元Q或原热梯度。
    transportTrace=axisymmetricFVTemperatureTrace(mesh,P(:,4),gx(:,3),gr(:,3),facePhi(bf,3));
    facePhi(:,3)=transportTrace.faceTemperature;
    faceTransport=axisymmetricFVTransport(facePhi(:,3),gas);
    faceMu=faceTransport.mu;faceK=faceTransport.k;
    if turbulent
        if axisSymmetry.active
            sa=axisymmetricSADiscretization(mesh,P,boundary,right,flux(:,1),gx,gr,gas, ...
                saState.nuTilde,saState.wallDistance,struct('face',facePhi(:,3),'source',transportTrace.sourceTemperature),axisSymmetry);
        else
            sa=axisymmetricSADiscretization(mesh,P,boundary,right,flux(:,1),gx,gr,gas, ...
                saState.nuTilde,saState.wallDistance,struct('face',facePhi(:,3),'source',transportTrace.sourceTemperature));
        end
        faceMu=sa.effectiveFaceMu;faceK=sa.effectiveFaceK;
    end
    txx=2*faceMu.*(Gx(:,1)-div/3);trr=2*faceMu.*(Gr(:,2)-div/3);
    txr=faceMu.*(Gr(:,1)+Gx(:,2));
else
    txx=2*gas.mu*(Gx(:,1)-div/3);trr=2*gas.mu*(Gr(:,2)-div/3);txr=gas.mu*(Gr(:,1)+Gx(:,2));
end
traction=[txx.*f.normal(:,1)+txr.*f.normal(:,2),txr.*f.normal(:,1)+trr.*f.normal(:,2)];
if variableTransport
    heat=-faceK.*(Gx(:,3).*f.normal(:,1)+Gr(:,3).*f.normal(:,2));
else
    heat=-gas.k*(Gx(:,3).*f.normal(:,1)+Gr(:,3).*f.normal(:,2));
end
viscous=[zeros(nf,1),traction,sum(facePhi(:,1:2).*traction,2)-heat];
viscousWorkVelocityRequested=isfield(boundary,'viscousWorkVelocity');
for side=1:4
    viscousWorkVelocityRequested=viscousWorkVelocityRequested||isfield(boundary.(names{side}),'viscousWorkVelocity');
end
if viscousWorkVelocityRequested
    % 原facePhi仍用于应力/hoop/温度，原WLS功保留诊断；只替换内部能量通量的功速度。
    [workVelocity,viscousWorkVelocity]=axisymmetricViscousWorkVelocity(mesh,P(:,2:3),boundary,facePhi(:,1:2),traction);
    if strcmp(viscousWorkVelocity.policy,'arithmetic')
        vids=viscousWorkVelocity.faceID;
        viscous(vids,4)=viscousWorkVelocity.workFluxPerAreaWm2-heat(vids);
    end
end
flux=flux-viscous;flux(axisFace,:)=0;integrated=flux.*f.areaPerRadian;

%% 三、源项在平面质心求积，区别于r加权体质心；轴单元质心仍严格r>0
delta=mesh.planarCentroid-mesh.volumeCentroid;
pA=P(:,5)+gx(:,4).*delta(:,1)+gr(:,4).*delta(:,2);
pressureQuadratureRequested=isfield(boundary,'pressureSourceQuadrature');
for side=1:4
    assert(~isfield(boundary.(names{side}),'pressureSourceQuadrature'), ...
        'AxisFV:PressureQuadraturePolicy','压力源积分是全域离散合同，不能写入单条边界后被忽略。');
end
if pressureQuadratureRequested
    [pA,pressureSourceQuadrature]=axisymmetricPressureSourceQuadrature(mesh,P(:,5),gx(:,4),gr(:,4), ...
        boundary.pressureSourceQuadrature);
end
vA=P(:,3)+gx(:,2).*delta(:,1)+gr(:,2).*delta(:,2);
hoopA=vA./mesh.planarCentroid(:,2);divA=gx(:,1)+gr(:,2)+hoopA;
if variableTransport
    sourceTemperature=transportTrace.sourceTemperature;
    sourceTransport=axisymmetricFVTransport(sourceTemperature,gas);
    sourceMu=sourceTransport.mu;if turbulent,sourceMu=sa.effectiveSourceMu;end
    ttt=2*sourceMu.*(hoopA-divA/3);
else
    ttt=2*gas.mu*(hoopA-divA/3);
end
source=zeros(nc,4);source(:,3)=(pA-ttt).*A;
assert(all(isfinite(pA)&pA>0)&&all(isfinite(source),'all')&&all(isfinite(integrated),'all'), ...
    'AxisFV:PressureSource','平面求积压力/几何源/通量无效，未裁状态。');
rate=(-C*integrated+source)./V;
net=sum(integrated(bf,:),1);src=sum(source,1);defect=sum(rate.*V,1)+net-src;

%% 四、提供完整功率、稳定频率和未获物理验收的状态
waveRight=right;waveRight(riemannFace,:)=riemannState(riemannFace,:);
speed=max(abs(sum(left(:,2:3).*f.normal,2))+left(:,6),abs(sum(waveRight(:,2:3).*f.normal,2))+waveRight(:,6));
if interiorFluxRequested
    % 保留原两侧声学频率下界；Roe界更大时还须包围实际新通量波速，不能低估CFL。
    speed(inside)=max(speed(inside),max(abs(interiorConvectiveFlux.waveSpeedsMPerS),[],2));
    if strcmp(arithmetic,'two-double-euler'),speed(inside)=speed(inside)+eps(speed(inside));end
end
if variableTransport
    % 频率取面及两侧温度的最大系数，不把参考mu/k误用为当前扩散系数。
    cellTransport=axisymmetricFVTransport(P(:,4),gas);
    rightTransport=axisymmetricFVTransport(right(:,4),gas);
    frequencyMu=max([faceTransport.mu,cellTransport.mu(owner),rightTransport.mu],[],2);
    frequencyK=max([faceTransport.k,cellTransport.k(owner),rightTransport.k],[],2);
    if turbulent
        rightMuT=sa.cellMuT(owner);rightMuT(inside)=sa.cellMuT(other(inside));rightMuT(bf)=sa.faceMuT(bf);
        frequencyMu=max([sa.effectiveFaceMu,sa.effectiveCellMu(owner),rightTransport.mu+rightMuT],[],2);
        effectiveCellK=cellTransport.k+gas.cp*sa.cellMuT/sa.cellModel.constants.turbulentPrandtl;
        frequencyK=max([sa.effectiveFaceK,effectiveCellK(owner), ...
            rightTransport.k+gas.cp*rightMuT/sa.cellModel.constants.turbulentPrandtl],[],2);
    end
    nu=max(4*frequencyMu./(3*min(left(:,1),right(:,1))),frequencyK./(gas.cv*min(left(:,1),right(:,1))));
else
    nu=max(4*gas.mu./(3*min(left(:,1),right(:,1))),gas.k./(gas.cv*min(left(:,1),right(:,1))));
end
convective=(mesh.cellFaceMultiplicity*(f.areaPerRadian.*speed))./V;
if variableTransport
    diffusive=2*(mesh.cellFaceMultiplicity*(nu.*f.areaPerRadian./f.normalDistance))./V+ ...
        4*cellTransport.mu./(P(:,1).*mesh.volumeCentroid(:,2).^2);
else
    diffusive=2*(mesh.cellFaceMultiplicity*(nu.*f.areaPerRadian./f.normalDistance))./V+ ...
        4*gas.mu./(P(:,1).*mesh.volumeCentroid(:,2).^2);
end
if turbulent
    % 仅是当前伪步估计，不宣称任意非线性状态的无条件稳定性。
    diffusive=diffusive+4*sa.cellMuT./(P(:,1).*mesh.volumeCentroid(:,2).^2);
    scalarNu=sa.diffusionCoefficient./min(left(:,1),right(:,1));
    sa.diffusiveFrequency=2*(mesh.cellFaceMultiplicity*(scalarNu.*f.areaPerRadian./f.normalDistance))./V;
    sa.sourceFrequency=sa.cellModel.constants.cb1*abs(sa.cellModel.sTilde)+ ...
        2*sa.cellModel.constants.cw1*sa.cellModel.fw.*sa.nuTilde./sa.wallDistance.^2;
    diffusive=max(diffusive,sa.diffusiveFrequency+sa.sourceFrequency);
end
w=find(wall);walls=struct('faceID',w,'ownerCell',owner(w),'centroid',f.centroid(w,:), ...
    'areaPerRadian',f.areaPerRadian(w),'normal',f.normal(w,:),'temperatureInput',wallInput(w,3), ...
    'velocityInput',wallInput(w,1:2),'pressure',left(w,5),'heatFluxIntoWall',heat(w), ...
    'viscousTraction',traction(w,:),'mechanicalPowerIntoFluidPerRadian',sum(wallInput(w,1:2).*traction(w,:),2).*f.areaPerRadian(w), ...
    'totalEnergyFluxOutPerRadian',integrated(w,4));
info=struct('primitive',P,'faceFluxIntegralPerRadian',integrated,'netBoundaryFlux',net, ...
    'sourceIntegral',src,'cellSourceIntegral',source,'conservationDefect',defect, ...
    'wall',walls,'convectiveFrequency',convective,'diffusiveFrequency',diffusive, ...
    'cellVelocityGradientX',gx(:,1:2),'cellVelocityGradientR',gr(:,1:2),'sourceHoopStress',ttt, ...
    'cellPressureAtPlanarCentroid',pA,'faceBoundaryKind',kind, ...
    'method','Axisymmetric per-radian FV; first-order HLLC; WLS/nonorthogonal Newton+Fourier; pressure/hoop source', ...
    'temperatureModel','calorically perfect; prescribed constant viscosity and Prandtl', ...
    'wallTemperatureValueKind','prescribed','physicalSolution',false,'eligibleForSurfaceCoverage',false);
if wallVelocityRequested&&~isempty(wallVelocityGradient.faceID)
    info.wallVelocityGradient=wallVelocityGradient;
    info.method=[info.method '; explicit static-wall Dirichlet-tangent Newton velocity gradient; no state clipping'];
end
if pressureQuadratureRequested
    info.pressureSourceQuadrature=pressureSourceQuadrature;
    info.method=[info.method '; explicit pressure-source quadrature=' pressureSourceQuadrature.policy ...
        '; unchanged cell state and pressure flux; no post-clipping'];
end
if interiorFluxRequested
    info.interiorConvectiveFlux=interiorConvectiveFlux;
    if ~strcmp(interiorPolicy,'hllc')
        info.method=strrep(info.method,'first-order HLLC',['first-order HLLE (' interiorPolicy ')']);
    end
    info.method=[info.method '; explicit interior Euler flux=' interiorPolicy ...
        '; mass-consistent SA advection; unchanged open-boundary flux; no state clipping'];
end
if arithmeticRequested
    info.interiorFluxArithmetic=arithmetic;
    info.method=[info.method '; explicit interior Euler arithmetic=' arithmetic '; viscous/boundary/source arithmetic unchanged'];
end
if viscousWorkVelocityRequested&&strcmp(viscousWorkVelocity.policy,'arithmetic')
    info.viscousWorkVelocity=viscousWorkVelocity;
    info.method=[info.method '; explicit first-order arithmetic interior viscous work velocity; no state clipping'];
end
if wallThermalRequested
    info.wallThermalGradient=wallThermalGradient;
    if ~isempty(wallThermalGradient.faceID)
        info.method=[info.method '; explicit Dirichlet-tangent Fourier wall gradient; signed q, no clipping or heat multiplier'];
    end
end
if variableTransport
    info.temperatureModel='calorically perfect; local Sutherland laminar viscosity; k=cp*mu(T)/constant Pr';
    info.method=[info.method '; local Sutherland laminar transport; cellwise Barth-Jespersen bounded constitutive face/source temperatures'];
    info.transport=struct('model',gas.transport,'faceTemperature',facePhi(:,3), ...
        'faceMu',faceTransport.mu,'faceK',faceTransport.k,'sourceTemperature',sourceTemperature, ...
        'sourceMu',sourceTransport.mu,'cellMu',cellTransport.mu,'cellK',cellTransport.k, ...
        'frequencyFaceMu',frequencyMu,'frequencyFaceK',frequencyK,'reconstruction',transportTrace);
    info.wall.viscosity=faceTransport.mu(w);info.wall.thermalConductivity=faceTransport.k(w);
end
if turbulent
    info.sa=sa;
    info.method=[info.method '; SA-noft2-nonnegative-AJS2012-v1 effective stress/heat; Prt=.90'];
    info.temperatureModel=[info.temperatureModel '; SA eddy viscosity and turbulent conductivity in interior; exact molecular wall k(Tw)'];
end
if any(riemannFace)
    info.method=[info.method '; open-face HLLE reservoir flux'];
    info.boundaryRiemannSolver='HLLE';
    info.boundaryConvectiveMethod='HLLE two-wave numerical flux of cell and prescribed exterior reservoir; positive intermediate density/internal energy checked; not physical flux of an assumed star state';
    info.exteriorReservoirValueKind='prescribed-boundary-input-not-solved-external-environment';
    info.faceConvectiveFluxPerRadian=convectiveFlux.*f.areaPerRadian;
    info.faceViscousFluxPerRadian=viscous.*f.areaPerRadian;
    info.faceWaveSpeed=speed;info.faceIncoming=incoming;info.faceNormalZeroGradient=normalZero;
    info.faceViscousGradientTrace=right;info.faceRiemannExteriorState=riemannState;
end

if axisSymmetry.active
    af=axisSymmetry.faceID;ao=axisSymmetry.ownerCell;ai=axisSymmetry.boundaryIndex;
    info.axisSymmetry=struct('kind',axisSymmetry.kind,'policy',axisSymmetry.policy, ...
        'faceID',af,'boundaryIndex',ai,'ownerCell',ao,'physicalPointXR',axisSymmetry.physicalPointXR, ...
        'ghostPointXR',axisSymmetry.ghostPointXR,'ghostDisplacementXR',axisSymmetry.ghostDisplacementXR, ...
        'physicalDisplacementXR',axisSymmetry.physicalDisplacementXR, ...
        'ghostPrimitiveRhoUxUrT',[P(ao,1:2),-P(ao,3),P(ao,4)], ...
        'ghostVelocityTemperaturePressure',ghostPhi,'coefficientPrimitiveRhoUxUrT',trace(ai,1:4), ...
        'coefficientPressurePa',trace(ai,5), ...
        'coefficientStateKind','owner-even-projection-at-zero-measure-axis-not-a-solved-face-trace', ...
        'coefficientStateIsPhysicalFaceSolution',false,'physicalTraceReconstructed',false, ...
        'rawTangentialPrimitiveRhoUxUrT',rawAxisPrimitive, ...
        'rawTangentialPressurePa',rawAxisPrimitive(:,1).*gas.R.*rawAxisPrimitive(:,4), ...
        'rawTangentialStateAdmissible',all(isfinite(rawAxisPrimitive),2)&rawAxisPrimitive(:,1)>0& ...
            rawAxisPrimitive(:,4)>=gas.temperatureRange(1)&rawAxisPrimitive(:,4)<=gas.temperatureRange(2), ...
        'rawTangentialStateUsed',false,'cellGradientXAtAxisOwner',gx(ao,:), ...
        'cellGradientRAtAxisOwner',gr(ao,:),'cellDensityGradientXAtAxisOwner',densityGX(ao), ...
        'cellDensityGradientRAtAxisOwner',densityGR(ao),'faceGradientX',Gx(af,:), ...
        'faceGradientR',Gr(af,:),'hoopLimit',hoop(af),'cellStencilUnchanged',true, ...
        'physicalFaceGeometryChanged',false,'zeroRadiusRegularizedWithEpsilon',false, ...
        'stateClipped',false,'radialDerivativeModel','urOwner/rOwner first-order regular proxy; not high-order for x*r', ...
        'constitutiveTemperatureLimiterPolicy','owner T repeats cell bound; raw extrapolation excluded; all actual face targets retained in cellwise limiter, so relative to extrapolated bounds other face/source mu/k may change', ...
        'physicalSolution',false,'eligibleForSurfaceCoverage',false);
    info.method=[info.method '; opt-in mirror-parity axis WLS ghosts; true zero-area axis face; EOS-consistent owner-projection coefficient proxy; raw tangential extrapolation diagnostic only'];
end

    function F=physical(U,p,n)
        un=sum(p(:,2:3).*n,2);F=[p(:,1).*un,U(:,2:3).*un+p(:,5).*n,(U(:,4)+p(:,5)).*un];
    end
    function H=hlleBoundary(UL,UR,l,r,n)
        % 两波HLL/HLLE开放通量；统一使用该格式，不依星压在两种格式之间跳变。
        % 波速包围两侧声学特征；只检查实际取用的中间态，绝不裁rho或内能。
        normalL=sum(l(:,2:3).*n,2);normalR=sum(r(:,2:3).*n,2);
        low=min(normalL-l(:,6),normalR-r(:,6));high=max(normalL+l(:,6),normalR+r(:,6));
        FL=physical(UL,l,n);FR=physical(UR,r,n);H=FL;
        rightward=high<=0;if any(rightward),H(rightward,:)=FR(rightward,:);end
        middle=low<0&high>0;
        if any(middle)
            lo=low(middle);hi=high(middle);den=hi-lo;
            average=(hi.*UR(middle,:)-lo.*UL(middle,:)-(FR(middle,:)-FL(middle,:)))./den;
            energy=average(:,4)-sum(average(:,2:3).^2,2)./(2*average(:,1));
            assert(all(isfinite(average),'all')&&all(average(:,1)>0&energy>0), ...
                'AxisFV:HLLE','HLLE实际中间态非有限或非正，没有裁密度或内能。');
            H(middle,:)=(hi.*FL(middle,:)-lo.*FR(middle,:)+ ...
                lo.*hi.*(UR(middle,:)-UL(middle,:)))./den;
        end
        assert(all(isfinite(H),'all'),'AxisFV:HLLE','HLLE开放通量非有限。');
    end
    function F=hllc(UL,UR,l,r,n)
        unL=sum(l(:,2:3).*n,2);unR=sum(r(:,2:3).*n,2);
        SL=min(unL-l(:,6),unR-r(:,6));SR=max(unL+l(:,6),unR+r(:,6));
        % 全单向波区只采样原状态的物理通量；未取用的星态不属于验收对象。
        % 在真实Mach7显式Q1中，SL>0而未使用p*略负；旧无条件检查会误拒绝。
        assert(all(isfinite([SL,SR]),'all')&&all(SL<SR),'AxisFV:HLLC','HLLC端点波速无效。');
        FL=physical(UL,l,n);FR=physical(UR,r,n);F=FL;
        maskFarRight=SR<=0;
        if any(maskFarRight),F(maskFarRight,:)=FR(maskFarRight,:);end
        mixedMask=SL<0&SR>0;
        if ~any(mixedMask),return;end
        % 只有SL<0<SR的混合波区构造双星态，并保留全部正性/波序保护。
        UL=UL(mixedMask,:);UR=UR(mixedMask,:);l=l(mixedMask,:);r=r(mixedMask,:);n=n(mixedMask,:);
        unL=unL(mixedMask);unR=unR(mixedMask);SL=SL(mixedMask);SR=SR(mixedMask);
        dL=l(:,1).*(SL-unL);dR=r(:,1).*(SR-unR);
        SM=(r(:,5)-l(:,5)+dL.*unL-dR.*unR)./(dL-dR);
        starL=l(:,1).*(SL-unL)./(SL-SM);starR=r(:,1).*(SR-unR)./(SR-SM);
        pL=l(:,5)+dL.*(SM-unL);pR=r(:,5)+dR.*(SM-unR);
        USL=[starL,starL.*(l(:,2:3)+(SM-unL).*n), ...
            starL.*(UL(:,4)./l(:,1)+(SM-unL).*(SM+l(:,5)./dL))];
        USR=[starR,starR.*(r(:,2:3)+(SM-unR).*n), ...
            starR.*(UR(:,4)./r(:,1)+(SM-unR).*(SM+r(:,5)./dR))];
        eL=USL(:,4)-sum(USL(:,2:3).^2,2)./(2*starL);
        eR=USR(:,4)-sum(USR(:,2:3).^2,2)./(2*starR);
        pressureScale=max([abs(l(:,5)),abs(r(:,5)),abs(dL.*(SM-unL)),abs(dR.*(SM-unR))],[],2);
        assert(all(isfinite([USL,USR,SM]),'all')&&all(starL>0&starR>0&pL>0&pR>0), ...
            'AxisFV:HLLC','HLLC星态非有限/非正，不以裁值代替。');
        assert(all(SM>SL&SM<SR&eL>0&eR>0)&&all(abs(pL-pR)<=256*eps(pressureScale)), ...
            'AxisFV:HLLC','HLLC波序、实际星内能或左右星压力舍入合同不满足。');
        FL=FL(mixedMask,:);FR=FR(mixedMask,:);mixedFlux=FL;
        % 避免与外层边界循环b共享嵌套工作区，单面调用也不得改写外层下标。
        maskLeft=SM>=0;maskRight=SM<0;
        % 单个开放面的逻辑空索引可能为1x0；先判分支，保持标量/批量同一公式。
        if any(maskLeft),mixedFlux(maskLeft,:)=FL(maskLeft,:)+SL(maskLeft).*(USL(maskLeft,:)-UL(maskLeft,:));end
        if any(maskRight),mixedFlux(maskRight,:)=FR(maskRight,:)+SR(maskRight).*(USR(maskRight,:)-UR(maskRight,:));end
        F(mixedMask,:)=mixedFlux;
    end
end
