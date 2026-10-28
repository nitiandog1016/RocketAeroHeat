function [dQdt,info]=compressibleNavierStokesResidual3D(Q,grid,solid,gas,boundary)
%COMPRESSIBLENAVIERSTOKESRESIDUAL3D 三维可压缩黏性/导热实验数值库空间算子。
% 模块作用：在均匀笛卡尔结构网格上，以共享面通量计算dQ/dt，Q最后
% 一维依次为rho,rho*u,rho*v,rho*w,rho*E。全部采用SI，不推进时间。
% 气体：定gamma理想气体，常数mu>0、k>0，cp=gamma*R/(gamma-1)，
% k=mu*cp/Pr。输入temperatureRange为本次明确接受的模型温度范围；
% 不是NASA7或真实空气普适范围。无湍流、滑移、化学、辐射或烧蚀。
% 对流：Toro HLLC、Davis左右声速界；黏性：Newton应力/Stokes假设；
% 导热：Fourier定律。固体为阶梯壁，墙面直接施加无滑移给定Tw，
% 不构造负温度ghost，也不裁剪rho/T。不可用输入/中间态明确报错。
% 来源：NASA/TM-2014-218282式1--10；Toro2010 HLLC原作者讲义；
% TN-D7508印刷4页式1--5提供底区需要黏性/能量方程的物理依据，
% 本函数不是该报告的二维Telenin求解器，不能继承其精度。
%
% grid.spacing=[dx dy dz] [m]；solid同前三维logical，固体Q可为NaN。
% gas.gamma/R/mu/Pr/temperatureRange必填，gas.k可选但须与Pr一致。
% boundary.xmin/xmax/ymin/ymax/zmin/zmax必须明确给出type：
% periodic成对；supersonicInflow另给primitive=[rho u v w T]；
% supersonicOutflow要求当前流动全特征向外，法向原始量梯度为零；
% characteristic须给mode=riemannFarfield或pressureOutlet及明确储库；
% 后者另给pressure。亚声速回流只用储库入射特征，不强制正速度。
% wall另给temperature、velocity=[uw vw ww]，法向墙速须零。
% 有固体时boundary.solidWall同wall合同。给定Tw始终不是求解壁温。
% 输出dQdt固体内NaN；info含共享面通量/墙热流/力/净边界积分；
% eligibleForSurfaceCoverage=false，不代表迭代/网格/边界收敛的物理解。

%% 1. 检查几何、物性和边界输入，拒绝任何无来源隐式边界
names={'xmin','xmax';'ymin','ymax';'zmin','zmax'};
spacing=double(grid.spacing(:)');validateattributes(spacing,{'double'},{'finite','positive','numel',3});
shape=[size(Q,1),size(Q,2),size(Q,3)];n=prod(shape);
assert(size(Q,4)==5 && ndims(Q)<=4,'RocketAeroHeat:NSDraft:Shape','Q须为Nx*Ny*Nz*5。');
assert(islogical(solid) && isequal([size(solid,1),size(solid,2),size(solid,3)],shape), ...
    'RocketAeroHeat:NSDraft:Mask','solid须为同网格logical掩码。');
solid=solid(:);fluid=~solid;assert(any(fluid),'RocketAeroHeat:NSDraft:NoFluid','没有流体控制体。');
for key={'gamma','R','mu','Pr'}
    validateattributes(gas.(key{1}),{'numeric'},{'real','finite','scalar','positive'});
end
assert(gas.gamma>1,'RocketAeroHeat:NSDraft:Gamma','gamma须大于1。');
range=double(gas.temperatureRange(:)');
validateattributes(range,{'double'},{'finite','positive','numel',2});
assert(range(2)>range(1),'RocketAeroHeat:NSDraft:Range','温度接受范围须递增。');
cp=gas.gamma*gas.R/(gas.gamma-1);cv=gas.R/(gas.gamma-1);
k=gas.mu*cp/gas.Pr;
if isfield(gas,'k')
    validateattributes(gas.k,{'numeric'},{'real','finite','scalar','positive'});
    assert(abs(gas.k/k-1)<1e-12,'RocketAeroHeat:NSDraft:Pr','k与mu*cp/Pr不一致。');k=gas.k;
end
bc=cell(3,2);
for axis=1:3
    for side=1:2
        assert(isfield(boundary,names{axis,side}),'RocketAeroHeat:NSDraft:Boundary', ...
            '须明确给出%s边界。',names{axis,side});
        bc{axis,side}=checkBoundary(boundary.(names{axis,side}),axis,2*side-3);
    end
    periodic=strcmp(bc{axis,1}.type,'periodic');
    assert(periodic==strcmp(bc{axis,2}.type,'periodic'),'RocketAeroHeat:NSDraft:Periodic', ...
        '周期边界须成对。');
    if periodic
        sm=reshape(solid,shape);a={':',':',':'};b=a;a{axis}=1;b{axis}=shape(axis);
        assert(isequal(sm(a{:}),sm(b{:})),'RocketAeroHeat:NSDraft:PeriodicSolid', ...
            '本草稿要求周期端固体掩码对应一致。');
    end
end
if any(solid)
    assert(isfield(boundary,'solidWall'),'RocketAeroHeat:NSDraft:SolidWall','有固体时须给solidWall边界。');
    solidWall=checkBoundary(boundary.solidWall,0,0);
    assert(strcmp(solidWall.type,'wall'),'RocketAeroHeat:NSDraft:SolidWall','solidWall只能是wall。');
end
flat=reshape(double(Q),n,5);P=nan(n,8);
for id=find(fluid)'
    P(id,:)=primitive(flat(id,:));
end

%% 2. 流体中心的速度/温度梯度；壁值直接在半单元距离取样
grad=zeros(n,4,3);
for id=find(fluid)'
    [ii,jj,kk]=ind2sub(shape,id);sub=[ii,jj,kk];
    for axis=1:3
        [minus,dm]=sample(id,sub,axis,-1);[plus,dp]=sample(id,sub,axis,1);
        grad(id,:,axis)=(plus-minus)/(dm+dp);
    end
end

%% 3. 每个共享面只求一次通量，再以相反符号更新相邻流体
rate=zeros(n,5);rate(solid,:)=NaN;faceFlux=cell(1,3);faceKind=cell(1,3);
walls=struct('fluidCell',{},'axis',{},'normal',{},'area',{},'temperatureInput',{}, ...
    'velocityInput',{},'pressure',{},'traction',{},'heatFluxIntoWall',{},'totalEnergyFluxOut',{});
openFaces=struct('fluidCell',{},'axis',{},'normal',{},'area',{},'primitive',{}, ...
    'boundaryState',{},'fluxOut',{});
net=zeros(1,5);wallNet=zeros(1,5);openNet=zeros(1,5);fallbacks=0;
for axis=1:3
    fs=shape;fs(axis)=fs(axis)+1;count=prod(fs);ff=zeros(count,5);kind=zeros(count,1,'uint8');
    faceArea=prod(spacing)/spacing(axis);periodic=strcmp(bc{axis,1}.type,'periodic');
    for face=1:count
        [ii,jj,kk]=ind2sub(fs,face);sub=[ii,jj,kk];coord=sub(axis);
        if periodic && coord==shape(axis)+1,continue;end
        left=0;right=0;external=coord==1 || coord==shape(axis)+1;
        if coord>1,leftSub=sub;leftSub(axis)=coord-1;left=sub2ind(shape,leftSub(1),leftSub(2),leftSub(3));end
        if coord<=shape(axis),rightSub=sub;rightSub(axis)=coord;right=sub2ind(shape,rightSub(1),rightSub(2),rightSub(3));end
        if periodic && coord==1
            ls=sub;ls(axis)=shape(axis);left=sub2ind(shape,ls(1),ls(2),ls(3));external=false;
        end
        lf=left>0 && fluid(left);rf=right>0 && fluid(right);
        if ~lf && ~rf,continue;end
        if lf && rf
            F=riemann(flat(left,:),flat(right,:),P(left,:),P(right,:),axis);
            G=0.5*(reshape(grad(left,:,:),4,3)+reshape(grad(right,:,:),4,3));
            G(:,axis)=(P(right,2:5)-P(left,2:5))'/spacing(axis);
            velocity=.5*(P(left,2:4)+P(right,2:4));
            F=F-viscous(G,velocity,axis);
            kind(face)=1+uint8(periodic && coord==1);
        else
            if lf,id=left;sign=1;else,id=right;sign=-1;end
            if external
                side=1+(coord>1);b=bc{axis,side};
            else
                b=solidWall;
            end
            if strcmp(b.type,'wall')
                assert(abs(b.velocity(axis))<1e-14,'RocketAeroHeat:NSDraft:WallMotion', ...
                    '阶梯壁法向速度必须为零。');
                normal=zeros(1,3);normal(axis)=sign;
                G=zeros(4,3);G(:,axis)=sign*([b.velocity,b.temperature]-P(id,2:5))'/(spacing(axis)/2);
                visc=viscous(G,b.velocity,axis)*sign;
                outward=[0,P(id,6)*normal,0]-visc;
                F=sign*outward;kind(face)=3;
                walls(end+1)=struct('fluidCell',id,'axis',axis,'normal',normal,'area',faceArea, ...
                    'temperatureInput',b.temperature,'velocityInput',b.velocity,'pressure',P(id,6), ...
                    'traction',outward(2:4),'heatFluxIntoWall',-sign*k*G(4,axis), ...
                    'totalEnergyFluxOut',outward(5)); %#ok<AGROW>
                wallNet=wallNet+faceArea*outward;
            elseif strcmp(b.type,'supersonicInflow')
                given=b.primitive;pg=primitive(conservative(given));
                assert(sign*pg(axis+1)<-pg(7),'RocketAeroHeat:NSDraft:InflowMach', ...
                    'supersonicInflow须全部特征向内；不支持亚声速或背流替代。');
                G=reshape(grad(id,:,:),4,3);
                G(:,axis)=sign*(given(2:5)-P(id,2:5))'/(spacing(axis)/2);
                F=physicalFlux(conservative(given),pg,axis)-viscous(G,given(2:4),axis);
                kind(face)=4;openNet=openNet+sign*faceArea*F;
            elseif strcmp(b.type,'supersonicOutflow')
                assert(sign*P(id,axis+1)>P(id,7),'RocketAeroHeat:NSDraft:OutflowMach', ...
                    'supersonicOutflow当前不满足全特征向外；亚声速/回流须另行建边界。');
                G=reshape(grad(id,:,:),4,3);G(:,axis)=0;
                F=physicalFlux(flat(id,:),P(id,:),axis)-viscous(G,P(id,2:4),axis);
                kind(face)=5;openNet=openNet+sign*faceArea*F;
            elseif strcmp(b.type,'characteristic')
                normal=zeros(1,3);normal(axis)=sign;
                [given,characteristic]=resolveNSCharacteristicBoundary3D(P(id,1:5),b,normal,gas);
                pg=primitive(conservative(given));G=reshape(grad(id,:,:),4,3);
                if characteristic.isInflow
                    G(:,axis)=sign*(given(2:5)-P(id,2:5))'/(spacing(axis)/2);
                else
                    G(:,axis)=0;
                end
                F=physicalFlux(conservative(given),pg,axis)-viscous(G,given(2:4),axis);
                kind(face)=6;openNet=openNet+sign*faceArea*F;
                openFaces(end+1)=struct('fluidCell',id,'axis',axis,'normal',normal, ...
                    'area',faceArea,'primitive',given,'boundaryState',characteristic, ...
                    'fluxOut',sign*F); %#ok<AGROW>
            else
                error('RocketAeroHeat:NSDraft:Boundary','周期边界内部不应出现孤立流体。');
            end
            net=net+sign*faceArea*F;
        end
        ff(face,:)=F;
        if lf,rate(left,:)=rate(left,:)-F*faceArea/prod(spacing);end
        if rf,rate(right,:)=rate(right,:)+F*faceArea/prod(spacing);end
    end
    if periodic
        ff=reshape(ff,[fs,5]);kind=reshape(kind,fs);
        a={':',':',':',':'};b=a;a{axis}=shape(axis)+1;b{axis}=1;ff(a{:})=ff(b{:});
        aa={':',':',':'};bb=aa;aa{axis}=shape(axis)+1;bb{axis}=1;kind(aa{:})=kind(bb{:});
    else
        ff=reshape(ff,[fs,5]);kind=reshape(kind,fs);
    end
    faceFlux{axis}=ff;faceKind{axis}=kind;
end
dQdt=reshape(rate,[shape,5]);
integrated=sum(rate(fluid,:),1)*prod(spacing);
info=struct('status','experimental-spatial-operator-only','physicalSolution',false, ...
    'eligibleForSurfaceCoverage',false,'timeIntegratorImplemented',false,'gridConvergenceVerified',false, ...
    'experimentalValidation',false,'method','HLLC对流 + Newton黏性应力 + Fourier导热，共享面有限体积', ...
    'scope','常物性定比热层流NS空间残差；阶梯solid壁；未推进或收敛，不提供完整底热/三维外流解', ...
    'temperatureRange',range,'mu',gas.mu,'k',k,'Pr',gas.Pr,'cp',cp,'cv',cv, ...
    'faceFlux',{faceFlux},'faceKind',{faceKind},'wallFaces',walls, ...
    'primitive',reshape(P(:,1:5),[shape,5]),'netBoundaryFlux',net,'wallBoundaryFlux',wallNet, ...
    'openBoundaryFlux',openNet,'openFaces',openFaces, ...
    'volumeIntegratedResidual',integrated,'balanceResidual',integrated+net, ...
    'fluxUnits',{{'kg/(m^2*s)','Pa','Pa','Pa','W/m^2'}}, ...
    'residualUnits',{{'kg/(m^3*s)','N/m^3','N/m^3','N/m^3','W/m^3'}}, ...
    'boundaryIntegralUnits',{{'kg/s','N','N','N','W'}},'positivityClippingCount',0,'fluxFallbacks',fallbacks);

    function b=checkBoundary(b,axis,sign)
        b.type=validatestring(b.type,{'periodic','supersonicInflow','supersonicOutflow','wall','characteristic'});
        if strcmp(b.type,'wall')
            validateTemperature(b.temperature);validateattributes(b.velocity,{'numeric'},{'real','finite','numel',3});
            b.velocity=double(b.velocity(:)');
            if axis>0,assert(abs(b.velocity(axis))<1e-14,'RocketAeroHeat:NSDraft:WallMotion','外墙法向墙速须零。');end
        elseif strcmp(b.type,'supersonicInflow')
            validateattributes(b.primitive,{'numeric'},{'real','finite','numel',5});b.primitive=double(b.primitive(:)');
            pg=primitive(conservative(b.primitive));
            assert(sign*pg(axis+1)<-pg(7),'RocketAeroHeat:NSDraft:InflowMach','输入不是向内超声速。');
        elseif strcmp(b.type,'characteristic')
            assert(axis>0 && isfield(b,'mode') && isfield(b,'reservoirPrimitive'), ...
                'RocketAeroHeat:NSDraft:CharacteristicInput','特征外边界須明确mode/储库。');
            b.mode=validatestring(b.mode,{'riemannFarfield','pressureOutlet'});
            validateattributes(b.reservoirPrimitive,{'numeric'},{'real','finite','numel',5});
            b.reservoirPrimitive=double(b.reservoirPrimitive(:)');
            primitive(conservative(b.reservoirPrimitive));
            if strcmp(b.mode,'pressureOutlet')
                assert(isfield(b,'pressure'),'RocketAeroHeat:NSDraft:CharacteristicInput','须给出口压力。');
                validateattributes(b.pressure,{'numeric'},{'real','finite','scalar','positive'});
            end
        end
    end
    function validateTemperature(T)
        assert(isscalar(T) && isfinite(T) && isreal(T) && T>=range(1) && T<=range(2), ...
            'RocketAeroHeat:NSDraft:Temperature','静温/壁温/中间态超出显式温度范围；不裁温。');
    end
    function p=primitive(q)
        assert(all(isfinite(q)) && isreal(q) && q(1)>0,'RocketAeroHeat:NSDraft:State','流体保守量无效；不裁密度。');
        rho=q(1);u=q(2:4)/rho;internal=q(5)-.5*rho*sum(u.^2);
        assert(internal>0,'RocketAeroHeat:NSDraft:InternalEnergy','内能非正；不能以裁能量假造解。');
        T=internal/(rho*cv);validateTemperature(T);p=[rho,u,T,rho*gas.R*T,sqrt(gas.gamma*gas.R*T),internal];
    end
    function q=conservative(p)
        assert(p(1)>0,'RocketAeroHeat:NSDraft:State','边界密度必须正。');validateTemperature(p(5));
        q=[p(1),p(1)*p(2:4),p(1)*(cv*p(5)+.5*sum(p(2:4).^2))];
    end
    function [value,distance]=sample(id,sub,axis,sign)
        neighbor=sub;neighbor(axis)=neighbor(axis)+sign;inside=neighbor(axis)>=1 && neighbor(axis)<=shape(axis);
        if inside
            other=sub2ind(shape,neighbor(1),neighbor(2),neighbor(3));
            if fluid(other),value=P(other,2:5);distance=spacing(axis);return;else,b=solidWall;end
        else
            side=1+(sign>0);b=bc{axis,side};
            if strcmp(b.type,'periodic')
                neighbor(axis)=1+(sign<0)*(shape(axis)-1);other=sub2ind(shape,neighbor(1),neighbor(2),neighbor(3));
                value=P(other,2:5);distance=spacing(axis);return;
            end
        end
        if strcmp(b.type,'wall'),value=[b.velocity,b.temperature];distance=spacing(axis)/2;
        elseif strcmp(b.type,'supersonicInflow'),value=b.primitive(2:5);distance=spacing(axis)/2;
        elseif strcmp(b.type,'supersonicOutflow'),value=P(id,2:5);distance=spacing(axis);
        elseif strcmp(b.type,'characteristic')
            sampleNormal=zeros(1,3);sampleNormal(axis)=sign;
            [faceState,characteristic]=resolveNSCharacteristicBoundary3D(P(id,1:5),b,sampleNormal,gas);
            if characteristic.isInflow,value=faceState(2:5);distance=spacing(axis)/2;
            else,value=P(id,2:5);distance=spacing(axis);end
        else,error('RocketAeroHeat:NSDraft:Sample','未支持的边界取样。');end
    end
    function f=physicalFlux(q,p,axis)
        % 独立局部名避免MATLAB嵌套函数把外层实际边界法向改成+axis。
        u=p(axis+1);fluxNormal=zeros(1,3);fluxNormal(axis)=1;
        f=[q(1)*u,q(2:4)*u+p(6)*fluxNormal,(q(5)+p(6))*u];
    end
    function f=viscous(G,u,axis)
        velocityGradient=G(1:3,:);tau=gas.mu*(velocityGradient+velocityGradient' ...
            -(2/3)*trace(velocityGradient)*eye(3));traction=tau(:,axis)';
        f=[0,traction,dot(u,traction)+k*G(4,axis)];
    end
    function f=riemann(qL,qR,pL,pR,axis)
        % Toro HLLC三波合同；Davis波速界不构成多维正性/熵稳定保证。
        uL=pL(axis+1);uR=pR(axis+1);SL=min(uL-pL(7),uR-pR(7));SR=max(uL+pL(7),uR+pR(7));
        fL=physicalFlux(qL,pL,axis);fR=physicalFlux(qR,pR,axis);
        if SL>=0,f=fL;return;elseif SR<=0,f=fR;return;end
        denominator=pL(1)*(SL-uL)-pR(1)*(SR-uR);
        assert(isfinite(denominator) && denominator~=0,'RocketAeroHeat:NSDraft:HLLC','HLLC接触波分母退化。');
        SM=(pR(6)-pL(6)+pL(1)*uL*(SL-uL)-pR(1)*uR*(SR-uR))/denominator;
        assert(all(isfinite([SL,SM,SR])) && SM>SL && SM<SR, ...
            'RocketAeroHeat:NSDraft:HLLC','HLLC波序失效；未自动回退。');
        % 模块作用：Q*的正内能不保证HLLC动量通量中的星压力为正。
        % 左右Rankine-Hugoniot动量关系须分别给出有限正p*且舍入内一致。
        % 强膨胀反例可有rho*>0、T*>0却p*<0；明确拒绝，不裁压力回退。
        jumpL=pL(1)*(SL-uL)*(SM-uL);jumpR=pR(1)*(SR-uR)*(SM-uR);
        pressureStarL=pL(6)+jumpL;pressureStarR=pR(6)+jumpR;
        assert(all(isfinite([pressureStarL,pressureStarR])) && ...
            pressureStarL>0 && pressureStarR>0,'RocketAeroHeat:NSDraft:HLLCStarPressure', ...
            'HLLC星压力非正或非有限；拒绝强膨胀通量，不以Q*正内能替代压力检查。');
        pressureScale=max([1,abs(pL(6)),abs(pR(6)),abs(jumpL),abs(jumpR)]);
        assert(abs(pressureStarL-pressureStarR)<=128*eps(pressureScale), ...
            'RocketAeroHeat:NSDraft:HLLCStarPressureMismatch', ...
            '左右HLLC星压力超出浮点舍入一致性范围。');
        if SM>=0,qs=star(qL,pL,SL,SM,uL,axis);f=fL+SL*(qs-qL);
        else,qs=star(qR,pR,SR,SM,uR,axis);f=fR+SR*(qs-qR);end
        primitive(qs); % 检查实际中间内能/温度，不裁剪或隐藏不适用。
    end
    function qs=star(q,p,S,SM,u,axis)
        rs=p(1)*(S-u)/(S-SM);velocity=p(2:4);velocity(axis)=SM;
        specific=q(5)/p(1)+(SM-u)*(SM+p(6)/(p(1)*(S-u)));
        qs=[rs,rs*velocity,rs*specific];
    end
end
