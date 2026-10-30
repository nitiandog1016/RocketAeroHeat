function column=solveParabolicShockLayerColumn(upstream,gas,shape,q,eta,input)
%SOLVEPARABOLICSHOCKLAYERCOLUMN 共同抛物激波的质量坐标、熵输运和薄层反算。
% 模块一：由同一曲面求跃迁/曲率，正向质量坐标直接由均匀上游流量建立。
% 模块二：沿同一激波线回溯入口，使用其真实激波熵，不逐站重置总压。
% 模块三：TP-3018显式压力/法向速度近似，同cp恢复T/rho/切向速度。
% 模块四：保留非正交Jacobian，将质量积分反解为激波内法向距离。
% 范围：共同正定抛物激波、均匀上游、薄激波层；尚非完整区域外流解。
if nargin<6,input=struct;end
order=16;if isfield(input,'quadratureOrder'),order=input.quadratureOrder;end
validateattributes(order,{'double'},{'real','finite','scalar','integer','>=',4,'<=',64});
validateattributes(eta,{'double'},{'real','finite','vector','>=',0,'<=',1});eta=eta(:);
assert(numel(eta)>=2&&eta(1)==1&&eta(end)==0&&all(diff(eta)<0), ...
    'ShockLayer:MassCoordinate','质量比从激波1严格递减到体表0。');
sample=parabolicShockSample(shape,q);shock=shockSurfaceJumpState(upstream,gas,sample,1);
velocity=upstream.velocityVector;U=norm(velocity);ex=velocity/U;
assert(norm(shape.frame(:,1)-ex)<1e-11,'ShockLayer:Frame','激波迎风轴须沿实际来流速度。');
tangent=velocity-dot(velocity,shock.normal)*shock.normal;us=norm(tangent);
assert(us>1e-6*U,'ShockLayer:StagnationLimit','驻点极限须单独整体处理，不能以有限列代替。');
tangent=tangent/us;transverse=cross(shock.normal,tangent);
A=sample.firstDerivatives;dn=shock.partialDerivatives.normal;
coordinates=A\[tangent transverse];S=-[tangent transverse]'*dn*coordinates;
assert(norm(S-S','fro')<1e-10*max(1,norm(S,'fro')),'ShockLayer:Curvature','共同曲面曲率不对称。');
ks=S(1,1);kb=S(2,2);twist=(S(1,2)+S(2,1))/2;largest=max(eig(S));
assert(largest>0,'ShockLayer:Curvature','本分支只处理正定抛物激波。');
vs=dot(shock.velocityVector,shock.normal);sine=dot(ex,shock.normal);cosine=us/U;
g=shape.K*q;massPerWidth=upstream.density*U*norm(g)/trace(shape.K);
% SI形式：psi_s/h_beta= rho_inf*U_inf*|Kq|/tr(K)。速度式保留rho_inf*U_inf。
P1=us*massPerWidth*ks;
P2=.5*vs*massPerWidth*(sine/cosine)*(ks+kb);
C=massPerWidth/(upstream.density*U*cosine)*(ks+kb);
initial=surfaceFlowThermodynamicState(upstream.temperature,upstream.pressure,gas);
H0=initial.enthalpy+U^2/2;
[eigenvectors,eigenvalues]=eig(shape.K);lambda=diag(eigenvalues);power=trace(shape.K)/min(lambda);
state=atEta(eta);

%% 一、不同阶Gauss积分独立记录质量距离的分辨率，不把差值称严格误差界
fine=massIntegrals(2*order);coarse=massIntegrals(order);
zeta=[0;cumsum(fine)];zetaCoarse=[0;cumsum(coarse)];
distances=nan(size(eta));coarseDistances=distances;
for j=1:numel(eta)
    distances(j)=invertCapacity(zeta(j));coarseDistances(j)=invertCapacity(zetaCoarse(j));
end
indicator=abs(distances-coarseDistances);
assert(max(indicator)<2e-6*max(distances(end),1e-6), ...
    'ShockLayer:Quadrature','跨层质量积分未分辨，须加密而非放宽匹配。');

%% 二、真实非正交法线坐标；不以两条独立正交尺度替代Jacobian
a=1-distances*ks;b=1-distances*kb;crossTerm=-distances*twist;
J=a.*b-crossTerm.^2;h=sqrt(a.^2+crossTerm.^2);
assert(all(J>0)&&all(a>0),'ShockLayer:Caustic','法线坐标已聚焦或翻向，不能向体表外推。');
xi=(a.*tangent'+crossTerm.*transverse')./h;
physicalVelocity=state.streamwiseSpeed.*xi+state.normalSpeed.*shock.normal';
positions=shock.position'+distances.*shock.normal';
column=struct('method','共同抛物激波质量坐标 + TP-3018薄层压力/法向速度 + 同cp熵/焓恢复', ...
    'shape',shape,'coordinate',q,'shock',shock,'eta',eta,'state',state, ...
    'normalDistanceM',distances,'position',positions,'velocityVector',physicalVelocity, ...
    'massDistanceM',zeta,'coarseMassDistanceM',zetaCoarse,'coarseNormalDistanceM',coarseDistances, ...
    'quadratureDifferenceM',indicator,'quadratureDifferenceIsRigorousBound',false, ...
    'quadratureOrder',order,'curvatureMatrix',S,'jacobianRatio',J,'streamwiseScaleRatio',h, ...
    'offDiagonalMetricRatio',crossTerm.*(a+b),'streamwiseDirection',xi, ...
    'massPerShockWidthKgPerMS',massPerWidth,'pressureCoefficientsPa',[P1 P2], ...
    'normalVelocityCoefficient',C,'totalEnthalpyReference',H0,'entropySourceMap', ...
    'q_entry=exp(K*log(eta)/tr(K))*q；eta=0取同一曲面驻点激波熵极限', ...
    'massCoordinateMethod','div_q((rho_inf*U_inf/tr(K))*K*q)=rho_inf*U_inf；正向质量来自同一均匀上游', ...
    'normalDistanceSignConvention','n朝激波内侧；eta从1降到0；dn/deta=-psi_s*h_xi/(rho*u*J)', ...
    'stateLocation','由给定共同激波反算出的层内及候选生成体表，尚未匹配真实箭体', ...
    'fullEulerSolved',false,'shockShapeSolved',false,'bodyMatched',false,'boundaryLayerEdgeAvailable',false, ...
    'heatFluxSolved',false,'wallTemperatureSolved',false,'productionQualified',false,'newProductionSolvedAreaM2',0);

    function segments=massIntegrals(n)
        [nodes,weights]=gaussRule(n);targets=eta.^(1/power);num=numel(eta)-1;
        widths=(targets(1:end-1)-targets(2:end))/2;centers=(targets(1:end-1)+targets(2:end))/2;
        t=centers+widths.*nodes';allEta=t.^power;values=atEta(allEta(:));
        integrand=reshape(1./(values.density.*values.streamwiseSpeed),num,n).*power.*t.^(power-1);
        segments=massPerWidth*widths.*(integrand*weights);
    end

    function values=atEta(fraction)
        fraction=fraction(:);count=numel(fraction);sourceQ=zeros(2,count);
        nonzero=fraction>0;
        sourceQ(:,nonzero)=eigenvectors*((eigenvectors'*q).* ...
            exp(lambda*log(fraction(nonzero)')/trace(shape.K)));
        gradient=shape.K*sourceQ;normalSpeed=U./sqrt(1+sum(gradient.^2,1))';
        entry=jumpEntropy(normalSpeed);
        p=shock.pressure+P1*(fraction-1)-P2*(fraction.^2-1);
        v=vs*(1+C*(fraction-1));
        assert(all(p>0),'ShockLayer:Pressure','薄层近似出现非正压力，不能补值。');
        lo=repmat(gas.temperatureRangeK(1),count,1);hi=repmat(gas.temperatureRangeK(2),count,1);
        [~,~,slo]=thermal(lo,p);[~,~,shi]=thermal(hi,p);
        assert(all(entry.entropy>=slo & entry.entropy<=shi),'ShockLayer:Temperature','沿真实入口熵恢复的温度超出物性范围。');
        for k=1:52
            mid=(lo+hi)/2;[~,~,s]=thermal(mid,p);lower=s<entry.entropy;
            lo(lower)=mid(lower);hi(~lower)=mid(~lower);
        end
        T=(lo+hi)/2;[cp,hval,entropy]=thermal(T,p);rho=p./(gas.R*T);squared=2*(H0-hval)-v.^2;
        assert(all(squared>0),'ShockLayer:Energy','层内能量不足以得到正切向速度；不能裁剪负值。');
        u=sqrt(squared);energy=(hval+(u.^2+v.^2)/2-H0)/max(1,abs(H0));
        values=struct('pressure',p,'temperature',T,'density',rho,'streamwiseSpeed',u,'normalSpeed',v, ...
            'speed',sqrt(u.^2+v.^2),'specificHeat',cp,'entropy',entropy,'sourceEntropy',entry.entropy, ...
            'sourceShockCoordinate',sourceQ','sourceNormalMach',entry.normalMach,'sourceJumpEnergyResidual',entry.energyResidual, ...
            'energyResidual',energy,'entropyResidual',entropy-entry.entropy, ...
            'pressureMethod','TP-3018式17的SI形式；薄激波层显式压力近似，非完整Euler动量解', ...
            'temperatureMethod','同一流线来源的激波熵与当前压力恢复，不沿列冻结局部激波熵', ...
            'speedMethod','统一上游总焓扣除静焓和式18法向动能');
    end

    function entry=jumpEntropy(w)
        count=numel(w);T1=upstream.temperature;p1=upstream.pressure;rho1=upstream.density;
        Mach=w/sqrt(initial.gammaLocal*gas.R*T1);
        assert(all(Mach>1+1e-5),'ShockLayer:NormalMach','共同激波线上的来源点不满足压缩激波范围。');
        normalH=initial.enthalpy+w.^2/2;lo=1+min(1e-4,.25*(Mach-1));hi=repmat(1e4,count,1);
        for k=1:72
            r=(lo+hi)/2;temp=T1./r+(w.^2/gas.R).*(r-1)./r.^2;
            E=nan(count,1);above=temp>gas.temperatureRangeK(2);below=temp<gas.temperatureRangeK(1);
            E(above)=1;E(below)=-1;inside=~(above|below);
            [~,hh]=thermal(temp(inside),repmat(p1,nnz(inside),1));
            E(inside)=hh+.5*(w(inside)./r(inside)).^2-normalH(inside);
            positive=E>0;lo(positive)=r(positive);hi(~positive)=r(~positive);
        end
        ratio=(lo+hi)/2;p2=p1+rho1*w.^2.*(1-1./ratio);T2=p2./(rho1*ratio*gas.R);
        [cp,hh,s]=thermal(T2,p2);down=w./ratio./sqrt(cp./(cp-gas.R).*gas.R.*T2);
        residual=(hh+.5*(w./ratio).^2-normalH)./max(1,abs(normalH));
        assert(all(ratio>1 & T2>T1 & p2>p1 & down<1 & s>initial.entropy)&&max(abs(residual))<1e-9, ...
            'ShockLayer:CompressionRoot','同一激波来源的压缩根/熵增/能量未通过。');
        entry=struct('entropy',s,'normalMach',Mach,'energyResidual',residual);
    end

    function [cp,hval,s]=thermal(T,p)
        assert(all(isfinite(T))&&all(T>=gas.temperatureRangeK(1)&T<=gas.temperatureRangeK(2)), ...
            'ShockLayer:Temperature','批量物性不得外推温度。');
        if strcmp(gas.thermalPropertyMode,'caloricallyPerfect')
            cp=repmat(gas.gamma*gas.R/(gas.gamma-1),size(T));hval=cp.*T;
            s=cp.*log(T/gas.referenceTemperatureK)-gas.R*log(p/gas.referencePressurePa);
        else
            assert(strcmp(gas.thermalPropertyMode,'thermallyPerfect')&& ...
                strcmp(gas.thermalPropertyConvention,'continuous-project-cp-v1'),'ShockLayer:Gas','不能混用其他物性约定。');
            [cp,hval,standard]=surfaceFlowAirThermalProperties(T);
            [~,~,reference]=surfaceFlowAirThermalProperties(gas.referenceTemperatureK);
            s=standard-reference-gas.R*log(p/gas.referencePressurePa);
        end
        assert(all(isfinite(cp))&&all(cp>gas.R),'ShockLayer:Temperature','cp或批量物性失效。');
    end

    function n=invertCapacity(target)
        if target==0,n=0;return;end
        lo=0;hi=(1-1e-8)/largest;
        assert(capacity(hi)>target,'ShockLayer:Caustic','质量要求超过法线坐标第一聚焦点。');
        for k=1:55
            mid=(lo+hi)/2;if capacity(mid)<target,lo=mid;else,hi=mid;end
        end
        n=(lo+hi)/2;
    end
    function value=capacity(n)
        [x,w]=gaussRule(32);z=(x+1)*(n/2);a=1-z*ks;b=1-z*kb;d=-z*twist;
        value=(n/2)*sum(w.*(a.*b-d.^2)./sqrt(a.^2+d.^2));
    end
end

function [x,w]=gaussRule(n)
% 模块用途：固定Gauss-Legendre节点，积分阶数及差值作为数值分辨率单列。
persistent orders nodes weights
if isempty(orders),orders=[];nodes={};weights={};end
j=find(orders==n,1);if ~isempty(j),x=nodes{j};w=weights{j};return;end
beta=(1:n-1)'./sqrt(4*(1:n-1)'.^2-1);[V,D]=eig(diag(beta,1)+diag(beta,-1));
[x,ids]=sort(diag(D));w=2*V(1,ids)'.^2;orders(end+1)=n;nodes{end+1}=x;weights{end+1}=w;
end
