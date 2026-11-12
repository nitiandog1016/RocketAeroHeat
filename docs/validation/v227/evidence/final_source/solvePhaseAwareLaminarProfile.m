function r=solvePhaseAwareLaminarProfile(edge,station,Tw)
%SOLVEPHASEAWARELAMINARPROFILE 相态空气局部相似层流速度/总焓剖面。
% 模块用途：NASA TP-3271式6–7/17/25/54/56–61的局部相似近似，
% 使用同Cp积分、线性mu-k及变Pr；实际积分密度加权厚度，绝非厚度换算。
% 温度logit只约束数学量热范围，不裁状态。最终p/T联合域另严格审核。
% 这不是全三维边界层，也不提供转捩、分离、实验或CFD入口资格。
%% 1. 同模型状态、能量和沿程动量相容性
r=struct('available',false,'reason','输入未接受','physicalQualification',false, ...
    'heatFlux',NaN,'wallShearStress',NaN,'delta',NaN,'deltaStar',NaN, ...
    'thetaPhysical',NaN,'massThickness',NaN,'profile',struct,'edge',edge,'station',station, ...
    'wallTemperature',Tw,'odeResidual',NaN,'boundaryResidual',NaN);
try
    a=phaseAwareAirState('PT',edge.pressurePa,edge.temperatureK,'linearTable');
    phaseAwareAirState('PT',edge.pressurePa,Tw,'linearTable');
    U=edge.speedMps;Te=edge.temperatureK;pe=edge.pressurePa;beta=station.beta;
    assert(isscalar(U)&&isfinite(U)&&U>=0&&isscalar(Tw)&&Tw>60&&Tw<2000);
    assert(isscalar(beta)&&isfinite(beta)&&beta>=0&&beta<=1);
    c60=phaseAwareAirCaloric(60);offset=c60.h-60*c60.cp;
    H=edge.totalEnthalpyJkg-offset;
    assert(H>0&&abs(edge.totalEnthalpyJkg-a.enthalpyJkg-U^2/2)/H<1e-6);
    if U==0 && station.metricM==0 && station.xi==0
        grad=station.velocityGradientPerSecond;
        assert(grad>0&&abs(beta-.5)<1e-10&&station.pressureGradientPam==0);
        assert(abs(station.pressureCurvaturePam2+a.densityKgm3*grad^2)/abs(station.pressureCurvaturePam2)<1e-6);
        K=sqrt(a.dynamicViscosityPas/a.densityKgm3/(2*grad));
    else
        assert(U>0&&station.metricM>0&&station.xi>0);
        K=sqrt(2*station.xi)/(a.densityKgm3*U*station.metricM);
        b=2*station.xi*station.velocityGradientPerSecond/(a.densityKgm3*a.dynamicViscosityPas*U^2*station.metricM^2);
        assert(abs(beta-b)<1e-6*max(1,abs(b)));
        dp=station.pressureGradientPam;du=a.densityKgm3*U*station.velocityGradientPerSecond;
        assert(abs(dp+du)<=1e-6*max(abs([dp du])));
    end
catch err
    r.reason=['输入拒绝: ' err.message];return
end
r.thicknessScale=K;r.enthalpyOffset=offset;r.shiftedTotalEnthalpy=H;
%% 2. 固定五状态BVP；n包含总焓扩散和黏性耗散，Pr随温度变化
% y=[f,F,m,z,n]；m=lF'；g=(h+U²F²/2-offset)/H。
% n=(l/Pr)g'+(1-1/Pr)(U²/H)Fm。不能遗漏变Pr通量或黏性耗散。
try
    initial=bvpinit(linspace(0,12,81),@guess);
    opt=bvpset('RelTol',1e-7,'AbsTol',1e-9,'NMax',12000,'Stats','off','Vectorized','on');
    sol=bvp4c(@ode,@bc,initial,opt);
    eta=linspace(0,12,1201);[Y,Yp]=deval(sol,eta);[T,l,D,C,Pr,g,dz]=props(Y);
    final=phaseAwareAirState('PT',pe,T,'linearTable'); %#ok<NASGU>
    rhs=ode(eta,Y);r.odeResidual=max(abs(Yp-rhs)./(1+abs(rhs)),[],'all');
    r.boundaryResidual=max(abs(bc(Y(:,1),Y(:,end))));
    F=Y(2,:);m=Y(3,:);n=Y(5,:);Fp=m./l;
    gp=Pr./l.*(n-(1-1./Pr).*(U^2/H).*F.*m);
    r.minimumVelocityRatio=min(F);r.maximumVelocityRatio=max(F);
    assert(r.odeResidual<=2e-5&&r.boundaryResidual<=1e-7&&m(1)>0&&min(F)>=-1e-6, ...
        'RAH:PhaseProfile:Residual','剖面残差/边界/无逆流附着支检查失败。');
    fs=max(abs(Fp));gs=max(.005*max(abs(gp)),1e-8);
    defect=defects(Y);bad=find(defect>1,1,'last');
    assert(~isempty(bad)&&bad<numel(eta)&&eta(bad+1)<=11, ...
        'RAH:PhaseProfile:Edge','有限边缘尚未解析或外界缓冲不足。');
    ee=fzero(@(x)max(defects(deval(sol,x)))-1,[eta(bad) eta(bad+1)]);
    nodes=linspace(0,ee,601);states=deval(sol,nodes);[t,ell,ratio,~,~,gg]=props(states);
    physical=phaseAwareAirState('PT',pe,t,'linearTable');
    ints=integral(@integrands,0,ee,'ArrayValued',true,'RelTol',1e-8,'AbsTol',1e-10);
    r.delta=K*ints(1);r.massThickness=K*ints(2);r.thetaPhysical=K*ints(3);
    % MATLAB嵌套函数共享局部m；外缘搜索会改写它。壁值必须取已保存的壁列。
    r.deltaStar=r.delta-r.massThickness;r.heatFlux=a.dynamicViscosityPas*H/K*Y(5,1);
    r.wallShearStress=a.dynamicViscosityPas*U/K*Y(3,1);
    assert(all(isfinite([r.delta r.massThickness r.thetaPhysical r.heatFlux r.wallShearStress])) ...
        &&r.delta>0&&r.massThickness>0&&r.thetaPhysical>0,'RAH:PhaseProfile:Integral','物理积分未通过。');
    r.etaEdge=ee;r.profile=struct('eta',nodes,'states',states,'temperature',t, ...
        'speed',U*states(2,:),'density',physical.densityKgm3,'totalEnthalpy',offset+H*gg, ...
        'viscosity',physical.dynamicViscosityPas,'conductivity',physical.conductivityWmK, ...
        'physicalY',K*cumtrapz(nodes,ratio),'velocityDerivative',states(3,:)./ell);
    r.numericalProfile=struct('eta',sol.x,'states',sol.y);
    r.massDefinitionResidual=abs(r.delta-r.deltaStar-r.massThickness)/r.delta;
    r.heatFluxNumericallyResolved=abs(n(1))>10*(1e-9+1e-7*abs(n(1)));
    r.available=true;r.reason='已计算：同相态变Pr局部相似工程剖面；未获物理资格';
catch err
    r.reason=['剖面拒绝: ' err.identifier ' ' err.message];
end
    function dy=ode(~,y)
        [~,l,D,C,Pr,~,dz]=props(y);F=y(2,:);m=y(3,:);n=y(5,:);
        Fp=m./l;gp=Pr./l.*(n-(1-1./Pr).*(U^2/H).*F.*m);
        Tp=Pr./(l.*C).*(H*n-U^2*F.*m);
        dy=[F;Fp;-y(1,:).*Fp+beta*(F.^2-D);Tp./dz;-y(1,:).*gp];
    end
    function z=bc(ya,yb)
        z=[ya(1);ya(2);ya(4)-log((Tw-60)/(2000-Tw));yb(2)-1;yb(4)-log((Te-60)/(2000-Te))];
    end
    function y=guess(x)
        F=tanh(x);Fp=1-F^2;t=Tw+(Te-Tw)*F;
        y=[log(cosh(x));F;0;log((t-60)/(2000-t));0];
        [~,l,~,C,Pr]=props(y);m=l*Fp;
        n=l/Pr*C/H*(Te-Tw)*Fp+(U^2/H)*F*m;
        y(3)=m;y(5)=n;
    end
    function [T,l,D,C,Pr,g,dz]=props(y)
        T=60+1940./(1+exp(-y(4,:)));dz=(T-60).*(2000-T)/1940;
        assert(all(isfinite(dz)&dz>0),'RAH:PhaseProfile:Bound','温度变换到达表端；不裁状态或扩表。');
        c=phaseAwareAirCaloric(T);[mu,k]=transport(T);C=c.cp;Pr=C.*mu./k;
        D=T/Te;l=mu/a.dynamicViscosityPas./D;g=(c.h+.5*U^2*y(2,:).^2-offset)/H;
    end
    function d=defects(y)
        [t,l,~,~,Pr,g]=props(y);F=y(2,:);m=y(3,:);
        gp=Pr./l.*(y(5,:)-(1-1./Pr).*(U^2/H).*F.*m);
        d=max([abs(1-F)/.005;abs(1-g)/.005;abs(t/Te-1)/.005; ...
            abs(Te./t-1)/.005;abs(m./l)/fs/.005;abs(gp)/gs],[],1);
    end
    function v=integrands(x)
        y=deval(sol,x);[~,~,D]=props(y);v=[D;y(2);y(2)*(1-y(2))];
    end
end
function [mu,k]=transport(T)
% 模块用途：缓存同一Cp节点的已具名Lemmon线性输运，不延拓、不固定Pr。
persistent x mn kn
if isempty(x)
    d=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'data/thermo/phase_aware_dry_air_v1.json')));
    x=d.temperatureK(:);z=log(x/103.3);o=exp(.431-.4623*z+.08406*z.^2+.005341*z.^3-.00331*z.^4);
    mn=2.66958e-8*sqrt(28.9586*x)./(.36^2*o);tau=132.6312./x;
    kn=.001308*(mn*1e6)+.001405*tau.^(-1.1)-.001036*tau.^(-.3);
end
mu=interp1(x,mn,T,'linear');k=interp1(x,kn,T,'linear');
end
