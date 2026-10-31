function jump=shockLayerNormalJumpBatch(upstream,gas,normalSpeed)
%SHOCKLAYERNORMALJUMPBATCH 同一上游和物性下批量正常压缩跃迁及对入射法向速度的导数。
% 模块一：定比热闭式或同cp能量压缩根，冻结已收敛根，不取平凡密度比1。
% 模块二：由隐式Rankine-Hugoniot导数得到dS/dw1，轴线上也不除以半径。
% 模块三：核对熵增、压缩、法向亚声速及能量；不将跃迁结果赋给真实壁面。
assert(isreal(normalSpeed)&&all(isfinite(normalSpeed)&normalSpeed>0,'all'), ...
    'ShockPatch:NormalMach','入射法向速度须正且有限。');
first=surfaceFlowThermodynamicState(upstream.temperature,upstream.pressure,gas);
assert(abs(upstream.pressure-upstream.density*gas.R*upstream.temperature)<1e-10*upstream.pressure, ...
    'ShockPatch:Upstream','上游状态方程不一致。');
w=normalSpeed;R=gas.R;T1=upstream.temperature;p1=upstream.pressure;rho1=upstream.density;
M=w/sqrt(first.gammaLocal*R*T1);assert(all(M>1+1e-5,'all'), ...
    'ShockPatch:NormalMach','共同激波来源的法向Mach须大于1+1e-5。');
H=first.enthalpy+w.^2/2;g=first.gammaLocal;
r=(g+1)*M.^2./((g-1)*M.^2+2);iterations=0;
if strcmp(gas.thermalPropertyMode,'thermallyPerfect')
    lo=1+min(1e-4,.25*(M-1));hi=repmat(1e4,size(w));active=true(size(w));
    for iterations=1:50
        [residual,derivative]=energy(r);active=active&abs(residual)>16*eps*max(1,abs(H));
        if ~any(active,'all'),break;end
        positive=active&residual>0;negative=active&~positive;lo(positive)=r(positive);hi(negative)=r(negative);
        candidate=r-residual./derivative;
        bisect=active&(~isfinite(candidate)|derivative>=0|candidate<=lo|candidate>=hi);
        candidate(bisect)=(lo(bisect)+hi(bisect))/2;r(active)=candidate(active);
    end
    assert(~any(active,'all'),'ShockPatch:CompressionRoot','同cp非平凡压缩根未收敛。');
end
p=p1+rho1*w.^2.*(1-1./r);rho=rho1*r;T=p./(rho*R);v=w./r;
[cp,h,s0]=shockLayerThermalBatch(T,gas);S=s0-R*log(p/gas.referencePressurePa);
gamma=cp./(cp-R);error=(h+v.^2/2-H)./max(1,abs(H));
assert(all(r>1&p>p1&T>T1&v.^2<gamma.*R.*T&S>first.entropy,'all')&&max(abs(error),[],'all')<2e-11, ...
    'ShockPatch:CompressionRoot','压缩、熵增、能量或激波后法向亚声速未通过。');
dTdr=-T1./r.^2+(w.^2/R).*(2-r)./r.^3;dTdw=2*w/R.*(r-1)./r.^2;
Fr=cp.*dTdr-w.^2./r.^3;Fw=cp.*dTdw+w.*(1./r.^2-1);
assert(all(isfinite(Fr)&abs(Fr)>1e-10*max(1,w.^2),'all'), ...
    'ShockPatch:JumpDerivative','压缩跃迁导数退化。');
dr=-Fw./Fr;dp=rho1*(2*w.*(1-1./r)+w.^2./r.^2.*dr);
dT=T.*(dp./p-dr./r);dS=cp./T.*dT-R./p.*dp;
jump=struct('pressure',p,'temperature',T,'density',rho,'normalSpeed',v,'entropy',S, ...
    'densityRatio',r,'upstreamNormalMach',M,'energyResidual',error, ...
    'entropyDerivativeByNormalSpeed',dS,'densityRatioDerivativeByNormalSpeed',dr, ...
    'pressureDerivativeByNormalSpeed',dp,'temperatureDerivativeByNormalSpeed',dT,'rootIterations',iterations);

    function [E,dE]=energy(ratio)
        temp=T1./ratio+(w.^2/R).*(ratio-1)./ratio.^2;
        below=temp<gas.temperatureRangeK(1);above=temp>gas.temperatureRangeK(2);inside=~(below|above);
        E=nan(size(temp));dE=E;E(below)=-max(1,H(below));E(above)=max(1,H(above));
        [heat,enthalpy]=shockLayerThermalBatch(temp(inside),gas);
        E(inside)=enthalpy+.5*(w(inside)./ratio(inside)).^2-H(inside);
        dE(inside)=heat.*(-T1./ratio(inside).^2+(w(inside).^2/R).*(2-ratio(inside))./ratio(inside).^3) ...
            -w(inside).^2./ratio(inside).^3;
    end
end
