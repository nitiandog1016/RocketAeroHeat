function result=solvePhaseAwareEntropyIntegral(source,options)
%SOLVEPHASEAWAREENTROPYINTEGRAL 变熵外缘与原配套层流积分质量的联合工程计算。
% 模块一：核对同相态/来源/冷等温壁与有限球帽范围；输入不升格物理资格。
% 模块二：隐式梯形推进I，同时匹配Ψ=Cm*sqrt(2*C0*I)和外缘熵/速度。
% 模块三：按CR-111921热相关式返回q/厚度/全部入口，实际剖面仍未知。
% 来源：CR-112180式65—71；CR-111922印刷53/57页的CTHM/Y4/THM；
% CR-111921式146/149—159。这里的theta不是Thwaites参考厚度。
if nargin<2,options=struct;end
o=struct('wallTemperatureK',300,'stepDeg',.5,'metricScale',1,'mode','coupled');
names=fieldnames(options);for k=1:numel(names),assert(isfield(o,names{k}),'RAH:EntropyIntegral:Option','未知选项。');o.(names{k})=options.(names{k});end
validateattributes(o.wallTemperatureK,{'numeric'},{'scalar','finite','positive'});
assert(o.wallTemperatureK==300,'RAH:EntropyIntegral:Wall','本次完整源合同仅采用300 K冷等温壁，非等温/高壁温须另证。');
validateattributes(o.stepDeg,{'numeric'},{'scalar','finite','positive','<=',1});
validateattributes(o.metricScale,{'numeric'},{'scalar','finite','positive'});
assert(ismember(o.mode,{'coupled','normal-shock-reference'}),'RAH:EntropyIntegral:Mode','未知分支。');
required={'pressurePa','temperatureK','speedMps','totalEnthalpyJkg','stagnationPressurePa', ...
    'stagnationTemperatureK','stagnationEntropyJkgK','radiusM','sourceMassCoordinate', ...
    'sourceEntropyJkgK','entropyDerivativeByMassAtZero','targetAnglesRad','sourceRole','sourceSHA256'};
assert(all(isfield(source,required))&&strcmp(source.sourceRole,'engineering-candidate')&&strlength(string(source.sourceSHA256))==64, ...
    'RAH:EntropyIntegral:Source','缺少工程供体身份；不能传入假accepted-edge。');
mass=source.sourceMassCoordinate(:);entropy=source.sourceEntropyJkgK(:);
assert(numel(mass)==numel(entropy)&&mass(1)==0&&all(diff(mass)>0)&&all(diff(entropy)<=0)&&all(isfinite(entropy)), ...
    'RAH:EntropyIntegral:Source','供体质量/熵不完整或非单调。');
entropyPP=pchip(mass,entropy); % 同一供体插值只建一次，根迭代复用。
targets=source.targetAnglesRad(:)';assert(~isempty(targets)&&all(isfinite(targets))&&all(targets>=0&targets<=pi/3), ...
    'RAH:EntropyIntegral:Angle','入口角超出0—60度供体范围，禁止外推。');
p0=source.stagnationPressurePa;T0=source.stagnationTemperatureK;H=source.totalEnthalpyJkg;
R=source.radiusM;Uinf=source.speedMps;pinf=source.pressurePa;scale=o.metricScale;
assert(R>0&&p0>pinf&&Uinf>0&&source.entropyDerivativeByMassAtZero<0);
up=phaseAwareAirState('PT',pinf,source.temperatureK,'linearTable');
stag=phaseAwareAirState('PT',p0,T0,'linearTable');wall=phaseAwareAirState('PT',p0,o.wallTemperatureK,'linearTable');
assert(abs(up.enthalpyJkg+.5*Uinf^2-H)<1e-8*max(1,abs(H))&&abs(stag.enthalpyJkg-H)<1e-8*max(1,abs(H)), ...
    'RAH:EntropyIntegral:Energy','供体不是同一总焓。');
assert(abs(entropy(1)-stag.entropyJkgK)<1e-7,'RAH:EntropyIntegral:Source','零质量端不是同一正激波驻点熵。');
cm=5.27*.4696;C0=wall.densityKgm3*wall.dynamicViscosityPas*Uinf;
a0=sqrt(2*(p0-pinf)/(stag.densityKgm3*R^2));a=a0;
if strcmp(o.mode,'coupled')
    coefficient=-2*T0*source.entropyDerivativeByMassAtZero*cm*sqrt(wall.densityKgm3*wall.dynamicViscosityPas/2);
    fun=@(x)(x/a0)^2-1-coefficient*sqrt(x)/a0^2;
    assert(fun(16*a0)>0,'RAH:EntropyIntegral:Start','正支未在登记括区内，拒绝起步。');
    a=fzero(fun,[a0 16*a0],optimset('TolX',1e-9,'Display','off'));
end
c60=phaseAwareAirCaloric(60);offset=c60.h-60*c60.cp;Hwork=H-offset;zeta=(wall.enthalpyJkg-offset)/Hwork;
assert(Hwork>0&&zeta>0&&zeta<1,'RAH:EntropyIntegral:Caloric','冷壁工作焓比无效。');
G7=(1.116+.411*.5^.686)/(1+.527*.5^.686);
qref=.768*1.033/G7*wall.Prandtl^(-.6)*sqrt(a0)* ...
    (wall.densityKgm3*wall.dynamicViscosityPas)^.1*(stag.densityKgm3*stag.dynamicViscosityPas)^.4*(H-wall.enthalpyJkg);
start=.001;maxAngle=max(targets);assert(maxAngle>start);
angles=unique([0 start deg2rad(o.stepDeg):deg2rad(o.stepDeg):maxAngle targets maxAngle]);
assert(~any(angles>0&angles<start),'RAH:EntropyIntegral:Start','有限正则起步前的非零点未求解。');
N=numel(angles);I=nan(1,N);Psi=I;Te=I;U=I;rho=I;S=I;thetaM=I;deltaM=I;q=I;beta=I;proxy=I;implicit=I;dSds=I;
pressure=pinf+(p0-pinf)*cos(angles).^2;metric=scale*R*sin(angles);s=R*angles;
I(1)=0;Psi(1)=0;Te(1)=T0;U(1)=0;rho(1)=stag.densityKgm3;S(1)=entropy(1);
thetaM(1)=.4696*sqrt(wall.densityKgm3*wall.dynamicViscosityPas/(2*a))/stag.densityKgm3;deltaM(1)=5.27*thetaM(1);
beta(1)=a0^2/(2*a^2);q(1)=qref*sqrt(a/a0)*G7*(1+.527*beta(1)^.686)/(1.116+.411*beta(1)^.686);
proxy(1)=a0^2/a;implicit(1)=0;
I(2)=scale^2*a*s(2)^4/(4*Uinf);implicit(2)=NaN;
[state,gOld]=evaluate(I(2),2);assign(2,state);
% 根括区上端略内置以避免log/exp往返越过表端一ulp；不裁求得的物理态。
capacity=(scale*mass(end)/cm)^2/(2*C0)*(1-16*eps);rootOptions=optimset('TolX',1e-11,'Display','off');
failure='';lastSolved=2;
for k=3:N
    try
    ds=s(k)-s(k-1);previous=I(k-1);
    f=@(logI)residual(exp(logI),k,previous,gOld,ds);
    assert(f(log(previous))<0&&f(log(capacity))>0,'RAH:EntropyIntegral:Capacity','隐式根越出有限供体范围，停止下游。');
    I(k)=exp(fzero(f,[log(previous) log(capacity)],rootOptions));
    implicit(k)=abs(f(log(I(k))));assert(implicit(k)<=1e-8,'RAH:EntropyIntegral:Residual','隐式沿程积分未闭合。');
    [state,gOld]=evaluate(I(k),k);assign(k,state);
    lastSolved=k;
    catch err
        failure=[err.identifier '：' err.message];break
    end
end
known=1:lastSolved;
assert(all(isfinite([I(known) Psi(known) Te(known) U(known) rho(known) S(known) thetaM(known) deltaM(known) q(known) beta(known)])), ...
    'RAH:EntropyIntegral:Finite','已接受前缀存在非有限计算量。');
% β中的压力梯度代理与真正的变熵沿程导数必须区分，不拿前者冒充Euler闭合。
pp=pchip(s(known),S(known));[breaks,coefs,~,~]=unmkpp(pp);power=size(coefs,2)-1;
derivative=mkpp(breaks,coefs(:,1:end-1).*(power:-1:1));dSds(known)=ppval(derivative,s(known));
pressureSlope=-(p0-pinf)*sin(2*angles)/R;
trueMomentumDefect=-rho.*Te.*dSds;
[found,index]=ismember(targets,angles);assert(all(found));
result=struct('schema','rah-phase-entropy-integral-v1','available',isempty(failure),'failure',failure,'acceptedStations',lastSolved, ...
    'sourceRole',source.sourceRole,'sourceSHA256',source.sourceSHA256, ...
    'method','CR-112180质量匹配＋CR-111922壁面性质积分厚度＋CR-111921层流热相关式；隐式梯形', ...
    'scope','拟合旋转球帽0—60度、冷等温壁300K、小横流工程近似；无实际三维剖面/湍流/实验资格', ...
    'options',o,'angleRad',angles,'arcM',s,'metricM',metric,'integralI',I,'massCoordinate',Psi, ...
    'pressurePa',pressure,'edgeTemperatureK',Te,'edgeSpeedMps',U,'edgeDensityKgm3',rho,'edgeEntropyJkgK',S, ...
    'momentumThicknessM',thetaM,'massThicknessM',deltaM,'heatFluxWm2',q,'betaPressureProxy',beta, ...
    'velocityGradientPressureProxy',proxy,'entropyGradientJkgKm',dSds,'alongMomentumDefectPam',trueMomentumDefect, ...
    'pressureGradientPam',pressureSlope,'implicitRelativeResidual',implicit,'targetIndices',index,'targetAnglesRad',targets, ...
    'startAngleRad',start,'referenceVelocityGradient',a0,'jointLeadingVelocityGradient',a,'referenceStagnationHeatFlux',qref, ...
    'workingEnthalpyOffsetJkg',offset,'workingTotalEnthalpyJkg',Hwork,'wallEnthalpyRatio',zeta,'correlationPrandtl',wall.Prandtl, ...
    'maximumImplicitResidual',max(implicit(3:end)),'physicalQualification',false,'newQualifiedAreaM2',0, ...
    'correlationRangeExperimentallyValidated',false,'betaAboveOneStations',sum(beta(known)>1), ...
    'betaAlgebraicContract','有限非负；CR-111922印刷57页无上截断；不表示实验有效范围已证', ...
    'profileReadyForCFD',false,'actualProfileSolved',false,'transverseMomentumSolved',false,'fixedHeatAreaQualification',false);

    function [aState,g]=evaluate(integralValue,j)
        psi=cm*sqrt(2*C0*integralValue);sourcePsi=psi/scale;
        assert(sourcePsi>=0&&sourcePsi<=mass(end),'RAH:EntropyIntegral:SourceRange','质量匹配越出上游表，不外推。');
        if strcmp(o.mode,'normal-shock-reference'),ss=entropy(1);else,ss=ppval(entropyPP,sourcePsi);end
        wanted=ss+up.gasConstantJkgK*log(pressure(j)/100000);
        tt=fzero(@(x)entropyResidual(x,wanted),[60 T0],optimset('TolX',1e-10,'Display','off'));
        at=phaseAwareAirState('PT',pressure(j),tt,'linearTable');u2=2*(H-at.enthalpyJkg);
        assert(u2>0,'RAH:EntropyIntegral:Kinetic','外缘动能非正，禁止裁速度。');
        vv=sqrt(u2);g=pressure(j)/p0*(vv/Uinf)*metric(j)^2;
        aState=struct('T',tt,'U',vv,'rho',at.densityKgm3,'S',ss,'h',at.enthalpyJkg,'Psi',psi);
    end
    function value=entropyResidual(T,wanted)
        at=phaseAwareAirCaloric(T);value=at.s0-wanted;
    end
    function value=residual(current,j,last,gLast,ds)
        [~,gNow]=evaluate(current,j);balance=last+.5*ds*(gLast+gNow);value=(current-balance)/balance;
    end
    function assign(j,st)
        Psi(j)=st.Psi;Te(j)=st.T;U(j)=st.U;rho(j)=st.rho;S(j)=st.S;
        thetaM(j)=.4696*sqrt(2*C0*I(j))/(st.rho*st.U*metric(j));deltaM(j)=5.27*thetaM(j);
        dp=-(p0-pinf)*sin(2*angles(j))/R;proxy(j)=-dp/(st.rho*st.U);te=(st.h-offset)/Hwork;
        beta(j)=2*I(j)*(proxy(j)/Uinf)/(te*(pressure(j)/p0)*(st.U/Uinf)^2*metric(j)^2);
        % 原程序BBAR直接进入0.686次幂，无beta<=1限制。这里只检查代数可计算性；
        % 不复用旧BVP参数范围，也不把未证明的实验适用范围升级为合格。
        assert(isfinite(beta(j))&&beta(j)>=0,'RAH:EntropyIntegral:PressureRange','相关式压力梯度参数非有限或负，拒绝求值。');
        shape=G7*(1+.527*beta(j)^.686)/(1.116+.411*beta(j)^.686)*(1.1-.1625*te+.0625*te^2)*(.85+.15*te-zeta)/(1-zeta);
        q(j)=qref*(pressure(j)/p0)*(st.U/Uinf)*metric(j)*shape/sqrt(4*a0/Uinf*I(j));
    end
end
