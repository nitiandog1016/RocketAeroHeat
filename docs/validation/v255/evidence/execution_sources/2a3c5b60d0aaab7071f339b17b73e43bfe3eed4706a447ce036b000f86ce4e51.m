function [result,solution]=solveRegularizedEntropyIntegral(source,options)
%SOLVEREGULARIZEDENTROPYINTEGRAL 驻点正则化的变熵外缘—层流积分候选。
% 模块一：保持v254供体/压力/物性/热相关式，核验身份和原温域。
% 模块二：z=log(theta)、W=I/(K*theta^4)，成熟ode113自适应积分。
% 模块三：原生连续解直接评价每个目标，不用粗站线性插值填原面。
% 模块四：输出原生接受步、代数账与ODE连续缺陷；实际剖面/物理资格未解。
if nargin<2,options=struct;end
o=struct('relativeTolerance',1e-8,'absoluteTolerance',1e-10,'startAngleRad',1e-4,'metricScale',1);
for key=fieldnames(options)'
    assert(isfield(o,key{1}),'RAH:RegularIntegral:Option','未知选项。');o.(key{1})=options.(key{1});
end
validateattributes(o.relativeTolerance,{'numeric'},{'scalar','finite','positive','<=',1e-6});
validateattributes(o.absoluteTolerance,{'numeric'},{'scalar','finite','positive','<=',1e-8});
validateattributes(o.startAngleRad,{'numeric'},{'scalar','finite','positive','<=',1e-4});
validateattributes(o.metricScale,{'numeric'},{'scalar','finite','positive'});
assert(strcmp(source.sourceRole,'engineering-candidate')&&strlength(string(source.sourceSHA256))==64, ...
    'RAH:RegularIntegral:Source','必须保留工程供体身份和准确源。');
mass=source.sourceMassCoordinate(:);entropy=source.sourceEntropyJkgK(:);
assert(mass(1)==0&&numel(mass)==numel(entropy)&&all(diff(mass)>0)&&all(diff(entropy)<=0), ...
    'RAH:RegularIntegral:Mass','质量/熵供体无效。');
targets=source.targetAnglesRad(:)';assert(all(isfinite(targets)&targets>=0&targets<=pi/3), ...
    'RAH:RegularIntegral:Angle','超有限球帽角域，不外推。');
assert(~any(targets>0&targets<o.startAngleRad),'RAH:RegularIntegral:Seed','种子以内非零位置不伪造数值。');
p0=source.stagnationPressurePa;T0=source.stagnationTemperatureK;H=source.totalEnthalpyJkg;
pinf=source.pressurePa;Uinf=source.speedMps;R=source.radiusM;scale=o.metricScale;
up=phaseAwareAirState('PT',pinf,source.temperatureK,'linearTable');
stag=phaseAwareAirState('PT',p0,T0,'linearTable');wall=phaseAwareAirState('PT',p0,300,'linearTable');
assert(abs(up.enthalpyJkg+.5*Uinf^2-H)<=1e-8*abs(H)&&abs(stag.enthalpyJkg-H)<=16*eps*abs(H), ...
    'RAH:RegularIntegral:Energy','供体总焓/驻点不一致，拒绝稳定差量恢复。');
assert(abs(entropy(1)-stag.entropyJkgK)<=1e-7,'RAH:RegularIntegral:Entropy','供体驻点熵不同。');
pp=pchip(mass,entropy-entropy(1));[breaks,coefs]=unmkpp(pp);
dp=mkpp(breaks,coefs(:,1:3).*[3 2 1]);dS0=ppval(dp,0);
cm=5.27*.4696;C0=wall.densityKgm3*wall.dynamicViscosityPas*Uinf;
a0=sqrt(2*(p0-pinf)/(stag.densityKgm3*R^2));coef=-2*T0*dS0*cm*sqrt(wall.densityKgm3*wall.dynamicViscosityPas/2);
a=fzero(@(v)(v/a0)^2-1-coef*sqrt(v)/a0^2,[a0 16*a0],optimset('TolX',1e-10,'Display','off'));
K=scale^2*a*R^4/(4*Uinf);A=1-pinf/p0;
c60=phaseAwareAirCaloric(60);offset=c60.h-60*c60.cp;Hw=H-offset;zeta=(wall.enthalpyJkg-offset)/Hw;
G7=(1.116+.411*.5^.686)/(1+.527*.5^.686);
qref=.768*1.033/G7*wall.Prandtl^(-.6)*sqrt(a0)*(wall.densityKgm3*wall.dynamicViscosityPas)^.1* ...
    (stag.densityKgm3*stag.dynamicViscosityPas)^.4*(H-wall.enthalpyJkg);
root=optimset('TolX',1e-14,'Display','off');sampleAngles=unique([0 targets o.startAngleRad pi/3]);
stepZ=log(o.startAngleRad);stepW=1;solution=[];failure='';
% OutputFcn保存真实接受步；失败时也返回可取得前缀，不补一个成功终点。
try
    solution=ode113(@rhs,[log(o.startAngleRad) log(pi/3)],1, ...
        odeset('RelTol',o.relativeTolerance,'AbsTol',o.absoluteTolerance,'MaxStep',.3,'Refine',1,'OutputFcn',@record));
    assert(solution.x(end)==log(pi/3),'RAH:RegularIntegral:Terminal','原生求解未到60度终点。');
catch err
    failure=[err.identifier '：' err.message];
end
if isempty(solution)
    result=struct('available',false,'failure',failure,'acceptedLogAngle',stepZ,'acceptedRegularIntegral',stepW, ...
        'physicalQualification',false,'profileReadyForCFD',false);return;
end
N=numel(sampleAngles);I=nan(1,N);Psi=I;Te=I;U=I;rho=I;S=I;delta=I;q=I;beta=I;defect=I;momentum=I;energy=I;entropyError=I;
I(1)=0;Psi(1)=0;Te(1)=T0;U(1)=0;rho(1)=stag.densityKgm3;S(1)=entropy(1);
delta(1)=cm*sqrt(wall.densityKgm3*wall.dynamicViscosityPas/(2*a))/stag.densityKgm3;
beta(1)=a0^2/(2*a^2);q(1)=qref*sqrt(a/a0)*G7*(1+.527*beta(1)^.686)/(1.116+.411*beta(1)^.686);
defect(1)=NaN;momentum(1)=0;energy(1)=(stag.enthalpyJkg-H)/H;entropyError(1)=0;
[W,dW]=deval(solution,log(sampleAngles(2:end)));
for j=2:N
    t=sampleAngles(j);st=edge(t,W(j-1));I(j)=K*t^4*W(j-1);Psi(j)=st.Psi;Te(j)=st.T;U(j)=st.U;rho(j)=st.rho;S(j)=st.S;
    metric=scale*R*sin(t);p=st.p;delta(j)=Psi(j)/(st.rho*st.U*metric);
    proxy=2*(p0-pinf)*sin(t)*cos(t)/(R*st.rho*st.U);te=(H-st.dropH-offset)/Hw;
    beta(j)=2*I(j)*proxy*Uinf/(te*(p/p0)*st.U^2*metric^2);
    assert(isfinite(beta(j))&&beta(j)>=0,'RAH:RegularIntegral:Beta','压力梯度代理不可代数求值。');
    shape=G7*(1+.527*beta(j)^.686)/(1.116+.411*beta(j)^.686)*(1.1-.1625*te+.0625*te^2)*(.85+.15*te-zeta)/(1-zeta);
    q(j)=qref*(p/p0)*(st.U/Uinf)*metric*shape/sqrt(4*a0/Uinf*I(j));
    g=(p/p0)*(st.U/Uinf)*metric^2;derivative=K*t^3/R*(4*W(j-1)+dW(j-1));
    defect(j)=(derivative-g)/g;
    momentum(j)=-st.rho*st.T*ppval(dp,Psi(j)/scale)*Psi(j)*g/(2*I(j)*scale);
    check=phaseAwareAirState('PT',p,st.T,'linearTable');
    energy(j)=(check.enthalpyJkg+.5*st.U^2-H)/H;entropyError(j)=(check.entropyJkgK-st.S)/max(1,abs(entropy(1)));
end
[found,indices]=ismember(targets,sampleAngles);assert(all(found));
result=struct('schema','rah-regular-entropy-integral-v1','available',isempty(failure),'failure',failure,'sourceSHA256',source.sourceSHA256, ...
    'method','CR-112180质量匹配/CR-111922积分/CR-111921热相关式；logθ与I/θ⁴正则化＋ode113', ...
    'scope','拟合旋转球帽0—60度；Tw300K；小横流工程候选；实际剖面/三维动量/实验资格未知', ...
    'options',o,'angleRad',sampleAngles,'arcM',R*sampleAngles,'metricM',scale*R*sin(sampleAngles),'integralI',I,'massCoordinate',Psi, ...
    'pressurePa',p0*(1-A*sin(sampleAngles).^2),'edgeTemperatureK',Te,'edgeSpeedMps',U,'edgeDensityKgm3',rho,'edgeEntropyJkgK',S, ...
    'massThicknessM',delta,'momentumThicknessM',delta/5.27,'heatFluxWm2',q,'betaPressureProxy',beta,'alongMomentumDefectPam',momentum, ...
    'denseODERelativeDefect',defect,'maximumDenseODERelativeDefect',max(abs(defect(2:end))), ...
    'maximumEnergyRelativeError',max(abs(energy)),'maximumEntropyRelativeError',max(abs(entropyError)), ...
    'referenceVelocityGradient',a0,'jointLeadingVelocityGradient',a,'sourcePchipEntropySlopeAtZero',dS0, ...
    'sourceSlopeRelativeDifference',dS0/source.entropyDerivativeByMassAtZero-1,'workingEnthalpyOffsetJkg',offset, ...
    'targetIndices',indices,'targetAnglesRad',targets,'acceptedLogAngle',solution.x,'acceptedRegularIntegral',solution.y, ...
    'acceptedSteps',numel(solution.x)-1,'outputSamples',N,'nativeTerminalReached',solution.x(end)==log(pi/3), ...
    'physicalQualification',false,'newQualifiedAreaM2',0,'profileReadyForCFD',false,'actualProfileSolved',false, ...
    'nearZeroTreatment','θ=0解析极限；0<θ<种子不输出数值；固定半种子独立控制检验有限起步', ...
    'correlationRangeExperimentallyValidated',false,'transverseMomentumSolved',false);

    function st=edge(t,w)
        assert(isfinite(w)&&w>0,'RAH:RegularIntegral:Positive','正则积分量非正，不裁态。');
        psi=cm*sqrt(2*C0*K*w)*t^2;massSource=psi/scale;
        assert(massSource>=0&&massSource<=mass(end),'RAH:RegularIntegral:Capacity','有限质量供体超域。');
        dEntropy=ppval(pp,massSource);lp=log1p(-A*sin(t)^2);
        wanted=-dEntropy-up.gasConstantJkgK*lp;
        drop=fzero(@(v)thermalResidual(v,wanted),[0 T0-60],root);
        [dh,~]=phaseAwareAirCaloricDrop(T0,drop);tt=T0-drop;p=p0*exp(lp);
        assert(dh>0,'RAH:RegularIntegral:Kinetic','非零点动能非正。');
        checked=phaseAwareAirState('PT',p,tt,'linearTable'); % 每次真实试探同样检查联合相态域。
        st=struct('T',tt,'U',sqrt(2*dh),'p',p,'rho',checked.densityKgm3,'S',entropy(1)+dEntropy,'dropH',dh,'Psi',psi);
    end
    function r=thermalResidual(drop,wanted)
        [~,ss]=phaseAwareAirCaloricDrop(T0,drop);r=ss-wanted;
    end
    function dw=rhs(z,w)
        t=exp(z);st=edge(t,w);dw=4*((st.p/p0)*(st.U/(a*R*t))*(sin(t)/t)^2-w);
    end
    function stop=record(z,w,flag)
        stop=false;if isempty(flag),stepZ=[stepZ z(:)'];stepW=[stepW w(:)'];end
    end
end
