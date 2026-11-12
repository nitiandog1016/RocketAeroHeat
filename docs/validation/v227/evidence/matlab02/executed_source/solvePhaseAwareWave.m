function w=solvePhaseAwareWave(kind,p1,T1,U1,turnDeg)
%SOLVEPHASEAWAREWAVE 同相态空气的法向/弱斜激波和二维等熵膨胀。
% 模块用途：内部外流状态转移核；不解锥流、三维相互作用、脱体驻点或黏性层。
% 保持质量/法向动量/总焓；斜激波保持切向速度。最终p/T联合域严格审核。
% 复用原工程波系的守恒与Prandtl-Meyer结构，显式采用新Cp积分，旧默认不变。
if nargin<5,turnDeg=0;end
kind=validatestring(kind,{'normal','oblique','expansion'});
for value={p1,T1,U1},validateattributes(value{1},{'numeric'},{'real','finite','scalar','positive'});end
validateattributes(turnDeg,{'numeric'},{'real','finite','scalar','>=',0,'<',90});
a=phaseAwareAirState('PT',p1,T1,'linearTable');c1=phaseAwareAirCaloric(T1);R=c1.R;
w=struct('available',false,'kind',kind,'status','未求解','pressurePa',NaN,'temperatureK',NaN, ...
    'speedMps',NaN,'densityKgm3',NaN,'Mach',NaN,'turnDeg',turnDeg,'shockAngleDeg',NaN, ...
    'maximumTurnDeg',NaN,'massResidual',NaN,'momentumResidual',NaN,'energyResidual',NaN, ...
    'entropyChangeJkgK',NaN,'angleResidualRad',NaN,'normalSpeedBefore',NaN,'normalSpeedAfter',NaN, ...
    'tangentialSpeed',NaN,'model',a.model,'transportModel',a.transportModel,'energyReference',a.energyReference, ...
    'sourceSHA256',a.sourceSHA256,'physicalQualification',false, ...
    'method',['相态标准Air：' kind '；Rankine-Hugoniot或二维等熵简单波'], ...
    'scope','60—2000K与压力相态/低密度联合域；局部二维波系，不授三维外缘资格');
if U1<=a.soundSpeedMps,w.status='上游非超声速';return;end
if ~strcmp(kind,'normal')&&turnDeg==0
    finish(p1,T1,U1,0,0,0,0,0);w.status='零转角：严格保持上游';return;
end

%% 1. 法向或弱支斜激波：先量热极线，最终状态另作联合域判定
if ismember(kind,{'normal','oblique'})
    beta=pi/2;
    if strcmp(kind,'oblique')
        bmin=asin(c1.a/U1);
        persistent key peak maximum
        thisKey=[p1 T1 U1];
        if isempty(key)||~isequal(key,thisKey)
            % 不把无量热解的极线点当作零转角。完整极线缺失即具体拒绝。
            [~,~,ok]=normalCore(U1);
            if ~ok,w.status='完整激波极线超出量热域；未判定脱体';return;end
            [pk,negative]=fminbnd(@(b)-polar(b),bmin,pi/2,optimset('TolX',1e-10,'Display','off'));
            peak=pk;maximum=-negative;key=thisKey;
        end
        w.maximumTurnDeg=rad2deg(maximum);
        target=deg2rad(turnDeg);
        if target>maximum,w.status='指定转角超过附体弱支最大角，需脱体模型';return;end
        if target==maximum,beta=peak;
        else,beta=fzero(@(b)polar(b)-target,[bmin peak],optimset('TolX',1e-12,'Display','off'));end
    end
    un1=U1*sin(beta);ut=U1*cos(beta);[T2,un2,ok]=normalCore(un1);
    if ~ok,w.status='法向守恒压缩根不在量热表内';return;end
    p2=p1+a.densityKgm3*un1*(un1-un2);U2=hypot(un2,ut);
    rho2=p2/(R*T2);m=(rho2*un2-a.densityKgm3*un1)/(a.densityKgm3*un1);
    mom=(p2+rho2*un2^2-p1-a.densityKgm3*un1^2)/(p1+a.densityKgm3*un1^2);
    angle=0;if strcmp(kind,'oblique'),angle=beta-atan2(un2,ut)-deg2rad(turnDeg);end
    w.shockAngleDeg=rad2deg(beta);w.normalSpeedBefore=un1;w.normalSpeedAfter=un2;w.tangentialSpeed=ut;
    finish(p2,T2,U2,m,mom,angle,NaN,NaN);
    return;
end

%% 2. 简单膨胀波：dtheta=sqrt(M²-1)*dU/U，在同一总焓下积分
h0=c1.h+.5*U1^2;target=deg2rad(turnDeg);
maxAngle=angleAt(60);w.maximumTurnDeg=rad2deg(maxAngle);
if target>maxAngle,w.status='膨胀要求温度低于60K';return;end
T2=fzero(@(t)angleAt(t)-target,[60 T1],optimset('TolX',1e-10,'Display','off'));
c2=phaseAwareAirCaloric(T2);p2=p1*exp((c2.s0-c1.s0)/R);U2=sqrt(2*(h0-c2.h));
finish(p2,T2,U2,NaN,NaN,angleAt(T2)-target,NaN,NaN);

    function value=angleAt(t)
        value=integral(@integrand,t,T1,'AbsTol',1e-10,'RelTol',1e-9);
    end
    function f=integrand(t)
        c=phaseAwareAirCaloric(t);u2=2*(h0-c.h);M2=u2./c.a.^2;
        assert(all(M2(:)>1),'RAH:PhaseWave:Expansion','膨胀积分经过非超声速状态。');
        f=sqrt(M2-1).*c.cp./u2;
    end
    function value=polar(b)
        if b<=asin(c1.a/U1)+2*eps,value=0;return;end
        [~,u,valid]=normalCore(U1*sin(b));
        assert(valid,'RAH:PhaseWave:Polar','极线量热根缺失，禁止以零角替代。');
        value=b-atan2(u,U1*cos(b));
    end
    function [t,u,ok]=normalCore(v)
        % z=(Un1-Un2)/Un1。除去z=0平凡根，避免近Mach角寻根回到无激波。
        % T(z)为二次式。按60/2000K交点切成合法区间，不在超温区伪造残差符号。
        A=v^2/R;B=A-T1;t=NaN;u=NaN;ok=false;
        if v<=c1.a,return;end
        ze=(B+sqrt(B^2+4*A*(T1-60)))/(2*A);
        cuts=[0 ze];disc=B^2-4*A*(2000-T1);
        if disc>0,cuts=[cuts (B-sqrt(disc))/(2*A) (B+sqrt(disc))/(2*A)];end
        cuts=sort(unique(cuts(cuts>=0&cuts<=ze)));
        for kk=1:numel(cuts)-1
            lo=cuts(kk);hi=cuts(kk+1);tm=T1+B*(lo+hi)/2-A*((lo+hi)/2)^2;
            if tm>2000,continue;end
            fl=residual(lo);fh=residual(hi);
            if fl>0&&fh<=0
                for iteration=1:60
                    mid=(lo+hi)/2;fm=residual(mid);
                    if fm>0,lo=mid;else,hi=mid;end
                end
                z=(lo+hi)/2;t=T1+B*z-A*z^2;u=v*(1-z);
                ok=t>=60&&t<=2000&&z>0&&z<1;return;
            end
        end
        function f=residual(z)
            temp=T1+B*z-A*z^2;
            % 端点公式的数舍入仅允许数十ulp回到准确60/2000交点，不裁物理解。
            endpointTolerance=64*eps(max([T1 A 2000]));
            if (z==cuts(kk)||z==cuts(kk+1))&&abs(temp-60)<endpointTolerance,temp=60;end
            if (z==cuts(kk)||z==cuts(kk+1))&&abs(temp-2000)<endpointTolerance,temp=2000;end
            c=phaseAwareAirCaloric(temp);
            if abs(temp-T1)<1e-3,cm=phaseAwareAirCaloric((temp+T1)/2);meanCp=cm.cp;
            else,meanCp=(c.h-c1.h)/(temp-T1);end
            f=meanCp*(B-A*z)+.5*v^2*(z-2);
        end
    end
    function finish(p,t,u,m,mom,angle,unused1,unused2) %#ok<INUSD>
        try,b=phaseAwareAirState('PT',p,t,'linearTable');
        catch err
            if ismember(err.identifier,{'RAH:PhaseAir:Domain','RAH:PhaseAir:Temperature'}),w.status='波后p/T联合物性域拒绝';return;end
            rethrow(err);
        end
        energy=(b.enthalpyJkg+.5*u^2-c1.h-.5*U1^2)/max(1,abs(c1.h)+.5*U1^2);
        entropy=b.entropyJkgK-a.entropyJkgK;
        residuals=[energy angle m mom];residuals=residuals(isfinite(residuals));
        assert(max(abs(residuals))<=1e-8,'RAH:PhaseWave:Conservation','波系守恒/角残差未通过。');
        if strcmp(kind,'expansion'),assert(abs(entropy)<=1e-7,'RAH:PhaseWave:Entropy','膨胀熵差未过。');
        else,assert(entropy>=-1e-8,'RAH:PhaseWave:Entropy','压缩波熵减。');end
        w.available=true;w.status='同模型波后状态已求解；区域资格需另外验证';
        w.pressurePa=p;w.temperatureK=t;w.speedMps=u;w.densityKgm3=b.densityKgm3;w.Mach=u/b.soundSpeedMps;
        w.massResidual=m;w.momentumResidual=mom;w.energyResidual=energy;w.entropyChangeJkgK=entropy;w.angleResidualRad=angle;
    end
end
