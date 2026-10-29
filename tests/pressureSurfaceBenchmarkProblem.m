function [problem,spec]=pressureSurfaceBenchmarkProblem(kind,overrides)
%PRESSURESURFACEBENCHMARKPROBLEM 为流线内核构造有解析依据的独立表面验证场。
% 模块一：平面均匀流/可压缩自由旋涡检验真实横向压力转向。
% 模块二：圆柱、锥和球面提供精确光滑几何，检查联络及标签横向尺度。
% 模块三：所有p/T均明确为验证输入，不充当真实箭体已求外流或实验热流。
if nargin<2, overrides=struct; end
spec=struct('kind',kind,'radius',0.6,'coneAngleRad',pi/12,'startMeridianM',0.2, ...
    'startPolarRad',0.35,'directionAngleRad',0.25,'circulationM2S',1200, ...
    'vortexRadiusM',1.3,'speedMS',750,'sphereSpeedScaleMS',900, ...
    'totalTemperatureK',1100,'totalPressurePa',200000,'origin',[0;0;0], ...
    'rotation',eye(3),'planeMap',eye(2),'coordinateScale',[1;1],'focusDistanceM',0.8);
for key=fieldnames(overrides)'
    name=key{1};assert(isfield(spec,name));spec.(name)=overrides.(name);
end
gas=struct('gamma',1.4,'R',287,'temperatureRangeK',[80 2000], ...
    'referenceTemperatureK',300,'referencePressurePa',101325);
if strcmp(kind,'vortex'),finish=1.28;step=0.04;unit='m';
elseif any(strcmp(kind,{'plane','unsupported','focus'})),finish=1.2;step=0.04;unit='m';
elseif startsWith(kind,'sphere'),finish=0.42;step=0.02;unit='rad';
else,finish=1;step=0.04;unit='rad';end
problem=struct('sample',@sample,'seed',@seed,'labelCenter',0,'labelStep',0.005, ...
    'labelUnit',unit,'arcStationsM',linspace(0,finish,round(finish/step)+1)', ...
    'gas',gas,'sources',struct('geometry',['解析验证曲面：' kind], ...
    'pressure','给定解析验证压力，与同一p0/T0完全气体场一致；非箭体压力求解', ...
    'thermodynamics','单一总焓/单一熵的理想气体解析验证输入', ...
    'entryDirection','解析入口方向；不由投影自由来流或Dijkstra选路'));

    function entry=seed(beta)
        angle=spec.directionAngleRad;
        if any(strcmp(kind,{'plane','unsupported'}))
            t=[-cos(angle);sin(angle);0];b=cross([0;0;1],t);
            q=spec.planeMap\(beta*b(1:2));direction=spec.rotation*t;
        elseif strcmp(kind,'vortex')
            q=spec.planeMap\[spec.vortexRadiusM-beta;0];direction=spec.rotation*[0;1;0];
        elseif strcmp(kind,'focus')
            q=spec.planeMap\[0;-beta];direction=spec.rotation*[-1;beta/spec.focusDistanceM;0];
        elseif strcmp(kind,'cylinder')
            theta=beta;q=[0;theta]./spec.coordinateScale;
            direction=spec.rotation*[-cos(angle);-sin(angle)*sin(theta);sin(angle)*cos(theta)];
        elseif strcmp(kind,'cone')
            theta=beta;q=[spec.startMeridianM;theta]./spec.coordinateScale;
            direction=spec.rotation*[-cos(spec.coneAngleRad);sin(spec.coneAngleRad)*cos(theta);sin(spec.coneAngleRad)*sin(theta)];
        elseif startsWith(kind,'sphere')
            phi=spec.startPolarRad;theta=beta;q=[phi;theta]./spec.coordinateScale;
            meridian=[-sin(phi);cos(phi)*cos(theta);cos(phi)*sin(theta)];azimuth=[0;-sin(theta);cos(theta)];
            if strcmp(kind,'sphere_sine'),angle=0;end
            direction=spec.rotation*(cos(angle)*meridian+sin(angle)*azimuth);
        else,error('SurfaceFlow:Benchmark','未知验证场。');end
        entry=struct('coordinate',q,'direction',direction);
    end

    function value=sample(q)
        scale=spec.coordinateScale;u=scale.*q;U=spec.speedMS;dU=zeros(2,1);
        if any(strcmp(kind,{'plane','unsupported','focus','vortex'}))
            xy=spec.planeMap*q;r=[xy;0];A=[spec.planeMap;zeros(1,2)];B=zeros(3,3);
            if strcmp(kind,'vortex')
                radius=norm(xy);assert(radius>0.2,'SurfaceFlow:ChartDomain','旋涡轴附近不在本验证范围。');
                U=spec.circulationM2S/radius;dU=spec.planeMap'*(-spec.circulationM2S*xy/radius^3);
            end
        elseif strcmp(kind,'cylinder')
            theta=u(2);R=spec.radius;r=[-u(1);R*cos(theta);R*sin(theta)];
            A=[[-1;0;0],[0;-R*sin(theta);R*cos(theta)]];
            B=[zeros(3,2),[0;-R*cos(theta);-R*sin(theta)]];
        elseif strcmp(kind,'cone')
            s=u(1);theta=u(2);delta=spec.coneAngleRad;
            assert(s>0.02,'SurfaceFlow:ChartDomain','锥尖坐标奇点须另行启动。');
            R=s*sin(delta);r=[-s*cos(delta);R*cos(theta);R*sin(theta)];
            A=[[-cos(delta);sin(delta)*cos(theta);sin(delta)*sin(theta)],[0;-R*sin(theta);R*cos(theta)]];
            B=[zeros(3,1),[0;-sin(delta)*sin(theta);sin(delta)*cos(theta)],[0;-R*cos(theta);-R*sin(theta)]];
        elseif startsWith(kind,'sphere')
            phi=u(1);theta=u(2);R=spec.radius;
            assert(phi>0.03&&phi<pi-0.03,'SurfaceFlow:ChartDomain','球面极点不在本坐标图内。');
            r=R*[cos(phi);sin(phi)*cos(theta);sin(phi)*sin(theta)];
            A=R*[[-sin(phi);cos(phi)*cos(theta);cos(phi)*sin(theta)],[0;-sin(phi)*sin(theta);sin(phi)*cos(theta)]];
            B=[-r,R*[0;-cos(phi)*sin(theta);cos(phi)*cos(theta)],R*[0;-sin(phi)*cos(theta);-sin(phi)*sin(theta)]];
            if strcmp(kind,'sphere_sine'),U=spec.sphereSpeedScaleMS*sin(phi);dU=[spec.sphereSpeedScaleMS*cos(phi)*scale(1);0];end
        else,error('SurfaceFlow:Benchmark','未知验证场。');end
        if ~any(strcmp(kind,{'plane','unsupported','focus','vortex'}))
            A=A.*scale';B=B.*[scale(1)^2,scale(1)*scale(2),scale(2)^2];
        end
        cp=gas.gamma*gas.R/(gas.gamma-1);T=spec.totalTemperatureK-U^2/(2*cp);dT=-U*dU/cp;
        p=spec.totalPressurePa*(T/spec.totalTemperatureK)^(gas.gamma/(gas.gamma-1));
        dp=gas.gamma/(gas.gamma-1)*p/T*dT;rho=p/(gas.R*T);drho=rho*(dp/p-dT/T);
        attached=~strcmp(kind,'unsupported')||q(1)>=-0.45;
        value=struct('position',spec.origin+spec.rotation*r,'firstDerivatives',spec.rotation*A, ...
            'secondDerivatives',spec.rotation*B,'pressure',p,'density',rho,'temperature',T,'speed',U, ...
            'pressureGradient',dp,'densityGradient',drho,'temperatureGradient',dT,'speedGradient',dU, ...
            'isAttached',logical(attached),'isSmooth',true);
    end
end
