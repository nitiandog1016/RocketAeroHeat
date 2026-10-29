function [problem,spec]=labelledSurfaceBenchmarkProblem(kind,input)
%LABELLEDSURFACEBENCHMARKPROBLEM 非均匀入口与真实内部斜激波锚点的曲面验证。
% 模块一：平面剪切/转向场核对H/S横向输运，输入压力明确是验证条件。
% 模块二：圆柱熵剪切层使用真正内部弱斜激波锚点；不把它拼成三维激波面。
% 模块三：入口斜置及坐标变换检验沿程标签偏移和物理横向梯度的区别。
if nargin<2,input=struct;end
spec=struct('kind',kind,'mode','thermallyPerfect','pressurePa',30000,'pressureExponentPerM',0, ...
    'totalTemperatureK',1200,'totalPressurePa',220000,'temperatureSlopeK',600, ...
    'temperatureCurvatureK',200,'logPressureSlope',2.5,'logPressureCurvature',.5, ...
    'directionAngleRad',0,'entrySkew',.35,'planeMap',eye(2),'rotation',eye(3), ...
    'radiusM',.1,'startAxialM',.45,'labelCenter',0,'labelStep',.004,'finishArcM',.2, ...
    'shockCenterAngleDeg',20,'shockVariationDeg',8);
for key=fieldnames(input)',name=key{1};assert(isfield(spec,name));spec.(name)=input.(name);end
gas=createSurfaceFlowGas(spec.mode,287.05,1.4);
unit='m';if strcmp(kind,'cylinderShock'),unit='rad';end
sourceBank={};bankLabels=[];
problem=struct('sample',@sample,'seed',@seed,'entryTotalState',@totalState, ...
    'entryRecords',@records,'labelCenter',spec.labelCenter,'labelStep',spec.labelStep, ...
    'labelUnit',unit,'arcStationsM',linspace(0,spec.finishArcM,21)','gas',gas, ...
    'sources',struct('geometry',['解析光滑验证曲面：' kind], ...
    'pressure','共同给定静压/连续指数压力；本模块未求真实箭体压力', ...
    'thermodynamics','逐入口标签固定H0/S0；同cp积分恢复，不沿程重置总压', ...
    'entryDirection','声明的光滑附着切向入口；斜置入口包含沿程标签偏移', ...
    'entryProfile','明确给定非均匀总状态的剪切/转向验证输入'));
if strcmp(kind,'cylinderShock')
    problem.sources.entryProfile='程序内部NASA7二维弱斜激波20+8sin(beta)度的独立锚点，显式转换约定；非同一三维激波面';
end

    function entry=seed(beta)
        if strcmp(kind,'cylinderShock')
            q=[spec.startAxialM+spec.entrySkew*beta;beta];direction=[-1;0;0];
        else
            t=[-cos(spec.directionAngleRad);sin(spec.directionAngleRad);0];b=cross([0;0;1],t);
            q=spec.planeMap\(beta*(b(1:2)+spec.entrySkew*t(1:2)));direction=t;
        end
        entry=struct('coordinate',q,'direction',spec.rotation*direction);
    end

    function value=sample(q)
        if strcmp(kind,'cylinderShock')
            u=q(1);theta=q(2);R=spec.radiusM;
            r=[-u;R*cos(theta);R*sin(theta)];
            A=[[-1;0;0],[0;-R*sin(theta);R*cos(theta)]];
            B=[zeros(3,2),[0;-R*cos(theta);-R*sin(theta)]];
            p=spec.pressurePa;dp=zeros(2,1);
        else
            xy=spec.planeMap*q;r=[xy;0];A=[spec.planeMap;zeros(1,2)];B=zeros(3,3);
            p=spec.pressurePa*exp(spec.pressureExponentPerM*xy(1));
            dp=spec.planeMap'*[spec.pressureExponentPerM*p;0];
        end
        value=struct('position',spec.rotation*r,'firstDerivatives',spec.rotation*A, ...
            'secondDerivatives',spec.rotation*B,'pressure',p,'pressureGradient',dp,'isSmooth',true,'isAttached',true);
    end

    function value=records()
        value=sourceBank;
    end

    function total=totalState(beta)
        if strcmp(kind,'cylinderShock')
            found=find(bankLabels==beta,1);
            if ~isempty(found),total=sourceBank{found}.totalState;return;end
            upstream=struct('Mach',5,'R',gas.R,'gamma',1.4,'temperature',220,'pressure',1200, ...
                'density',1200/(gas.R*220),'speed',5*sqrt(1.4*gas.R*220));
            angle=spec.shockCenterAngleDeg+spec.shockVariationDeg*sin(beta);
            shock=solveThermallyPerfectObliqueShock(upstream,angle);
            assert(shock.available&&shock.attached,'SurfaceFlow:ShockAnchor','入口二维弱斜激波未解，不能填值。');
            total=surfaceFlowTotalStateFromAnchor(shock,gas);
            record=struct('label',beta,'upstream',upstream,'deflectionDeg',angle,'shock',shock,'totalState',total);
        else
            T0=spec.totalTemperatureK+spec.temperatureSlopeK*beta+spec.temperatureCurvatureK*beta^2;
            p0=spec.totalPressurePa*exp(spec.logPressureSlope*beta+spec.logPressureCurvature*beta^2);
            total=struct('available',true,'totalPressure',p0,'totalTemperature',T0, ...
                'thermalPropertyMode',spec.mode,'thermalPropertyConvention',gas.thermalPropertyConvention, ...
                'method','逐标签给定的光滑总温/总压验证剖面，非箭体外流已求入口');
            record=struct('label',beta,'totalState',total);
        end
        bankLabels(end+1)=beta;sourceBank{end+1}=record;
    end
end
