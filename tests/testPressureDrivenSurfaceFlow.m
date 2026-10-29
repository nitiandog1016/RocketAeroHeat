function report=testPressureDrivenSurfaceFlow(folder)
%TESTPRESSUREDRIVENSURFACEFLOW 验证压力转向、曲面联络、流管横向尺度与失败边界。
% 模块一：用直线、自由旋涡、柱/锥/球面解析结果核对，分开弧长与标签误差。
% 模块二：同一物理问题更换坐标尺度、仿射图及空间姿态，不能改变物理解路径。
% 模块三：拒绝坏EOS/导数/驻点/入口，保留附着区前缀并阻止聚焦后假结果。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupSurfacePressureFlow;
assert(~exist(folder,'dir'));mkdir(folder);
kinds={'plane','vortex','cylinder','cone','sphere_sine','sphere_geodesic'};cases=cell(numel(kinds),1);
for k=1:numel(kinds)
    [p,spec]=pressureSurfaceBenchmarkProblem(kinds{k});
    answer=marchPressureDrivenSurfaceFlow(p,struct('maximumArcStepM',0.005));
    [position,tangent,h]=reference(spec,p.arcStationsM);
    positionError=max(vecnorm(answer.center.position-position,2,2));
    tangentError=max(vecnorm(answer.center.tangent-tangent,2,2));
    metricError=max(abs(answer.metric-h))/max(h);
    assert(all(answer.validStation)&&positionError<2e-8&&tangentError<2e-8&&metricError<2e-4);
    assert(~answer.productionQualified&&~answer.physicalStateClipped&&answer.newProductionSolvedAreaM2==0);
    cases{k}=struct('kind',kinds{k},'specification',spec,'result',answer, ...
        'referencePosition',position,'referenceTangent',tangent,'referenceMetric',h, ...
        'positionErrorM',positionError,'tangentError',tangentError,'relativeMetricError',metricError);
end

%% 独立控制两类误差：RK4步长与相邻标签中心差分
rk=cell(3,1);label=cell(3,1);steps=[0.04 0.02 0.01];labels=[0.02 0.01 0.005];
for j=1:3
    [p,spec]=pressureSurfaceBenchmarkProblem('vortex');p.arcStationsM=[0;1.28];p.labelStep=0.001;
    a=marchPressureDrivenSurfaceFlow(p,struct('maximumArcStepM',steps(j)));
    expected=reference(spec,p.arcStationsM);errorValue=norm(a.center.position(end,:)-expected(end,:));
    rk{j}=struct('stepM',steps(j),'result',a,'positionErrorM',errorValue);
    [p,spec]=pressureSurfaceBenchmarkProblem('cone');p.labelStep=labels(j);
    a=marchPressureDrivenSurfaceFlow(p,struct('maximumArcStepM',0.02));[~,~,expectedH]=reference(spec,p.arcStationsM);
    label{j}=struct('labelStep',labels(j),'result',a,'relativeMetricError',max(abs(a.metric-expectedH))/max(expectedH));
end
rkError=cellfun(@(x)x.positionErrorM,rk);labelError=cellfun(@(x)x.relativeMetricError,label);
assert(all(rkError(1:2)>8*rkError(2:3))&&all(rkError(1:2)<24*rkError(2:3)));
assert(all(labelError(1:2)>3.9*labelError(2:3))&&all(labelError(1:2)<4.1*labelError(2:3)));

%% 重参数化及刚体姿态不变性
angle=0.4;rotation=[cos(angle) 0 sin(angle);0 1 0;-sin(angle) 0 cos(angle)];
[p,spec]=pressureSurfaceBenchmarkProblem('vortex',struct('planeMap',[2 0.35;0 0.5], ...
    'rotation',rotation,'origin',[-3;1;2]));
mapped=marchPressureDrivenSurfaceFlow(p,struct('maximumArcStepM',0.005));
base=cases{2}.result;expected=(rotation*base.center.position')'+spec.origin';
assert(max(vecnorm(mapped.center.position-expected,2,2))<2e-10);
assert(max(abs(mapped.metric-base.metric))<2e-9);
[p,spec]=pressureSurfaceBenchmarkProblem('sphere_geodesic',struct('coordinateScale',[0.7;2.5]));
scaled=marchPressureDrivenSurfaceFlow(p,struct('maximumArcStepM',0.005));base=cases{6}.result;
assert(max(vecnorm(scaled.center.position-base.center.position,2,2))<2e-10);
assert(max(abs(scaled.metric-base.metric))<2e-9);

%% 拒绝与前缀：无效输入不能变成零热流或新的区域覆盖
refusals={};[p,~]=pressureSurfaceBenchmarkProblem('plane');
for mode={'EOS','EOSGradient','Momentum','EnergyGradient','Stagnation','Temperature','Unsupported'}
    b=p;source=p.sample;b.sample=@(q)corrupt(source(q),mode{1});
    identifier=expectReject(b,struct);assert(strcmp(identifier,['SurfaceFlow:' mode{1}]));refusals{end+1}=identifier; %#ok<AGROW>
end
b=p;b.sources=rmfield(b.sources,'pressure');assert(strcmp(expectReject(b,struct),'SurfaceFlow:Source'));
b=p;source=p.seed;b.seed=@(beta)badSeed(source(beta));assert(strcmp(expectReject(b,struct),'SurfaceFlow:Seed'));
assert(strcmp(expectReject(p,struct('notAnOption',1)),'SurfaceFlow:Options'));
[p,~]=pressureSurfaceBenchmarkProblem('unsupported',struct('directionAngleRad',0));
prefix=marchPressureDrivenSurfaceFlow(p,struct('maximumArcStepM',0.01));
assert(any(prefix.validStation)&&any(~prefix.validStation));first=find(~prefix.validStation,1);
assert(all(prefix.validStation(1:first-1))&&all(~prefix.validStation(first:end)));
assert(all(isnan(prefix.metric(first:end)))&&all(isnan(prefix.center.pressure(first:end))));
assert(all(cellfun(@(r)strcmp(r.failureIdentifier,'SurfaceFlow:Unsupported'),prefix.rays)));
[p,~]=pressureSurfaceBenchmarkProblem('focus');p.arcStationsM=(0:0.04:1.2)';
focus=marchPressureDrivenSurfaceFlow(p,struct('maximumArcStepM',0.01));
assert(any(~focus.validStation));first=find(~focus.validStation,1);
assert(focus.arcLengthM(first)<=0.8+1e-12,'SurfaceFlow:FocusResolution','不能把解析聚焦点的差分残宽标为有效。');
assert(all(isnan(focus.metric(first:end)))&&all(isnan(focus.center.tangent(first:end,:)),'all'));
report=struct('status','passed-pressure-driven-surface-rays-analytic-geometry-coordinate-and-refusal-checks', ...
    'cases',numel(cases),'RKPositionErrorsM',rkError,'labelRelativeMetricErrors',labelError, ...
    'refusals',{refusals},'pressureAndThermalFieldsAreInputs',true, ...
    'fullExternalFlowOrHeatingSolved',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'pressure_surface_flow_evidence.mat'),'cases','rk','label','mapped','scaled','prefix','focus','report','-v7');disp(report);
end

function [position,tangent,metric]=reference(s,arc)
% 模块用途：世界坐标解析解，不调用被测曲面联络或压力RHS。
if strcmp(s.kind,'plane')
    t=[-cos(s.directionAngleRad),sin(s.directionAngleRad),0];position=arc*t;tangent=repmat(t,numel(arc),1);metric=ones(size(arc));
elseif strcmp(s.kind,'vortex')
    angle=arc/s.vortexRadiusM;position=s.vortexRadiusM*[cos(angle),sin(angle),zeros(size(arc))];
    tangent=[-sin(angle),cos(angle),zeros(size(arc))];metric=ones(size(arc));
elseif strcmp(s.kind,'cylinder')
    angle=s.directionAngleRad;theta=arc*sin(angle)/s.radius;
    position=[-arc*cos(angle),s.radius*cos(theta),s.radius*sin(theta)];
    tangent=[-cos(angle)*ones(size(arc)),-sin(angle)*sin(theta),sin(angle)*cos(theta)];metric=s.radius*cos(angle)*ones(size(arc));
elseif strcmp(s.kind,'cone')
    distance=s.startMeridianM+arc;position=[-distance*cos(s.coneAngleRad),distance*sin(s.coneAngleRad),zeros(size(arc))];
    tangent=repmat([-cos(s.coneAngleRad),sin(s.coneAngleRad),0],numel(arc),1);metric=distance*sin(s.coneAngleRad);
elseif strcmp(s.kind,'sphere_sine')
    phi=s.startPolarRad+arc/s.radius;position=s.radius*[cos(phi),sin(phi),zeros(size(arc))];
    tangent=[-sin(phi),cos(phi),zeros(size(arc))];metric=s.radius*sin(phi);
else
    phi=s.startPolarRad;psi=s.directionAngleRad;n0=[cos(phi),sin(phi),0];t0=cos(psi)*[-sin(phi),cos(phi),0]+sin(psi)*[0,0,1];
    angle=arc/s.radius;position=s.radius*(cos(angle)*n0+sin(angle)*t0);tangent=-sin(angle)*n0+cos(angle)*t0;
    normal=position/s.radius;J=[zeros(size(arc)),-position(:,3),position(:,2)];metric=sum(cross(normal,tangent,2).*J,2);
end
position=(s.rotation*position')'+s.origin';tangent=(s.rotation*tangent')';
end

function value=corrupt(value,mode)
% 模块用途：独立制造不相容输入，验证不能静默接受。
switch mode
    case 'EOS',value.pressure=value.pressure*1.1;
    case 'EOSGradient',value.pressureGradient(1)=value.pressureGradient(1)+1;
    case 'Momentum',value.speedGradient(1)=1;
    case 'EnergyGradient',value.temperatureGradient(1)=1;value.densityGradient=-value.density/value.temperature*value.temperatureGradient;
    case 'Stagnation',value.speed=0;
    case 'Temperature',value.temperature=50;
    case 'Unsupported',value.isSmooth=false;
end
end
function seed=badSeed(seed),seed.direction=[0;0;1];end
function identifier=expectReject(problem,options)
identifier='';try,marchPressureDrivenSurfaceFlow(problem,options);catch ex,identifier=ex.identifier;end
assert(~isempty(identifier),'SurfaceFlow:TestRefusal','错误输入未拒绝。');
end
