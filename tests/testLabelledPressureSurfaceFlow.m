function report=testLabelledPressureSurfaceFlow(folder)
%TESTLABELLEDPRESSURESURFACEFLOW 非均匀H/S、物理两向梯度及内部斜激波入口验证。
% 模块一：明确区分单条射线辅助导数与真实相邻标签梯度，原v95文件不修改。
% 模块二：检查斜置入口、两组间距收敛、转向及坐标变换，保留全部实际输出。
% 模块三：实际调用内部20度附近的弱斜激波生成熵差；不宣称整箭外流或热流已求。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupSurfacePressureFlow;
assert(~exist(folder,'dir'),'SurfaceFlow:TestOverwrite','必须使用新目录。');mkdir(folder);
cases={};checks={};rejections={};

%% 一、平面精确剪切解：H/S均随横向位置变化，入口同时存在沿程偏移
for mode={'caloricallyPerfect','thermallyPerfect'}
    [p,spec]=labelledSurfaceBenchmarkProblem('plane',struct('mode',mode{1}));
    r=marchLabelledPressureSurfaceFlow(p,struct('maximumArcStepM',.01));
    assert(all(r.validStation)&&max(abs(r.metric-1))<1e-11);
    assert(max(abs(r.labelAlongStreamShift-spec.entrySkew))<1e-11);
    assert(max(abs(r.fullLabelSeparationNorm-sqrt(1+spec.entrySkew^2)))<1e-11);
    assert(abs(r.entryEnthalpyLabelDerivative)>1e5&&abs(r.entryEntropyLabelDerivative)>1);
    assert(r.entryDerivativeResolution.enthalpyVariationResolved&&r.entryDerivativeResolution.entropyVariationResolved);
    assert(max(vecnorm(r.physicalGradients.entropyGradient,2,2))>1);
    assert(max(abs(r.physicalGradientRHSIdentityResidual))<1e-10);
    assert(max(abs(r.physicalGradients.temperatureGradient(:,1)))<1e-10);
    assert(max(abs(r.physicalGradients.temperatureGradient(:,2)))>1);
    for j=1:5
        assert(max(abs(r.rays{j}.totalEnthalpy-r.entryTotalEnthalpy(j)))<2e-6);
        assert(max(abs(r.rays{j}.entropy-r.entryEntropy(j)))<2e-9);
    end
    cases{end+1}=struct('kind',['shear_' mode{1}],'specification',spec,'result',r,'entryRecords',{p.entryRecords()}); %#ok<AGROW>
end

%% 二、非均匀总状态的真实转向；标签半距给出二阶横向导数收敛
for step=[.008 .004 .002]
    [p,spec]=labelledSurfaceBenchmarkProblem('plane',struct('pressureExponentPerM',1, ...
        'directionAngleRad',.3,'labelStep',step));
    r=marchLabelledPressureSurfaceFlow(p,struct('maximumArcStepM',.01));
    assert(all(r.validStation)&&max(abs(r.physicalGradientRHSIdentityResidual))<1e-8);
    assert(norm(r.center.tangent(end,:)-r.center.tangent(1,:))>1e-3);
    assert(max(vecnorm(r.physicalGradients.entropyGradient,2,2))>1);
    cases{end+1}=struct('kind','nonuniform_turn_plane','specification',spec,'result',r,'entryRecords',{p.entryRecords()}); %#ok<AGROW>
end

%% 三、斜坐标和刚体旋转后的同一物理梯度，不能使用坐标导数替代物理梯度
angle=.61;rotation=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
[p,spec]=labelledSurfaceBenchmarkProblem('plane',struct('planeMap',[2 .3;.1 .6], ...
    'rotation',rotation,'pressureExponentPerM',1,'directionAngleRad',.3,'labelStep',.004));
r=marchLabelledPressureSurfaceFlow(p,struct('maximumArcStepM',.01));reference=cases{4}.result;
assert(all(r.validStation));
assert(max(abs(r.center.position-reference.center.position*rotation'),[],'all')<1e-11);
assert(max(abs(r.physicalGradients.entropyGradient-reference.physicalGradients.entropyGradient*rotation'),[],'all')<2e-8);
cases{end+1}=struct('kind','transformed_nonuniform_turn','specification',spec,'result',r,'entryRecords',{p.entryRecords()});

%% 四、真实内部20度附近弱斜激波锚点，映射到圆柱轴向熵剪切层
for center=[0 pi/4]
    [p,spec]=labelledSurfaceBenchmarkProblem('cylinderShock',struct('labelCenter',center,'labelStep',.004,'entrySkew',.03));
    r=marchLabelledPressureSurfaceFlow(p,struct('maximumArcStepM',.01));
    assert(all(r.validStation)&&max(abs(r.metric-spec.radiusM))<3e-7);
    assert(abs(r.entryEntropyLabelDerivative)>1&&max(vecnorm(r.physicalGradients.entropyGradient,2,2))>10);
    assert(~r.entryDerivativeResolution.enthalpyVariationResolved&&r.entryDerivativeResolution.entropyVariationResolved);
    assert(max(abs(r.physicalGradientRHSIdentityResidual))<1e-8);
    records=p.entryRecords();assert(numel(records)==5);
    for j=1:5
        assert(records{j}.shock.available&&records{j}.shock.attached);
        assert(~records{j}.totalState.sourcePhysicalStateChanged);
    end
    cases{end+1}=struct('kind','internal_oblique_shock_entropy_cylinder','specification',spec,'result',r,'entryRecords',{records}); %#ok<AGROW>
end

%% 五、均匀入口退化到原v95；不是替代非均匀剖面测试
[p,spec]=labelledSurfaceBenchmarkProblem('plane',struct('temperatureSlopeK',0,'temperatureCurvatureK',0, ...
    'logPressureSlope',0,'logPressureCurvature',0,'pressureExponentPerM',1,'directionAngleRad',.3));
r=marchLabelledPressureSurfaceFlow(p,struct('maximumArcStepM',.01));
old=bindIsentropicSurfacePressure(p,p.entryTotalState(0));old=marchPressureDrivenSurfaceFlow(old,struct('maximumArcStepM',.01));
for name={'position','tangent','pressure','temperature','density','speed','totalEnthalpy','entropy'}
    assert(isequaln(r.center.(name{1}),old.center.(name{1})),name{1});
end
assert(isequaln(r.metric,old.metric));
checks{end+1}=struct('uniformLimitExactlyEqualsV95',true);

%% 六、来源、旧总状态约定、法向入口和失效前缀必须明确拒绝
bad=p;bad.sources=rmfield(bad.sources,'entryProfile');
rejections{end+1}=refused(@()marchLabelledPressureSurfaceFlow(bad),'SurfaceFlow:Source');
bad=p;rawTotal=p.entryTotalState(0);rawTotal=rmfield(rawTotal,'thermalPropertyConvention');bad.entryTotalState=@(~)rawTotal;
rejections{end+1}=refused(@()marchLabelledPressureSurfaceFlow(bad),'SurfaceFlow:ThermalConvention');
bad=p;bad.seed=@(~)struct('coordinate',[0;0],'direction',[0;0;1]);
rejections{end+1}=refused(@()marchLabelledPressureSurfaceFlow(bad),'SurfaceFlow:Seed');
bad=p;baseTotal=p.entryTotalState(0);bad.entryTotalState=@(beta)jumpTotal(beta,baseTotal);
rejections{end+1}=refused(@()marchLabelledPressureSurfaceFlow(bad),'SurfaceFlow:EntryDerivativeResolution');
[prefixProblem,prefixSpec]=labelledSurfaceBenchmarkProblem('plane');baseSample=prefixProblem.sample;
prefixProblem.sample=@(q)limited(baseSample(q),q);
prefix=marchLabelledPressureSurfaceFlow(prefixProblem,struct('maximumArcStepM',.01));
assert(any(prefix.validStation)&&any(~prefix.validStation));first=find(~prefix.validStation,1);assert(all(~prefix.validStation(first:end)));
assert(all(isnan(prefix.center.temperature(first:end)))&&all(isnan(prefix.metric(first:end))));
assert(all(isnan(prefix.physicalGradients.entropyGradient(first:end,:)),'all'));
report=struct('status','passed-nonuniform-entry-HS-physical-gradient-and-internal-shock-anchor-surface-flow', ...
    'cases',numel(cases),'checks',{checks},'rejections',{rejections},'allNumericalResultsAreConditional',true, ...
    'auxiliaryDerivativesNotPresentedAsPhysical',true,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'labelled_pressure_surface_flow_evidence.mat'),'cases','report','prefix','prefixSpec','-v7');
fid=fopen(fullfile(folder,'labelled_surface_flow_report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(report);
end

function sample=limited(sample,q)
sample.isAttached=q(1)>-.085;
end
function total=jumpTotal(beta,total)
total.totalPressure=total.totalPressure*exp(.02*(beta>0));
total.method='有意构造未分辨入口熵跳跃，必须拒绝而不能当作光滑横向梯度';
end
function name=refused(action,expected)
name='';try,action();catch ex,name=ex.identifier;end
assert(strcmp(name,expected),'SurfaceFlow:TestRejection','实际拒绝%s，期望%s。',name,expected);
end
