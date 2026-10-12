function report=diagnoseSphericalNoseInterface(inputFolder,outputFolder)
%DIAGNOSESPHERICALNOSEINTERFACE 用已保存算例检查球鼻模型的接口跳变与热流敏感性。
% 模块作用：仅在本函数的试验副本中替换球鼻外缘，再跑原流线/热流。
% 不修改默认求解器。覆盖率增长不等于接口正确或热流精度通过。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end
% 在同一解析位置比较两种状态，排除相邻离散面的几何步长影响。
% 球鼻 theta=60 度对应局部切楔转角 30 度；当前诊断固定 M4/20 km。
evalc('referenceCfg=createCaseConfig();');
referenceCfg.flight.Mach=4;referenceCfg.flight.altitude=20000;
referenceCfg.flight.alpha=20;referenceCfg.flight.beta=0;
referenceFree=calcFreestream(referenceCfg.flight);
normal=solveThermallyPerfectNormalShock(referenceFree);
wedge=solveThermallyPerfectObliqueShock(referenceFree,30);
assert(normal.available && wedge.available);
p=referenceFree.pressure+(normal.stagnationPressure-referenceFree.pressure)*cosd(60)^2;
state=recoverThermallyPerfectEdgeState(p,normal.stagnationPressure, ...
    normal.stagnationTemperature,referenceFree.R);
exact=struct('polarAngleDeg',60,'spherePressure',p,'wedgePressure',wedge.pressure, ...
    'sphereSpeed',state.speed,'wedgeSpeed',wedge.speed, ...
    'normalTotalPressure',normal.stagnationPressure,'wedgeTotalPressure',wedge.totalPressure, ...
    'speedJumpOverFreestream',(wedge.speed-state.speed)/referenceFree.speed);
writeJSON(fullfile(outputFolder,'exact_interface.json'),exact);
report=struct([]);
for kind={'hemisphereCylinder','bluntedCone'}
    for resolution=[24 48]
        name=sprintf('%s_n%d_a20',kind{1},resolution);
        saved=load(fullfile(inputFolder,[name '.mat']),'mesh','cfg');
        mesh=saved.mesh;cfg=saved.cfg;fs=calcFreestream(cfg.flight);
        sphere=solveSphericalNoseEdge(mesh,fs,cfg);
        row=struct('name',name,'available',sphere.available, ...
            'sphereFaces',nnz(sphere.algebraicallyValid), ...
            'baselineConnected',nnz(mesh.surfaceStreamlines.connectedFace), ...
            'trialConnected',NaN,'baselineHeatArea',mesh.heatingCoverage.summary.computedAreaFraction, ...
            'trialHeatArea',NaN,'baselineSeparationOnsets',nnz(mesh.heating.separationTrace.onsetFace), ...
            'trialSeparationOnsets',NaN,'interfaceEdges',0, ...
            'maxInterfaceSpeedStepOverFreestream',NaN,'maxInterfacePressureStepOverStagnation',NaN, ...
            'sameFacePressureDiscrepancyMax',NaN,'sameFaceSpeedDiscrepancyMax',NaN, ...
            'status',sphere.status);
        if sphere.available
            trial=mesh;e=mesh.boundaryLayerEdge;mask=sphere.algebraicallyValid;
            old=e;both=mask & old.algebraicallyValid;
            if any(both)
                row.sameFacePressureDiscrepancyMax=max(abs(sphere.pressure(both)./old.pressure(both)-1));
                row.sameFaceSpeedDiscrepancyMax=max(abs(sphere.speed(both)./old.speed(both)-1));
            end
            for key={'pressure','temperature','density','speed','Mach','speedOfSound', ...
                    'totalPressure','totalTemperature','velocityVector','incidenceNormal'}
                e.(key{1})(mask,:)=sphere.(key{1})(mask,:);
            end
            e.dynamicViscosity(mask)=1.458e-6*e.temperature(mask).^(3/2)./(e.temperature(mask)+110.4);
            e.algebraicallyValid(mask)=true;e.physicallyUnresolved(mask)=false;
            e.faceMethodID(mask)=7;e.faceMethodNames{8}=sphere.method;
            e.thermalPropertyModeByFace(mask)=2;e.thermalExpansionOutOfRangeFace(mask)=false;
            e.shockAngleDeg(mask)=NaN;e.incidenceAngleDeg(mask)=90-sphere.polarAngleDeg(mask);
            e.method=[old.method '；试验球鼻替换（未验证接口）：' sphere.method];
            e.scope=[old.scope '；' sphere.scope];
            trial.boundaryLayerEdge=e;trial.externalFlow.edge=e;trial.externalFlow.method=e.method;
            trial=traceSurfaceStreamlines(trial,cfg);trial=calcSurfaceHeating(trial,fs,cfg);
            row.trialConnected=nnz(trial.surfaceStreamlines.connectedFace);
            row.trialHeatArea=trial.heatingCoverage.summary.computedAreaFraction;
            row.trialSeparationOnsets=nnz(trial.heating.separationTrace.onsetFace);
            F=mesh.F;n=size(F,1);
            [~,~,group]=unique(sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2),'rows');
            counts=accumarray(group,1);faces=repmat((1:n)',3,1);
            a=accumarray(group,faces,[],@min);b=accumarray(group,faces,[],@max);
            crossing=counts==2 & xor(mask(a),mask(b)) & ...
                e.algebraicallyValid(a) & e.algebraicallyValid(b);
            a=a(crossing);b=b(crossing);row.interfaceEdges=numel(a);
            row.maxInterfaceSpeedStepOverFreestream=maxOrNaN(abs(e.speed(a)-e.speed(b))/fs.speed);
            row.maxInterfacePressureStepOverStagnation=maxOrNaN(abs(e.pressure(a)-e.pressure(b))/sphere.normalShock.stagnationPressure);
            % 不导出带临时混合来源的 mesh，避免试验字段进入正式查看器。
            save(fullfile(outputFolder,[name '_sphere_diagnostic.mat']),'sphere','row');
        end
        if isempty(report),report=row;else,report(end+1)=row;end %#ok<AGROW>
        fprintf('%s: sphere %d, connected %d -> %g, heat area %.2f -> %.2f%%, separation %d -> %g, max interface dU/Uinf %.4f.\n', ...
            name,row.sphereFaces,row.baselineConnected,row.trialConnected,100*row.baselineHeatArea, ...
            100*row.trialHeatArea,row.baselineSeparationOnsets,row.trialSeparationOnsets,row.maxInterfaceSpeedStepOverFreestream);
    end
end
writeJSON(fullfile(outputFolder,'sphere_interface.json'),report);
fprintf('Exact same-position interface: dU/Uinf %.4f, sphere/wedge pressure %.2f/%.2f Pa.\n', ...
    exact.speedJumpOverFreestream,exact.spherePressure,exact.wedgePressure);
end
function value=maxOrNaN(values)
% 无有效接口时保持缺失，不能用零假装跳变已通过。
value=NaN;if ~isempty(values),value=max(values);end
end
function writeJSON(path,data)
% 每个报告单独关闭文件；NaN 按 JSON null 保存，避免虚构零误差。
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'无法写入接口诊断报告');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(data,'PrettyPrint',true));
end
