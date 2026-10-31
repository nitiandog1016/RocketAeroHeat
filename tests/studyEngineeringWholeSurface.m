function studyEngineeringWholeSurface(out)
%STUDYENGINEERINGWHOLESURFACE 重算真实20度箭体当前工况，绑定区域缺口。
% 模块一：只载入小型原几何/自由来流输入，不继承旧压力、热流或矩阵缓存。
% 模块二：执行与主程序一致的完整单工况工程链，报告真实阶段。
% 模块三：保存一份当前结果和逐面元审计；未解和给定状态保持原身份。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(root,'outputs','v104_validation','baseline01');end
assert(~exist(out,'dir'),'研究输出须使用新目录，不能改写历史证据。');mkdir(out);
input=fullfile(root,'tests','data','shock_layer_entry', ...
    'actual_rocket_geometry_and_freestream.mat');
loaded=load(input,'result','fs');mesh=loaded.result;
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
% 小型输入为节省空间可存整数；GUI/计算入口的角度统一转为double。
cfg.flight.Mach=double(loaded.fs.Mach);cfg.flight.altitude=double(loaded.fs.altitude);
cfg.flight.alpha=double(loaded.fs.alpha);cfg.flight.beta=double(loaded.fs.beta);cfg.flight.Twall=300;
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
[current,fs,report]=runAeroHeatSingleCase(mesh,cfg,@progress);
assert(isequal(current.F,mesh.F) && isequal(current.V,mesh.V));
assert(fs.Mach==5 && fs.alpha==20 && fs.altitude==30000);
coverage=current.heatingCoverage;area=double(current.localGeometry.area(:));
summary=report;summary=rmfield(summary,'cfg');
summary.heatingCoverage=coverage.summary;
summary.reasonNames=coverage.reasonNames;summary.areaByReason=coverage.areaByReason;
summary.countByReason=coverage.countByReason;
summary.sourceInput=input;summary.thermalModel=cfg.heating.airThermalModel;
summary.heatMethod=current.heating.method;
summary.completePhysicalSolution=false;
fid=fopen(fullfile(out,'summary.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(summary,'ConvertInfAndNaN',true));fclose(fid);
rows=table((1:size(current.F,1))',current.localGeometry.center(:,1), ...
    current.localGeometry.center(:,2),current.localGeometry.center(:,3),area, ...
    current.pressure.staticPressure,current.pressure.fieldStatus, ...
    current.boundaryLayerEdge.temperature,current.boundaryLayerEdge.speed, ...
    current.heating.heatFlux,current.heating.recoveryTemperature, ...
    current.heating.wallTemperatureByFace,current.heating.faceMethodID, ...
    coverage.reasonCode,'VariableNames',{'face','x_m','y_m','z_m','area_m2', ...
    'pressure_Pa','pressure_status','edge_temperature_K','edge_speed_m_s', ...
    'heat_flux_W_m2','recovery_temperature_K','prescribed_wall_temperature_K', ...
    'heat_method_id','heat_reason_code'});
writetable(rows,fullfile(out,'faces.csv'));
save(fullfile(out,'current_case.mat'),'current','fs','cfg','report','-v7');
fprintf('Current engineering baseline: %d faces, heat area %.12f%%, p computed %.9f%%, p prescribed %.9f%%.\n', ...
    size(current.F,1),100*report.heatComputedAreaFraction, ...
    100*report.pressureComputedAreaFraction,100*report.pressurePrescribedAreaFraction);
for j=1:numel(coverage.reasonNames)
    fprintf('Reason %d: %d faces; %.12f m2; %s\n', ...
        j-1,coverage.countByReason(j),coverage.areaByReason(j),coverage.reasonNames{j});
end
end

function progress(event)
% 模块四：打印阶段实际完成事件，不把执行进度当成物理覆盖率。
if strcmp(event.status,'completed')
    fprintf('ENGINEERING_PHASE %s %.3fs\n',event.phase,event.elapsedSeconds);
end
end
