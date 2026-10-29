function report=exportThermalPressureSurfaceFlowStudy(evidenceFile,folder)
%EXPORTTHERMALPRESSURESURFACEFLOWSTUDY 将已通过的同物性流线结果接入统一工作台。
% 模块一：只读取实际保存的核验结果，不为看图重新积分或覆盖旧MAT。
% 模块二：两组平面验证与完整锥柱底图分开；只叠加真实中心流管站点。
% 模块三：压力输入、条件反算Te/rho/U、未知Tr/Tw/热流各自保留身份。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
assert(~exist(folder,'dir'),'SurfaceFlow:StudyOverwrite','导出必须使用新目录。');
saved=load(evidenceFile,'cases','report');
assert(isfield(saved,'report')&&strcmp(saved.report.status, ...
    'passed-pressure-state-adapter-variable-cp-and-internal-cone-surface-flow'), ...
    'SurfaceFlow:StudyInput','只接受已通过的实际核验结果。');
cases=saved.cases;plane=find(cellfun(@(c)strcmp(c.kind,'NASA7_exponential_pressure_plane'),cases));
cone=find(cellfun(@(c)strcmp(c.kind,'internal_NASA7_Taylor_Maccoll_cone'),cases));
assert(numel(plane)==2&&numel(cone)==2,'SurfaceFlow:StudyInput','缺少对应的实际保存案例。');
groups={plane(1),plane(2),cone};ids={'NASA7_hot_plane','NASA7_cold_plane','NASA7_internal_cone'};
captions={'NASA7变压平面：跨1000K','低温变压平面：跨150K','内部NASA7切锥：零攻角15°'};
scopes={'解析连续压力验证片；不是完整箭体外流，检验跨1000K的条件反算与动量输运', ...
    '解析连续压力验证片；80–150K为原振动冻结近似，检验跨150K的条件反算', ...
    '完整闭合锥柱底图；锥段压力来自实际内部Taylor–Maccoll，仅两组保存流管，柱段/底面未求'};
meshes=cell(3,1);overlays=cell(3,1);viewCases=cell(3,1);
for k=1:3
    selected=cases(groups{k});first=selected{1}.result;coupling=first.thermalPressureCoupling;
    if k<3
        meshes{k}=struct('V',[-.65 -.08 0;.05 -.08 0;.05 .23 0;-.65 .23 0],'F',[1 2 3;1 3 4]);
        pressureRole='连续指数压力解析验证输入；未求真实箭体压力';
    else
        meshes{k}=buildCanonicalBody('coneCylinder',struct('radius',.1,'cylinderLength',.6, ...
            'coneAngleDeg',15,'resolution',32));
        meshes{k}.V(:,1)=meshes{k}.V(:,1)-max(meshes{k}.V(:,1));
        pressureRole='程序内部M5零攻角15° NASA7 Taylor–Maccoll已算锥面压力；作为流线输入';
    end
    geometry=geometryAudit(meshes{k},k==3);tubes=cell(numel(selected),1);curves=tubes;
    drift=0;RHSCalls=0;stations=0;
    for j=1:numel(selected)
        result=selected{j}.result;center=result.center;
        assert(all(result.validStation)&&~result.productionQualified&& ...
            strcmp(result.gas.thermalPropertyConvention,'continuous-project-cp-v1'), ...
            'SurfaceFlow:StudyInput','流管未完成声明范围或物性约定错误。');
        assert(isequal(result.thermalPressureCoupling.totalPressure,coupling.totalPressure)&& ...
            isequal(result.thermalPressureCoupling.totalTemperature,coupling.totalTemperature));
        tubes{j}=struct('labelCenter',result.labelCenter,'arcLengthM',result.arcLengthM, ...
            'position',center.position,'tangent',center.tangent,'metric',result.metric, ...
            'metricRatio',result.metricRatio,'metricDifferenceIndicator',result.metricDifferenceIndicator, ...
            'pressureInputPa',center.pressure,'temperatureRecoveredK',center.temperature, ...
            'densityRecoveredKgM3',center.density,'speedRecoveredMS',center.speed);
        curves{j}=struct('sourceTubeIndexOneBased',j,'points',center.position,'normals',center.normal);
        for ray=1:numel(result.rays)
            drift=max(drift,max(result.rays{ray}.unitTangentError));RHSCalls=RHSCalls+result.rays{ray}.RHSCalls;
        end
        RHSCalls=RHSCalls+result.entryValidationRHSCalls;stations=stations+numel(result.arcLengthM);
    end
    sources=first.sources;sources.thermodynamicConvention='NASA7原cp逐值保留，1000K高段焓/熵常数连续匹配；不修改生产物性';
    if k==3
        converted=selected{1}.convertedTotalState;
        sources.anchorConversion=sprintf('原物理锥面p/T/rho/U不变；总压约定差%.12gPa，总温约定差%.12gK，非沿程重置', ...
            converted.totalPressureConventionDifference,converted.totalTemperatureConventionDifference);
    end
    viewCases{k}=struct('id',ids{k},'caption',captions{k},'scope',scopes{k}, ...
        'geometryScope',scopes{k},'geometry',geometry,'labelUnit',first.metricUnit(3:end), ...
        'metricUnit',first.metricUnit,'methodScope',[first.scope '；Te为固定总状态与压力的同物性条件反算，未求Tr/Tw/热流'], ...
        'sources',sources,'totalTemperatureK',coupling.totalTemperature,'totalPressurePa',coupling.totalPressure, ...
        'maximumRelativeMetricDifference',first.options.maximumRelativeMetricDifference, ...
        'tubes',{tubes},'curves',{curves},'summary',struct('validStations',stations, ...
        'requestedStations',stations,'maximumUnitTangentError',drift,'RHSCalls',RHSCalls), ...
        'pressureSourceRole',pressureRole,'pressureIsInput',true,'thermodynamicStateIsInput',false, ...
        'thermodynamicStateIsConditionallyRestored',true,'thermodynamicConvention','continuous-project-cp-v1', ...
        'thermalConventionCaption','原NASA7/低温冻结cp，1000K处h/s连续匹配', ...
        'heatFluxSolved',false,'wallTemperatureSolved',false,'productionQualified',false,'newProductionSolvedAreaM2',0);
    overlays{k}=struct('kind','pressure-surface-flow-research','curves',{curves}, ...
        'thermodynamicStateIsConditionallyRestored',true, ...
        'initialCameraAzimuthRad',pi/2,'initialCameraElevationRad',.35,'productionQualified',false);
end
mkdir(folder);workbench=fullfile(folder,'workbench');geometryFolder=fullfile(workbench,'geometry');
exportRocketTestGeometryView(meshes,geometryFolder,captions,overlays);
copyfile(fullfile(root,'web','pressure-surface-geometry.js'),fullfile(geometryFolder,'rocket-geometry-bootstrap.js'));
stamp(fullfile(geometryFolder,'index.html'),geometryFolder,{'geometry.css','geometry.js','rocket-geometry-data.js','rocket-geometry-bootstrap.js'});
for name={'rocket-sensitivity.css','pressure-surface-study.css','pressure-surface-study.js'}
    copyfile(fullfile(root,'web',name{1}),fullfile(workbench,name{1}));
end
copyfile(fullfile(root,'web','pressure-surface-study.html'),fullfile(workbench,'index.html'));
data=struct('kind','pressure-surface-flow-study-v1','cases',{viewCases}, ...
    'geometryURL',['geometry/index.html?v=' fileSHA(fullfile(geometryFolder,'index.html'))], ...
    'sourceEvidenceFile',strrep(evidenceFile,'\','/'),'sourceEvidenceSHA256',fileSHA(evidenceFile), ...
    'productionQualified',false,'newProductionSolvedAreaM2',0);
text=jsonencode(data);writeUtf8(fullfile(workbench,'pressure-surface-flow-data.json'),text);
writeUtf8(fullfile(workbench,'pressure-surface-flow-data.js'),['window.PRESSURE_SURFACE_FLOW_STUDY=' text ';' newline]);
stamp(fullfile(workbench,'index.html'),workbench,{'rocket-sensitivity.css','pressure-surface-study.css','pressure-surface-study.js','pressure-surface-flow-data.js'});
report=struct('status','exported-saved-conditional-NASA7-surface-flow-and-internal-cone-study', ...
    'models',3,'conditionalTubes',sum(cellfun(@numel,groups)),'actualStations',sum(cellfun(@(c)c.summary.validStations,viewCases)), ...
    'sourceEvidenceSHA256',fileSHA(evidenceFile),'newSurfaceIntegrations',0,'newActualNSEvaluations',0, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'workbench',fullfile(workbench,'index.html'));
writeUtf8(fullfile(folder,'thermal_study_export.json'),jsonencode(report));disp(report);
end

function result=geometryAudit(mesh,closed)
% 模块用途：完整箭体逐边/定向/正体积核对；平面验证片明确另列。
V=mesh.V;F=mesh.F;edge=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];[~,~,group]=unique(sort(edge,2),'rows');
count=accumarray(group,1);orientation=accumarray(group,sign(edge(:,2)-edge(:,1)));
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
if closed,assert(all(count==2)&&all(orientation==0)&&volume>0);end
result=struct('vertices',size(V,1),'faces',size(F,1),'closed',closed, ...
    'boundaryEdges',nnz(count==1),'nonManifoldEdges',nnz(count>2),'signedVolumeM3',volume);
end
function stamp(file,folder,names)
% 模块用途：资源URL绑定实际字节，避免旧缓存或误改JSON下载路径。
html=fileread(file);
for k=1:numel(names)
    name=names{k};digest=fileSHA(fullfile(folder,name));pattern=[regexptranslate('escape',name) '(?:\?[^" ]*)?(?=[" ])'];
    html=regexprep(html,pattern,[name '?v=' digest(1:12)]);
end
writeUtf8(file,html);
end
function digest=fileSHA(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),bytes=fread(fid,1024^2,'*uint8');if ~isempty(bytes),md.update(typecast(bytes,'int8'));end,end
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
function writeUtf8(file,value)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
