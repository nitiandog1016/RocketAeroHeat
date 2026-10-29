function report=exportLabelledPressureSurfaceFlowStudy(evidenceFile,folder)
%EXPORTLABELLEDPRESSURESURFACEFLOWSTUDY 将非均匀入口及真实横向梯度接回统一工作台。
% 模块一：只读取已保存的八组核验结果，不重算流线或斜激波。
% 模块二：展示两个平面和完整锥柱底图，锥段/底面仍明确未求。
% 模块三：入口标签、总焓、熵、物理梯度与逐变量方法均进入探针和CSV。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
assert(~exist(folder,'dir'),'SurfaceFlow:StudyOverwrite','导出须使用新目录。');
saved=load(evidenceFile,'cases','report');
assert(strcmp(saved.report.status,'passed-nonuniform-entry-HS-physical-gradient-and-internal-shock-anchor-surface-flow') ...
    &&saved.report.cases==8,'SurfaceFlow:StudyInput','缺少通过的非均匀入口核验。');
groups={2,4,[7 8]};captions={'非均匀H/S：精确平面剪切','非均匀H/S：曲面Euler转向','内部斜激波锚点：圆柱熵剪切'};
ids={'labelled_NASA7_shear','labelled_NASA7_turn','internal_shock_entropy_cylinder'};
scopes={'给定共同静压和非均匀入口；精确稳态Euler剪切验证，不是箭体已求外流', ...
    '给定连续指数压力与非均匀入口；核对转向及物理横向梯度，非完整外流解', ...
    '独立二维弱斜激波锚点构造的圆柱熵剪切验证；不是20度整箭或同一三维激波面，锥段/底面未求'};
meshes=cell(3,1);overlays=meshes;viewCases=meshes;
for k=1:3
    selected=saved.cases(groups{k});first=selected{1}.result;
    if k<3
        meshes{k}=struct('V',[-.45 -.12 0;.1 -.12 0;.1 .2 0;-.45 .2 0],'F',[1 2 3;1 3 4]);
    else
        meshes{k}=buildCanonicalBody('coneCylinder',struct('radius',.1,'cylinderLength',.6,'coneAngleDeg',15,'resolution',32));
        meshes{k}.V(:,1)=meshes{k}.V(:,1)-max(meshes{k}.V(:,1));
    end
    tubes=cell(numel(selected),1);curves=tubes;drift=0;RHSCalls=0;stations=0;T0=[];p0=[];
    for j=1:numel(selected)
        r=selected{j}.result;c=r.center;g=r.physicalGradients;
        assert(all(r.validStation)&&~r.productionQualified&&~r.auxiliaryRayDerivativesArePhysical ...
            &&~r.entryTotalStateResetAlongRay&&~r.fullExternalContinuitySolved,'SurfaceFlow:StudyInput','身份或有效前缀不符合条件研究。');
        assert(all(isfinite(g.entropyGradient),'all')&&all(isfinite(g.totalEnthalpyGradient),'all'), ...
            'SurfaceFlow:StudyInput','不能把缺失或辅助延拓导数导出为物理梯度。');
        entry=r.entryTotalStates{3};
        tubes{j}=struct('labelCenter',r.labelCenter,'arcLengthM',r.arcLengthM, ...
            'position',c.position,'tangent',c.tangent,'metric',r.metric,'metricRatio',r.metricRatio, ...
            'metricDifferenceIndicator',r.metricDifferenceIndicator,'pressureInputPa',c.pressure, ...
            'temperatureRecoveredK',c.temperature,'densityRecoveredKgM3',c.density,'speedRecoveredMS',c.speed, ...
            'inletTotalTemperatureK',entry.totalTemperature,'inletTotalPressurePa',entry.totalPressure, ...
            'inletDerivativeResolution',r.entryDerivativeResolution, ...
            'totalEnthalpyJPerKg',c.totalEnthalpy,'entropyJPerKgK',c.entropy, ...
            'totalEnthalpyGradient',g.totalEnthalpyGradient,'entropyGradient',g.entropyGradient, ...
            'temperatureGradient',g.temperatureGradient,'entropyGradientNormJPerKgKM',vecnorm(g.entropyGradient,2,2));
        curves{j}=struct('sourceTubeIndexOneBased',j,'points',c.position,'normals',c.normal);
        for ray=1:5
            drift=max(drift,max(r.rays{ray}.unitTangentError));RHSCalls=RHSCalls+r.rays{ray}.RHSCalls;
            T0(end+1)=r.entryTotalStates{ray}.totalTemperature;p0(end+1)=r.entryTotalStates{ray}.totalPressure; %#ok<AGROW>
        end
        RHSCalls=RHSCalls+5;stations=stations+numel(r.arcLengthM);
    end
    sources=first.sources;sources.transverseGradients=[first.gradientMethod '；标签二阶差分，需收敛核对，非严格误差界'];
    sources.derivativeRole='单射线保存导数仅为固定标签辅助延拓；中心探针仅导出相邻标签重建的物理两向梯度';
    sources.thermodynamicConvention='NASA7原cp保持，1000K高段h/s积分常数连续匹配；每标签总状态只冻结一次';
    viewCases{k}=struct('id',ids{k},'caption',captions{k},'scope',scopes{k},'geometryScope',scopes{k}, ...
        'geometry',geometryAudit(meshes{k},k==3),'labelUnit',first.labelUnit,'metricUnit',first.metricUnit, ...
        'methodScope',[first.scope '；每标签H0/S0守恒，中心两向梯度由物理横向h恢复；Tr/Tw/热流未求'], ...
        'sources',sources,'maximumRelativeMetricDifference',first.options.maximumRelativeMetricDifference, ...
        'tubes',{tubes},'curves',{curves},'summary',struct('validStations',stations,'requestedStations',stations, ...
        'maximumUnitTangentError',drift,'RHSCalls',RHSCalls), ...
        'pressureSourceRole','给定共同静压/连续指数压力验证输入；内部斜激波只提供入口总状态锚点，未求整箭压力', ...
        'pressureIsInput',true,'thermodynamicStateIsInput',false,'thermodynamicStateIsConditionallyRestored',true, ...
        'inletVariesByLabel',true,'entryTotalTemperatureRangeK',[min(T0) max(T0)], ...
        'entryTotalPressureRangePa',[min(p0) max(p0)], ...
        'thermodynamicConvention','continuous-project-cp-v1','thermalConventionCaption','逐标签H/S；原cp及连续积分约定', ...
        'heatFluxSolved',false,'wallTemperatureSolved',false,'productionQualified',false,'newProductionSolvedAreaM2',0);
    overlays{k}=struct('kind','pressure-surface-flow-research','curves',{curves}, ...
        'thermodynamicStateIsConditionallyRestored',true,'initialCameraAzimuthRad',pi/2, ...
        'initialCameraElevationRad',.35,'productionQualified',false);
end
mkdir(folder);web=fullfile(folder,'workbench');geometryFolder=fullfile(web,'geometry');
exportRocketTestGeometryView(meshes,geometryFolder,captions,overlays);
copyfile(fullfile(root,'web/pressure-surface-geometry.js'),fullfile(geometryFolder,'rocket-geometry-bootstrap.js'));
stamp(fullfile(geometryFolder,'index.html'),geometryFolder,{'geometry.css','geometry.js','rocket-geometry-data.js','rocket-geometry-bootstrap.js'});
for name={'rocket-sensitivity.css','pressure-surface-study.css','pressure-surface-study.js'}
    copyfile(fullfile(root,'web',name{1}),fullfile(web,name{1}));
end
copyfile(fullfile(root,'web/pressure-surface-study.html'),fullfile(web,'index.html'));
data=struct('kind','pressure-surface-flow-study-v1','cases',{viewCases}, ...
    'geometryURL',['geometry/index.html?v=' fileSHA(fullfile(geometryFolder,'index.html'))], ...
    'sourceEvidenceFile',strrep(evidenceFile,'\','/'),'sourceEvidenceSHA256',fileSHA(evidenceFile), ...
    'productionQualified',false,'newProductionSolvedAreaM2',0);
text=jsonencode(data);writeUtf8(fullfile(web,'pressure-surface-flow-data.json'),text);
writeUtf8(fullfile(web,'pressure-surface-flow-data.js'),['window.PRESSURE_SURFACE_FLOW_STUDY=' text ';' newline]);
stamp(fullfile(web,'index.html'),web,{'rocket-sensitivity.css','pressure-surface-study.css','pressure-surface-study.js','pressure-surface-flow-data.js'});
report=struct('status','exported-saved-labelled-HS-and-physical-gradient-surface-flow-study','models',3, ...
    'conditionalTubes',4,'actualStations',sum(cellfun(@(c)c.summary.validStations,viewCases)), ...
    'sourceEvidenceSHA256',fileSHA(evidenceFile),'newSurfaceIntegrations',0,'newActualNSEvaluations',0, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'workbench',fullfile(web,'index.html'));
writeUtf8(fullfile(folder,'labelled_study_export.json'),jsonencode(report));disp(report);
end

function result=geometryAudit(mesh,closed)
% 模块用途：完整底图闭合/定向/体积单列，不把闭合几何当作热流已解面积。
V=mesh.V;F=mesh.F;e=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];[~,~,g]=unique(sort(e,2),'rows');
count=accumarray(g,1);orientation=accumarray(g,sign(e(:,2)-e(:,1)));
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
if closed,assert(all(count==2)&&all(orientation==0)&&volume>0);end
result=struct('vertices',size(V,1),'faces',size(F,1),'closed',closed, ...
    'boundaryEdges',nnz(count==1),'nonManifoldEdges',nnz(count>2),'signedVolumeM3',volume);
end
function stamp(file,folder,names)
html=fileread(file);
for j=1:numel(names)
    name=names{j};hash=fileSHA(fullfile(folder,name));pattern=[regexptranslate('escape',name) '(?:\?[^" ]*)?(?=[" ])'];
    html=regexprep(html,pattern,[name '?v=' hash(1:12)]);
end
writeUtf8(file,html);
end
function value=fileSHA(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),b=fread(fid,1024^2,'*uint8');if ~isempty(b),md.update(typecast(b,'int8'));end,end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
function writeUtf8(file,value)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
