function report=testAerothermalEvidenceWorkbench(out,holdSeconds)
%TESTAEROTHERMALEVIDENCEWORKBENCH 用既存原面结果验收新增区域页，不重新运行CFD。
% 模块一：工程与研究各沿用准确缓存；只将同V/F的几何标签用于区域查看。
% 模块二：实例化主工作台使用的相同组件，核对原值、选区、探针、截面与导出。
% 模块三：验证更新清理旧回执，最后允许实际鼠标查看；不冒称重跑旧五步数值链。
if nargin<2,holdSeconds=0;end
assert(~isfolder(out));mkdir(out);root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
engineeringFile=fullfile(root,'outputs/v104_validation/baseline02/current_case.mat');
researchFile=fullfile(root,'outputs/v167_validation/main01/result.mat');
d=load(engineeringFile,'current');current=d.current;d=load(researchFile,'r');r=d.r;
assert(size(current.F,1)==19752 && isequal(current.V,r.V) && isequal(current.F,r.F));
assert(r.receipt.exitCode==0 && r.thermalStateAudit.interfacePassed && ...
    r.convergenceAudit.checkpointIdentityConfirmed && ~r.productionQualified);
assert(r.convergenceAudit.checkpointIteration==199 && r.convergenceAudit.completedUpdatesBeforeRun==200);
current.evidenceIdentity=struct('origin','engineering','sourceTraceable',true, ...
    'sourceFrame',struct('kind','saved-engineering-surface','path',engineeringFile,'sha256',fileSHA(engineeringFile)));
r.geometryRegion=current.localGeometry.region;
r.geometryRegionSource='与研究V/F逐项相同的工程原面几何标签；没有转移工程流动资格';

fig=uifigure('Name','RocketAeroHeat · 原面区域结果复核 v176','Position',[40 50 1420 900]);
guard=onCleanup(@()closeIfPresent(fig));
layout=uigridlayout(fig,[2 1]);layout.RowHeight={30,'1x'};layout.Padding=[4 4 4 4];
note=uilabel(layout,'Text','已保存原面结果复核 · 没有重新求解 · 工程值与研究值分开 · 物理验证资格另列');
note.Layout.Row=1;tabs=uitabgroup(layout);tabs.Layout.Row=2;
engineeringTab=uitab(tabs,'Title','工程原面结果');researchTab=uitab(tabs,'Title','内嵌原面研究');
engineeringTabs=uitabgroup(engineeringTab,'Units','normalized','Position',[0 0 1 1]);
engineeringUI=createAerothermalWebTab(engineeringTabs,current);
engineeringTabs.SelectedTab=engineeringUI.evidenceTab;tabs.SelectedTab=engineeringTab;
ev=engineeringUI.evidenceView;waitView(ev,19752);ed=ev.getData();
assert(isequal(ed.faces,current.F-1) && isequal(ed.vertices,current.V));
assert(max(abs(ed.faceArea-r.faceArea))<1e-12);
assert(abs(sum(ed.faceArea)-1.9646292565715551)<1e-11);
assert(ed.evidenceLedger.joint.pressureHeatFluxQualifiedArea==0);
ep=field(ed,'pressure');eq=field(ed,'heatFlux');
assert(isequaln(ep.values,current.pressure.staticPressure(:)/1000));
assert(isequaln(eq.values,current.heating.heatFlux(:)/1000));
assert(ep.areaSummary.prescribedArea>0 && eq.areaSummary.engineeringArea>0);
assert(eq.areaSummary.researchArea==0 && eq.areaSummary.qualifiedArea==0);
ev.selectField('pressure');waitView(ev,19752);ev.exportTo(fullfile(out,'engineering_pressure.csv'));
ev.selectField('heatFlux');ev.setRange(-100,100);waitView(ev,19752);
region=ed.evidenceLedger.regions(end);ev.selectRegion(region.name);waitView(ev,19752);
selected=ev.getState();assert(strcmp(selected.regionName,region.name));
assert(isequal(selected.ranges.heatFlux,[-100 100]));
assert(selected.rendered.geometry.selectedFaceCount==numel(region.faceIDs));
ev.inspectFace(12000);ev.setSectionX(-1.5);waitView(ev,19752);
selected=ev.getState();assert(selected.selectedFaceID==12000);
assert(selected.rendered.geometry.inspectedFace==11999);
section=ev.getSection();assert(numel(section.segments)==92);
ev.exportTo(fullfile(out,'engineering_heatFlux.csv'));ev.exportSectionTo(fullfile(out,'engineering_section.csv'));
writeOracleData(fullfile(out,'engineering_data.json'),ed);
% 工况更新与旧帧反例只作用于测试页面；原缓存和物理值保持不变。
old=ev.getState();oldMessage=old.rendered;changed=current;
changed.evidenceIdentity.sourceFrame.kind='same-values-new-review-source';
updateAerothermalWebTab(engineeringUI,changed);waitView(ev,19752);
fresh=ev.getState();assert(fresh.revision>old.revision && fresh.selectedFaceID==0 && isempty(fresh.regionName));
ev.html.Data=oldMessage;ev.html.DataChangedFcn(ev.html,[]);after=ev.getState();
assert(after.revision==fresh.revision && after.selectedFaceID==0,'旧帧不能恢复旧选择');
updateAerothermalWebTab(engineeringUI,current);waitView(ev,19752);

researchUI=createSU2ResearchView(researchTab,r);tabs.SelectedTab=researchTab;
researchUI.tabs.SelectedTab=researchUI.evidenceTab;rv=researchUI.evidenceView;waitView(rv,19752);
rd=rv.getData();rp=field(rd,'pressure');rq=field(rd,'heatFlux');
assert(isequal(rd.vertices,r.V) && isequal(rd.faces,r.F-1));
assert(isequaln(rp.values,r.pressureRaw) && isequaln(rq.values,r.heatFluxRaw));
assert(rp.areaSummary.engineeringArea==0 && rq.areaSummary.engineeringArea==0);
assert(rq.areaSummary.researchArea>0 && rq.areaSummary.qualifiedArea==0);
assert(rd.evidenceLedger.joint.pressureHeatFluxQualifiedArea==0);
rv.selectField('pressure');waitView(rv,19752);rv.exportTo(fullfile(out,'research_pressure.csv'));
rv.selectField('heatFlux');rv.setRange(-1e5,1e5);rv.selectRegion(rd.evidenceLedger.regions(end).name);waitView(rv,19752);
rv.inspectFace(12000);rv.setSectionX(-1.5);waitView(rv,19752);
selected=rv.getState();section=rv.getSection();
assert(isequal(selected.ranges.heatFlux,[-1e5 1e5]));
assert(numel(section.segments)==92);
rv.exportTo(fullfile(out,'research_heatFlux.csv'));rv.exportSectionTo(fullfile(out,'research_section.csv'));
rv.selectField('wallSolved');waitView(rv,19752);unknown=field(rv.getData(),'wallSolved');
assert(all(isnan(unknown.values)) && all(unknown.availability==0));
rv.exportTo(fullfile(out,'research_wallSolved.csv'));
rv.selectField('wallInput');waitView(rv,19752);given=field(rv.getData(),'wallInput');
assert(all(given.availability==3) && given.areaSummary.qualifiedArea==0);
writeOracleData(fullfile(out,'research_data.json'),rd);
rv.selectField('heatFlux');rv.setRange(-1e5,1e5);rv.selectRegion(rd.evidenceLedger.regions(end).name);waitView(rv,19752);
report=struct('passed',true,'componentIntegration',true,'fullWorkflowNumericalRerun',false, ...
    'originalFaces',19752,'totalArea',sum(ed.faceArea),'engineeringHeatArea',eq.areaSummary.engineeringArea, ...
    'engineeringPressurePrescribedArea',ep.areaSummary.prescribedArea,'researchHeatArea',rq.areaSummary.researchArea, ...
    'qualifiedPressureHeatIntersection',0,'sectionSegments',92,'sameOriginalGeometry',true, ...
    'unchangedPhysicalValues',true,'regionSelectionRendered',true,'probeRendered',true, ...
    'fixedRange',true,'unknownSolvedTw',true,'givenTwNotQualified',true, ...
    'updatedSourceClearedSelection',true,'staleMessageRejected',true,'newCFDTasks',0, ...
    'engineeringSourceSHA256',fileSHA(engineeringFile),'researchSourceSHA256',fileSHA(researchFile), ...
    'nativeMouseObserved',false);
assert(ancestor(ev.html,'figure')==fig && ancestor(rv.html,'figure')==fig);
report.singleWindowHost=true;
writeJSON(fullfile(out,'report.json'),report);
fprintf('RAH_V176_UI_READY: %s\n',out);
% 实际UI观察等待与求解/测试计算分列；关闭信号只控制本次创建的窗口。
observed=tic;closeFile=fullfile(out,'close_ui.request');
while toc(observed)<holdSeconds && isgraphics(fig) && ~isfile(closeFile),drawnow;pause(.1);end
report.observationWaitSeconds=toc(observed);writeJSON(fullfile(out,'report.json'),report);
end

function f=field(data,key)
f=data.evidenceLedger.fields(strcmp({data.evidenceLedger.fields.key},key));assert(isscalar(f));
end
function waitView(view,n)
t=tic;while toc(t)<35
    drawnow;s=view.html.UserData;
    if ~isempty(s.error),error('RocketAeroHeat:EvidenceTest:HTML','%s',s.error.message);end
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(view.getSection()),break;end
    pause(.05);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision,'区域页面没有当前绘制回执');
assert(s.rendered.geometry.faces==n && s.rendered.geometry.webglError==0 && s.rendered.geometry.visiblePixelSamples>0);
end
function writeJSON(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,'ConvertInfAndNaN',true));
end
function writeOracleData(path,data)
% 独立复核只需一份原几何/原场；去掉重复副本，重复原因采用字典无损编码。
l=data.evidenceLedger;
names=intersect(fieldnames(l),{'vertices','faces','faceArea','region'});
if ~isempty(names),l=rmfield(l,names);end
copy={'values','status','legacyStatus','methodNames','methodIndex','scope','scopeNames','scopeIndex','missingReason','history'};
names=intersect(fieldnames(l.fields),copy);if ~isempty(names),l.fields=rmfield(l.fields,names);end
for i=1:numel(l.fields)
    [reasonNames,~,index]=unique(l.fields(i).reason,'stable');
    l.fields(i).reasonNames=reasonNames;l.fields(i).reasonIndex=index-1;
end
l.fields=rmfield(l.fields,'reason');data.evidenceLedger=l;
writeJSON(path,data);
end
function value=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function closeIfPresent(fig)
if isgraphics(fig),delete(fig);end
end
