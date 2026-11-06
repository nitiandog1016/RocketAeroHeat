function report=testAerothermalEvidenceProbe(out,holdSeconds)
%TESTAEROTHERMALEVIDENCEPROBE 只复核实际发现的探针可读性修正，不重算或重导大CSV。
% 模块一：读取同一已审核原面研究缓存，实例化主工作台的共享区域页。
% 模块二：核对绘制实际使用的固定色标、原面探针及周向交线，保留实际窗口观察。
if nargin<2,holdSeconds=0;end
assert(~isfolder(out));mkdir(out);root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
d=load(fullfile(root,'outputs/v167_validation/main01/result.mat'),'r');r=d.r;
d=load(fullfile(root,'outputs/v104_validation/baseline02/current_case.mat'),'current');
assert(isequal(r.V,d.current.V)&&isequal(r.F,d.current.F));r.geometryRegion=d.current.localGeometry.region;
fig=uifigure('Name','RocketAeroHeat · 面元查询复核 v176','Position',[40 50 1420 900]);
guard=onCleanup(@()closeIfPresent(fig));
v=createAerothermalEvidenceView(fig,buildSU2EvidenceViewerData(r));waitView(v);
v.selectField('heatFlux');v.setRange(-1e5,1e5);v.inspectFace(9327);waitView(v);
s=v.getState();assert(isequal(s.ranges.heatFlux(:),[-1e5;1e5]));
assert(isequal(s.rendered.geometry.fieldRange(:),[-1e5;1e5]));
assert(s.rendered.geometry.inspectedFace==9326 && s.selectedFaceID==9327);
assert(s.rendered.geometry.selectedFaceCount>0 && numel(v.getSection().segments)==92);
report=struct('passed',true,'originalFaces',19752,'actualRenderedRange',s.rendered.geometry.fieldRange, ...
    'selectedFaceID',s.selectedFaceID,'sectionSegments',92,'newCFDTasks',0,'repeatedExports',false, ...
    'nativeMouseObserved',false,'scope','探针折叠修改后的共享区域页；不重跑完整五步或已过数值/导出');
writeJSON(fullfile(out,'report.json'),report);fprintf('RAH_V176_PROBE_READY: %s\n',out);
t=tic;while toc(t)<holdSeconds && isgraphics(fig) && ~isfile(fullfile(out,'close_ui.request')),drawnow;pause(.1);end
report.observationWaitSeconds=toc(t);s=v.getState();report.finalSelectedFaceID=s.selectedFaceID;
report.finalRange=s.ranges.heatFlux;writeJSON(fullfile(out,'report.json'),report);
end
function waitView(view)
t=tic;while toc(t)<30
    drawnow;s=view.html.UserData;
    if ~isempty(s.error),error('RocketAeroHeat:EvidenceProbe:HTML','%s',s.error.message);end
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(view.getSection())&&s.rendered.geometry.visiblePixelSamples>0,break;end
    pause(.05);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision);
assert(s.rendered.geometry.faces==19752 && s.rendered.geometry.webglError==0 && s.rendered.geometry.visiblePixelSamples>0);
end
function writeJSON(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value));
end
function closeIfPresent(fig)
if isgraphics(fig),delete(fig);end
end
