function check_surface_ui(out)
% 模块一：真实打开双支管归档，核对离线 WebGL 绘制、全原面与五区探针。
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supports';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'web'),'%s',fig.UserData.status.Text);
s=v.report;assert(s.coverage.fullSurfaceComputed && s.faces==32480 && s.jointExperimentalAreaM2==0);
html=v.web.html;t=tic;while isempty(html.UserData.rendered) && toc(t)<65,drawnow;pause(.15);end
gl=html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==s.faces);
e=v.web.evidenceView;e.selectField('heatFlux');
for j=1:numel(s.regions)
    id=s.regions(j).peakOutFace;e.inspectFace(id,[]);u=e.getState();assert(u.selectedFaceID==id);
end
% 模块二：在相同结果页面选择相对压力，核对真实截面与全部壁面对照导出。
e.selectField('pressureGauge');x=.00137;e.setSectionX(x);t=tic;
while toc(t)<55
    drawnow;cut=e.getSection();
    if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressureGauge')&&abs(cut.x-x)<1e-13,break;end
    pause(.12);
end
assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressureGauge')&&abs(cut.x-x)<1e-13);
e.exportSectionTo(fullfile(out,'section.csv'));w=e.exportWallPairTo(fullfile(out,'wall.csv'));
assert(w.rows==s.faces && w.joint.pressureHeatFluxQualifiedArea==0);
v.tabs.SelectedTab=v.surface.tab;drawnow;
assert(height(v.surface.table.Data)==5 && numel(v.surface.line.YData)==7);v.surface.exportTo(fullfile(out,'stations.csv'));
page=exportAerothermalWeb(v.mesh,fullfile(out,'workbench'));
assert(isfile(page) && isempty(fig.UserData.workflowState.rawMesh));
proof=struct('frame',s.frameSHA256,'page',page,'faces',s.faces,'areaM2',s.areaM2,'pixels',gl.visiblePixelSamples, ...
    'webglError',gl.webglError,'allRegionProbes',[s.regions.peakOutFace],'sectionX',x,'wallW',s.thermal.early.wallW);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);done=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
