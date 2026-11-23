function check_sector_ui(out)
% 模块一：从实际侧栏进入周向加密银球，检查三维渲染、七站和五区。
fig=rocketAeroHeatWorkbench();done=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supportDirectionAngular';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'comparison'),'%s',fig.UserData.status.Text);
html=v.web.html;t=tic;while isempty(html.UserData.rendered) && toc(t)<60,drawnow;pause(.12);end
gl=html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==32480);
p=v.mesh.shortSupportAngular;v.tabs.SelectedTab=v.comparison.tab;drawnow;
assert(height(v.comparison.regions.Data)==5 && height(v.comparison.stations.Data)==7);
assert(isequal(v.comparison.stations.Data.nu128,p.fineNu(:)) && isequal(v.comparison.stations.Data.nu64,p.coarseNu(:)));
% 模块二：同一原面上核对峰值、背风点、拒判原因与截面，随后导出供独立核对的离线页。
e=v.web.evidenceView;v.tabs.SelectedTab=v.web.evidenceTab;
d=e.getData();assert(strcmp(d.defaultField,'heatFlux') && strcmp(d.evidenceIdentity.sourceFrame.sha256,p.fineFrame));
assert(d.wallCoverage.jointComputedFaceCount==32480 && abs(d.wallCoverage.wholeSurfaceHeatPowerW+p.fineW)<1e-10);
ids=[25055,p.finePeak.faceID,p.direction.stations(6).faceID,p.direction.stations(7).faceID];
for key={'heatFlux','pressureGauge','supportAngularShearPa','supportAngularSpeedMps','supportAngularStatus'}
    e.selectField(key{1});
    for j=ids,e.inspectFace(j,[]);state=e.getState();assert(state.selectedFaceID==j && strcmp(state.fieldKey,key{1}));end
end
field=d.fields(strcmp({d.fields.key},'supportAngularStatus'));
j=find(contains(field.missingReason,'棱边'),1);assert(isnan(field.values(j)));
e.inspectFace(j,[]);state=e.getState();assert(state.selectedFaceID==j);
e.selectField('heatFlux');x=.00137;e.setSectionX(x);t=tic;
while toc(t)<45
    drawnow;cut=e.getSection();
    if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-13,break;end
    pause(.1);
end
assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-13);
e.exportSectionTo(fullfile(out,'section.csv'));
page=exportAerothermalWeb(v.mesh,fullfile(out,'workbench'));
proof=struct('faces',gl.faces,'pixels',gl.visiblePixelSamples,'oldPeak',ids(1),'peak',ids(2),'stations',7, ...
    'regions',5,'frame',p.fineFrame,'jointComputedAreaM2',d.wallCoverage.jointComputedAreaM2,'sectionX',x,'page',page);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
