function check_time_ui(out,branch)
% 模块一：从实际侧栏载入目标物理时间层，在三维页与能量账核对时间、面数和身份。
fig=rocketAeroHeatWorkbench();done=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value=['supportTime' upper(branch(1)) branch(2:end)];
c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'timePair'),'%s',fig.UserData.status.Text);
html=v.web.html;timer=tic;while isempty(html.UserData.rendered) && toc(timer)<60,drawnow;pause(.12);end
gl=html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==32480);
v.tabs.SelectedTab=v.timePair.tab;drawnow;
assert(height(v.timePair.ledger.Data)==3 && height(v.timePair.stations.Data)==7);
p=v.mesh.shortSupportTime;e=v.report.energy;
assert(strcmp(e.frame,p.([branch 'Frame'])) && e.elapsedS==p.elapsedS && abs(v.report.outW-v.report.coverage.wholeSurfaceHeatPowerW*(-1))<1e-10);
assert(max(abs(v.timePair.ledger.Data.storage/1e6-[p.fullEnergy.storageW;p.firstHalfEnergy.storageW;p.halfEnergy.storageW]))<1e-15);
% 模块二：检查真实读到的面号、压力热流、反向区与灰色拒判面，导出同刻截面。
ev=v.web.evidenceView;v.tabs.SelectedTab=v.web.evidenceTab;data=ev.getData();assert(strcmp(data.defaultField,'heatFlux'));
assert(strcmp(data.evidenceIdentity.sourceFrame.sha256,e.frame) && data.resultCondition.elapsedPhysicalTimeS==e.elapsedS);
key=['supportTime' branch 'Status'];ff=data.fields(strcmp({data.fields.key},key));
unknown=find(contains(ff.missingReason,'棱边'),1);assert(isnan(ff.values(unknown)));
for field={'pressureGauge','heatFlux',key}
    ev.selectField(field{1});
    for id=[v.report.peakOutFace,p.halfStations(6).faceID,p.halfStations(7).faceID,unknown]
        ev.inspectFace(id,[]);state=ev.getState();assert(state.selectedFaceID==id && strcmp(state.fieldKey,field{1}));
    end
end
ev.selectField('heatFlux');x=.00137;ev.setSectionX(x);timer=tic;
while toc(timer)<45
    drawnow;cut=ev.getSection();
    if isstruct(cut)&&isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-13,break;end
    pause(.1);
end
assert(isstruct(cut)&&isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-13);
ev.exportSectionTo(fullfile(out,'section.csv'));page=exportAerothermalWeb(v.mesh,fullfile(out,'workbench'));
proof=struct('frame',e.frame,'oldFrame',e.old,'elapsedS',e.elapsedS,'dtDU',e.dt,'faces',gl.faces,'pixels',gl.visiblePixelSamples, ...
    'wallW',e.wallW,'storageW',e.storageW,'page',page);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
