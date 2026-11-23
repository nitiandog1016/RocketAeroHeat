function check_ui(out)
% 模块一：用实际侧栏切换归档和冷态新算，检查两个结果都具有各自的来源绑定。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
for action={'archive','solve'}
    mode=action{1};folder=fullfile(out,mode);mkdir(folder);
    c=fig.UserData.controls;c.shortSphereMode.Value=mode;c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
    v=fig.UserData.shortSphereReview;assert(~isempty(v) && strcmp(v.report.action,mode),'%s',fig.UserData.status.Text);
    s=v.report;stations=v.stations;assert(height(stations)==7 && s.coverage.fullSurfaceComputed && s.zeroVerifiedAreaM2==0);
    assert(size(v.table.Data,1)==7 && abs(v.table.Data(6,5)-100*stations.relativeError(6))<.001);
    html=v.web.html;timer=tic;while isempty(html.UserData.rendered) && toc(timer)<45,drawnow;pause(.1);end
    g=html.UserData.rendered;assert(g.webglError==0 && g.visiblePixelSamples>0 && g.faces==s.faces);
    % 模块二：检查前后半球探针、真实截线、全表面 CSV 与独立测站导出。
    e=v.web.evidenceView;e.selectField('heatFlux');[~,j]=max(v.mesh.shortSphere.angleDeg);e.inspectFace(j,[]);
    u=e.getState();assert(u.selectedFaceID==j);
    e.selectField('pressureGauge');[~,j]=min(v.mesh.shortSphere.angleDeg);e.inspectFace(j,[]);u=e.getState();assert(u.selectedFaceID==j);
    x=.00137;e.setSectionX(x);timer=tic;
    while toc(timer)<45
        drawnow;cut=e.getSection();
        if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressureGauge')&&abs(cut.x-x)<1e-13,break;end
        pause(.1);
    end
    assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressureGauge')&&abs(cut.x-x)<1e-13);
    e.exportSectionTo(fullfile(folder,'section.csv'));w=e.exportWallPairTo(fullfile(folder,'wall.csv'));assert(w.rows==s.faces && w.joint.pressureHeatFluxQualifiedArea==0);
    v.exportTo(fullfile(folder,'stations.csv'));page=exportAerothermalWeb(v.mesh,fullfile(folder,'workbench'));
    proof.(mode)=struct('frame',s.frameSHA256,'reaction',s.reactionSHA256,'page',page,'job',s.jobDirectory,'faces',s.faces, ...
        'wholeOutW',s.wholeOutW,'stationRMSE',s.stationRMSE,'glPixels',g.visiblePixelSamples,'glError',g.webglError);
end
assert(~isempty(proof.solve.job) && isfile(fullfile(proof.solve.job,'completed.json')) && isempty(proof.archive.job));
assert(isempty(fig.UserData.workflowState.rawMesh)); % 独立对照试件保持独立状态。
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
