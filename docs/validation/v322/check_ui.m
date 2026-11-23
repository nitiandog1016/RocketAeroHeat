function check_ui(out,kind)
% 模块一：通过普通科研侧栏验证真实端点与两套下载包。首个同刻走实际计算入口。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
if strcmp(kind,'compute'),modes={'wakeStart'};else,modes={'wakeCoarse','wakeFine','archive'};end
for item=modes
    mode=item{1};dest=fullfile(out,mode);mkdir(dest);c=fig.UserData.controls;
    c.shortSphereMode.Value=mode;c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
    v=fig.UserData.shortSphereReview;assert(~isempty(v) && strcmp(v.report.action,mode),'%s',fig.UserData.status.Text);
    s=v.report;d=v.web.evidenceView.getData();h=v.web.html;timer=tic;
    while isempty(h.UserData.rendered)&&toc(timer)<40,drawnow;pause(.1);end
    g=h.UserData.rendered;assert(g.webglError==0 && g.visiblePixelSamples>0 && g.faces==4464);
    if strcmp(mode,'archive')
        assert(~isfield(v,'wakeTab') && ~isfield(d.resultCondition,'elapsedPhysicalTimeS'));continue;
    end
    assert(isfield(v,'wakeTab') && s.wakeComparison.withinTolerance && strcmp(s.frameSHA256,v.mesh.evidenceIdentity.sourceFrame.sha256));
    assert(d.resultCondition.elapsedPhysicalTimeS==s.elapsedPhysicalTimeS && s.physicalStep*s.dtDU==s.timeDU);
    assert(contains(d.resultCondition.timeCaption,'扰动') && height(v.stations)==7 && s.coverage.fullSurfaceComputed);
    if strcmp(kind,'compute'),assert(s.timeDU==.25 && isfile(fullfile(s.jobDirectory,'result_002.json')));
    else,assert(s.timeDU==8 && isempty(s.jobDirectory));end
    v.tabs.SelectedTab=v.wakeTab;drawnow;assert(numel(v.wakeVelocity.Children.XData)==s.physicalStep+1);
    v.tabs.SelectedTab=v.web.tab;drawnow;
    % 模块二：核对同一物理层的探针、截线、原面热功率和测站导出。
    e=v.web.evidenceView;e.selectField('heatFlux');e.inspectFace(200,[]);u=e.getState();assert(u.selectedFaceID==200);
    x=.00137;e.setSectionX(x);timer=tic;
    while toc(timer)<40
        drawnow;cut=e.getSection();
        if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-12);
    e.exportSectionTo(fullfile(dest,'section.csv'));w=e.exportWallPairTo(fullfile(dest,'wall.csv'));
    assert(w.rows==4464 && w.joint.pressureHeatFluxQualifiedArea==0);v.exportTo(fullfile(dest,'stations.csv'));
    page=exportAerothermalWeb(v.mesh,fullfile(dest,'workbench'));
    proof.(mode)=struct('timeS',s.elapsedPhysicalTimeS,'timeDU',s.timeDU,'dtDU',s.dtDU,'frame',s.frameSHA256,'previous',s.previousFrameSHA256, ...
        'page',page,'job',s.jobDirectory,'powerW',s.wholeOutW,'pixels',g.visiblePixelSamples,'error',g.webglError);
end
assert(isempty(fig.UserData.workflowState.rawMesh));
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
