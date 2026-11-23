function check_resolution_ui(out)
% 模块一：按原侧栏打开两种步长，检查空间剖面和三维面元来自同一帧。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
proof=struct();
for entry={'wakeCoarse','wakeFine','archive'}
    mode=entry{1};c=fig.UserData.controls;
    c.shortSphereMode.Value=mode;c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
    v=fig.UserData.shortSphereReview;
    assert(~isempty(v) && strcmp(v.report.action,mode),'%s',fig.UserData.status.Text);
    if strcmp(mode,'archive')
        assert(~isfield(v,'resolution') && ~isfield(v.report,'resolution'));
        proof.steadyExcluded=true;continue;
    end
    d=v.report.resolution;page=v.resolution;
    assert(strcmp(d.frameSHA256,v.mesh.evidenceIdentity.sourceFrame.sha256) && d.timeS==v.report.elapsedPhysicalTimeS);
    v.tabs.SelectedTab=page.tab;drawnow;assert(numel(page.choice.Items)==6);
    for k=1:6
        page.choice.Value=k;page.choice.ValueChangedFcn(page.choice,[]);drawnow;
        section=d.sections(k);band=section.band;
        assert(isequal(page.curve.XData(:),section.radiusD(:)) && isequal(page.curve.YData(:),section.downstreamSpeedOverU(:)));
        assert(isequal(page.marks.XData,[band.r10D band.r90D]) && contains(page.note.Text,sprintf('%d 条独立射线',band.innerRays)));
        assert(contains(page.choice.Items{k},sprintf('%g',section.xD)));
    end
    % 模块二：切换剖面后回到表面页和测站导出，核对现有的面、时刻和七行数据。
    v.tabs.SelectedTab=v.web.tab;drawnow;e=v.web.evidenceView;
    e.selectField('heatFlux');e.inspectFace(2976,[]);state=e.getState();
    assert(state.selectedFaceID==2976 && v.report.coverage.fullSurfaceComputed);
    path=fullfile(out,[mode '.csv']);v.exportTo(path);rows=readtable(path,'TextType','string');
    assert(height(rows)==7 && all(rows.frameSHA256==string(d.frameSHA256)));
    assert(max(abs(rows.computedNu(end-1:end)-[d.stations.computedNu]'))<1e-10);
    proof.(mode)=struct('frame',d.frameSHA256,'seconds',d.timeS,'sections',numel(d.sections), ...
        'probeFace',state.selectedFaceID,'stationRows',height(rows));
end
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
