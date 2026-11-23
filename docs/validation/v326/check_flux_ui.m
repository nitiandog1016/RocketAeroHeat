function check_flux_ui(out)
% 模块一：从统一工作台选择银球装置，在同一页检查三种结果身份。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supports';
c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'fluidGrid'),'%s',fig.UserData.status.Text);
page=v.fluidGrid;v.tabs.SelectedTab=page.tab;drawnow;
for entry={'structured','flux','hybrid'}
    selection=entry{1};page.method.Value=selection;page.method.ValueChangedFcn(page.method,[]);
    if strcmp(selection,'flux')
        grid=page.report.hybrid;f=page.flux.frames(end);q=page.flux.first;
        assert(strcmp(page.point.Visible,'on') && contains(page.stats.Text,num2str(q.worstID)));
        assert(contains(page.outcome.Text,page.flux.method) && contains(page.outcome.Text,sprintf('%.4g',f.minScalar)));
        assert(contains(page.outcome.Tooltip,sprintf('%.3g K',page.flux.scalarUnitK*f.minScalar)));
        assert(page.point.XData==q.xyzD(2) && abs(page.point.YData-hypot(q.xyzD(1),q.xyzD(3)))<1e-14);
        proof.projectedPointD=[page.point.XData page.point.YData];proof.firstUpdate=q.update;
    else
        grid=page.report.(selection);f=grid.stages(end).frames(end);
        assert(strcmp(page.point.Visible,'off') && isempty(page.outcome.Tooltip));
    end
    assert(contains(page.outcome.Text,'0 m²') && contains(page.outcome.Text,sprintf('%d 个节点低于',f.below)));
    % 模块二：对照实际截面边，检查近壁、外域与装置视图切换。
    for level={'joint','outer','body'}
        page.zoom.Value=level{1};page.zoom.ValueChangedFcn(page.zoom,[]);drawnow;
        source=grid.section.pointsD;edges=grid.section.edges;
        assert(isequal(page.line.XData(1:3:end)',source(edges(:,1),1)));
        assert(isequal(page.line.YData(2:3:end)',source(edges(:,2),2)));
    end
    proof.(selection)=struct('mesh',grid.meshSHA256,'frame',f.sha256,'label',page.outcome.Text);
end
v.tabs.SelectedTab=v.geometryTab;drawnow;
assert(height(v.geometry.stations)==7 && isempty(fig.UserData.workflowState.rawMesh));
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);
closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
