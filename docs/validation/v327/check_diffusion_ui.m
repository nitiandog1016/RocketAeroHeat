function check_diffusion_ui(out)
% 模块一：通过统一侧栏打开双支管档案，从对流回放切换到扩散精度页。
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supports';
c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
review=fig.UserData.shortSphereReview;page=review.fluidGrid;review.tabs.SelectedTab=page.tab;
select('flux');assert(strcmp(page.point.Visible,'on') && strcmp(page.edgeMark.Visible,'off'));
select('diffusion');d=page.diffusion;w=d.star.worst;
assert(~d.passed && strcmp(page.edgeMark.Visible,'on') && strcmp(page.point.Visible,'off'));
assert(contains(page.stats.Text,num2str(w.edge)) && contains(page.stats.Text,sprintf('%.2f%%',100*w.normalizedError)));
assert(contains(page.outcome.Text,'解析检验未通过') && contains(page.outcome.Text,'0 m²'));
assert(contains(page.outcome.Tooltip,sprintf('%.6g',w.exact)) && contains(page.outcome.Tooltip,sprintf('%.6g',w.limited)));
% 模块二：用归档端点检查实际画线，确认切换范围后仍为同一条边。
ends=squeeze(d.star.segmentEndsD(d.star.edgeIDs==w.edge,:,:));ref=page.report.hybrid;
for level={'outer','joint','body'}
    page.zoom.Value=level{1};page.zoom.ValueChangedFcn(page.zoom,[]);drawnow;
    assert(isequal(page.edgeMark.XData(:),ends(:,2)) && ...
        all(abs(page.edgeMark.YData(:)-hypot(ends(:,1),ends(:,3)))<1e-14));
    assert(isequal(page.line.XData(1:3:end)',ref.section.pointsD(ref.section.edges(:,1),1)));
end
proof=struct('edge',w.edge,'endpointsD',ends,'meshSHA256',d.meshSHA256,'label',page.outcome.Text, ...
    'details',page.outcome.Tooltip,'stats',page.stats.Text,'acceptedAreaM2',page.flux.acceptedAreaM2);
for mode={'structured','hybrid','flux'}
    select(mode{1});assert(strcmp(page.edgeMark.Visible,'off'));
end
review.tabs.SelectedTab=review.geometryTab;drawnow;assert(height(review.geometry.stations)==7);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);
closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
    function select(mode)
        page.method.Value=mode;page.method.ValueChangedFcn(page.method,[]);drawnow;
    end
end
