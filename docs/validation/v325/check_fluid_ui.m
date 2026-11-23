function check_fluid_ui(out)
% 模块一：从统一侧栏选择双支管，实际切换两套流体体网和近壁视图。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supports';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'fluidGrid'),'%s',fig.UserData.status.Text);
page=v.fluidGrid;v.tabs.SelectedTab=page.tab;drawnow;
for entry={'structured','hybrid'}
    kind=entry{1};page.method.Value=kind;page.method.ValueChangedFcn(page.method,[]);d=page.report.(kind);
    assert(contains(page.outcome.Text,'0 m²') && contains(page.outcome.Text,sprintf('%d 个节点低于',d.stages(end).frames(end).below)));
    for level={'outer','joint','body'}
        page.zoom.Value=level{1};page.zoom.ValueChangedFcn(page.zoom,[]);drawnow;
        expected=d.section.pointsD(d.section.edges(:,1),1);
        assert(isequal(page.line.XData(1:3:end)',expected) && all(isnan(page.line.XData(3:3:end))));
    end
    proof.(kind)=struct('mesh',d.meshSHA256,'sectionEdges',size(d.section.edges,1),'frame',d.stages(end).frames(end).sha256);
end
% 模块二：在同一页重新生成混合体网；对照本次磁盘结果，确认未启动物理求解。
page.build.ButtonPushedFcn(page.build,[]);job=page.message.UserData;
assert(~isempty(job) && isfile(fullfile(job,'prepared.json')),'%s',page.message.Text);
p=jsondecode(fileread(fullfile(job,'prepared.json')));assert(strcmp(p.files.mesh_su2,proof.hybrid.mesh));
assert(~isfile(fullfile(job,'endpoints.json')) && isempty(fig.UserData.workflowState.rawMesh));
proof.rebuiltJob=job;proof.rebuiltMesh=p.files.mesh_su2;proof.transferFaces=p.geometry.transfer.targetFaces;
v.tabs.SelectedTab=v.geometryTab;drawnow;assert(height(v.geometry.stations)==7);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
