function check_support_ui(out)
% 模块一：经普通侧栏生成带支管几何，等候三维像素及真实接头线段的回执。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supports';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'scene'),'%s',fig.UserData.status.Text);
g=v.geometry;n=size(g.F,1);count=size(g.joinPositive,1)+size(g.joinNegative,1);
base=rendered(v.scene.html);assert(base.faces==n && base.inspectedFace==-1 && base.sectionSegments==count);
proof=struct('faces',n,'joinSegments',count,'pixels',base.visiblePixelSamples,'stationFaces',zeros(7,1));
for k=1:7
    v.choice.Value=k;v.choice.ValueChangedFcn(v.choice,[]);r=rendered(v.scene.html);
    assert(r.inspectedFace==g.stations.faceID(k)-1 && contains(v.info.Text,sprintf('原球面号 %d',g.stations.parentID(k))));
    assert(~isfield(v.scene.html.UserData.payload,'faceField'));
    proof.stationFaces(k)=r.inspectedFace+1;
end
% 模块二：检查实际导出、完全遮盖的原面及返回原稳态测站的行为。
v.tabs.SelectedTab=v.areaTab;drawnow;proof.affectedRows=size(v.ledger.Data,1);
assert(proof.affectedRows==nnz(g.parents.coveredM2>0));
v.exportTo(fullfile(out,'supported.zip'));
c.shortSphereMode.Value='archive';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
old=fig.UserData.shortSphereReview;assert(strcmp(old.report.action,'archive') && height(old.stations)==7 && ~isfield(old,'geometry'));
assert(isempty(fig.UserData.workflowState.rawMesh));proof.steadyReturned=true;
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function result=rendered(html)
timer=tic;
while isempty(html.UserData.rendered) && isempty(html.UserData.error) && toc(timer)<30,drawnow;pause(.1);end
result=html.UserData.rendered;
assert(isempty(html.UserData.error) && ~isempty(result) && result.visiblePixelSamples>0 && result.webglError==0);
end
