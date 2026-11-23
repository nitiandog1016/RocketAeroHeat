function check_direction_ui(out,g)
% 模块一：从真实五步侧栏载入银球工况，确认三维页就绪及七站表的方向和单位。
fig=rocketAeroHeatWorkbench();done=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value=['supportDirection' upper(g(1)) g(2:end)];
c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'direction'),'%s',fig.UserData.status.Text);
html=v.web.html;t=tic;
while isempty(html.UserData.rendered) && toc(t)<60,drawnow;pause(.12);end
gl=html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==32480);
v.tabs.SelectedTab=v.direction.tab;drawnow;t=v.direction.table.Data;
assert(height(t)==7 && isnan(t.tau10(1)) && t.fine(1)=="方向退化" && t.fine(6)=="顺向" && t.coarse(7)=="反向");
assert(abs(t.tau10(6)-.02874429172147)<1e-12 && abs(t.u5(7)+.00910772672552)<1e-12);
% 模块二：逐档对照同面选择、来源说明、缺值和切向模，保留全表压力热流的台账。
v.tabs.SelectedTab=v.web.tab;e=v.web.evidenceView;d=e.getData();w=v.mesh.shortSupportDirection;
tag=['support' upper(g(1)) g(2:end)];
for suffix={'ShearPa','SpeedMps','ShearMagnitudePa','SpeedMagnitudeMps','Status'}
    col=suffix{1};key=[tag col];field=d.fields(strcmp({d.fields.key},key));values=w.([g col]);
    if strcmp(col,'Status'),values(values==0)=NaN;end
    assert(isequaln(field.values(:),values(:)) && isequal(field.status(:),double(isfinite(values))));
    expected=max(abs([w.(['fine' col]);w.(['coarse' col])]),[],'omitnan');
    assert(d.comparisonColorRanges.(key)(2)==expected);
    e.selectField(key);e.inspectFace(25055,[]);s=e.getState();assert(strcmp(s.fieldKey,key) && s.selectedFaceID==25055);
end
field=d.fields(strcmp({d.fields.key},[tag 'Status']));
assert(field.values(970)==1 && field.values(918)==-1 && isnan(field.values(164)) && contains(field.missingReason{164},'轴点'));
j=find(contains(field.missingReason,'棱边'),1);assert(isnan(field.values(j)));
e.inspectFace(j,[]);s=e.getState();assert(s.selectedFaceID==j);
key=[tag 'ShearPa'];e.selectField(key);x=.00137;e.setSectionX(x);timer=tic;
while toc(timer)<45
    drawnow;cut=e.getSection();
    if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,key)&&abs(cut.x-x)<1e-13,break;end
    pause(.10);
end
assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,key)&&abs(cut.x-x)<1e-13);
e.exportSectionTo(fullfile(out,'shear_section.csv'));
e.selectField([tag 'Status']);e.inspectFace(918,[]);s=e.getState();assert(s.selectedFaceID==918);
page=exportAerothermalWeb(v.mesh,fullfile(out,'workbench'));
proof=struct('faces',gl.faces,'pixels',gl.visiblePixelSamples,'directionFields',5,'stations',height(t), ...
    'sectionX',x,'faceIDs',[164 918 970 25055],'page',page,'grid',g,'frame',v.report.frameSHA256);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
