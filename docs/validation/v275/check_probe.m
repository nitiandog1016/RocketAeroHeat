function check_probe(base)
%CHECK_PROBE 回读实际点击计算所得快照，检查结果宿主接收的选面编号与字段。
run=strtrim(fileread(fullfile(base,'gui01','run_folder.txt')));
saved=load(fullfile(run,'result.mat'),'current');
fig=uifigure('Visible','off');cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
tabs=uitabgroup(fig);ui=createAerothermalWebTab(tabs,saved.current);
started=tic;
while toc(started)<12 && isempty(ui.html.UserData.rendered),drawnow;pause(.1);end
state=ui.html.UserData;d=state.payload;
p=d.fields(strcmp({d.fields.key},'pressureGauge'));q=d.fields(strcmp({d.fields.key},'heatFlux'));
ids=[1 7000 size(d.faces,1)];
for id=ids
    ui.html.Data=struct('kind','selection','face',id);
    ui.html.DataChangedFcn(ui.html,[]);picked=ui.html.UserData.selectedFace;
    assert(picked==id && isequal(d.faces(picked,:)+1,saved.current.F(id,:)));
    assert(p.values(picked)==saved.current.hybridWall.pressureGaugePa(id) && ...
        abs(1000*q.values(picked)-saved.current.hybridWall.heatFluxWm2(id))<1e-14);
end
ui.html.Data=struct('kind','selection','face',ids(end)+1);ui.html.DataChangedFcn(ui.html,[]);
assert(ui.html.UserData.selectedFace==ids(end));
receipt=struct('faceIDs',ids,'matched',true,'outOfRangeRejected',true,'hostReady',state.ready, ...
    'renderedReceipt',state.rendered,'mouseDriven',false);
fid=fopen(fullfile(base,'probe.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(receipt));
end
