function check_workflow(out)
% 检查主窗口入口、独立工况的表格和三维页、来源身份，以及返回已有输入。
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
fig.Visible='off';before=fig.UserData.workflowState;
b=fig.UserData.controls.lewisCylinder;b.ButtonPushedFcn(b,[]);drawnow;
view=fig.UserData.lewisCylinderReview;
assert(~isempty(view),fig.UserData.status.Text);
assert(isequal(view.tabs.SelectedTab,view.comparisonTab) && isequaln(before,fig.UserData.workflowState));
tableUI=findobj(view.comparisonTab,'Type','uitable');
assert(isscalar(tableUI) && isequal(tableUI.Data(:,4),round(view.stations.figureNu,2)));
writetable(view.stations,fullfile(out,'stations.csv'));
view.tabs.SelectedTab=view.surface.tab;drawnow;
html=view.surface.html;d=html.UserData.payload;
assert(strcmp(d.defaultField,'nusselt') && ~d.researchIdentity.productionQualified);
f=d.fields(strcmp({d.fields.key},'nusselt'));idx=find(f.status==1,1);
assert(nnz(f.status)==220 && numel(f.values)==2160);
assert(all(d.fields(strcmp({d.fields.key},'pressure')).status==0));
assert(all(d.fields(strcmp({d.fields.key},'heatFlux')).status==0));
t=tic;while toc(t)<20 && isempty(html.UserData.rendered),drawnow;pause(.1);end
rendered=html.UserData.rendered;
assert(~isempty(rendered)&&rendered.faces==2160 && rendered.visiblePixelSamples>0 && rendered.webglError==0);
for id=[1 idx]
    html.Data=struct('kind','selection','face',id);html.DataChangedFcn(html,[]);
    assert(html.UserData.selectedFace==id);
end
b=fig.UserData.workflowButtons(1);b.ButtonPushedFcn(b,[]);drawnow;
assert(isempty(fig.UserData.lewisCylinderReview)&&~isgraphics(view.tabs)&&isequaln(before,fig.UserData.workflowState));
[mesh,stations,summary]=runLewisCylinderBenchmark();
assert(isequaln(stations,view.stations) && isequaln(summary,view.report));
% 圆柱适用域边界、负角对称、驻点、跨过端界及空工况量分别核对。
assert(isfinite(calcLewisCylinderNu(55,0))&&isnan(calcLewisCylinderNu(55,.01)));
assert(calcLewisCylinderNu(-50,1)==calcLewisCylinderNu(50,1));
assert(isnan(calcLewisCylinderNu(-60,0)));
page=exportAerothermalWeb(mesh,fullfile(out,'workbench'));assert(isfile(page));
receipt=struct('selectFaceEvents',[1 idx],'mouseDriven',false,'rendered',rendered, ...
    'snapshotPreserved',true,'stations',height(stations),'report',summary);
fid=fopen(fullfile(out,'workflow.json'),'w','n','UTF-8');assert(fid>=0);
g=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(receipt));
end
