function check_workflow(out)
% 主窗口显式计算同一圆柱的给定工况，验证边界改动撤回旧场及原实测页身份。
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
fig.Visible='off';before=fig.UserData.workflowState;
button=fig.UserData.controls.lewisCylinder;button.ButtonPushedFcn(button,[]);drawnow;
v=fig.UserData.lewisCylinderReview;assert(~isempty(v),fig.UserData.status.Text);
archive=v.surface.html.UserData.payload;
assert(all(archive.fields(strcmp({archive.fields.key},'pressure')).status==0));
v.tabs.SelectedTab=v.caseTab;drawnow;
press(v.calculate);active=v.caseHost.UserData;
assert(~isempty(active),v.caseStatus.Text);
d=active.web.html.UserData.payload;coverage=active.report.wallCoverage;
assert(coverage.jointComputedFaceCount==220 && coverage.calculatedHeatPowerW<0 && ...
    isnan(coverage.wholeSurfaceHeatPowerW) && strcmp(d.defaultField,'heatFlux'));
key={d.fields.key};pf=d.fields(strcmp(key,'pressure'));qf=d.fields(strcmp(key,'heatFlux'));
idx=find(qf.status==1,1);assert(isequal(pf.status,qf.status));
assert(all(d.fields(strcmp(key,'wallInput')).status==2) && all(d.fields(strcmp(key,'wallSolved')).status==0));
assert(d.resultCondition.wallUniformValueK==303.15 && d.resultCondition.freestreamPressurePa==101325);
t=tic;while toc(t)<20 && isempty(active.web.html.UserData.rendered),drawnow;pause(.1);end
rendered=active.web.html.UserData.rendered;
assert(~isempty(rendered)&&rendered.faces==2160&&rendered.visiblePixelSamples>0&&rendered.webglError==0);
for id=[1 idx]
    active.web.html.Data=struct('kind','selection','face',id);
    active.web.html.DataChangedFcn(active.web.html,[]);assert(active.web.html.UserData.selectedFace==id);
end
assert(strcmp(v.trialButton.Enable,'on'));press(v.trialButton);
trialPage=active.web.html.UserData.payload;keys={trialPage.fields.key};
candidate=trialPage.fields(strcmp(keys,'trialHeatFlux'));
assert(isscalar(candidate)&&all(isfinite(candidate.values))&&strcmp(candidate.computedStateLabel,'未收敛试算'),v.caseStatus.Text);
assert(isequaln(trialPage.wallCoverage,d.wallCoverage)&&isequaln(trialPage.fields(strcmp(keys,'heatFlux')),qf));
change(v.inputs(3),283.15);assert(isempty(v.caseHost.UserData)&&~isgraphics(active.web.tab));
press(v.calculate);cold=v.caseHost.UserData;assert(cold.report.wallCoverage.calculatedHeatPowerW>0&&strcmp(v.trialButton.Enable,'off'));
try
    attachLewisTunnelTrial(cold.mesh);error('归档误配到冷壁工况');
catch ex
    assert(strcmp(ex.identifier,'RAH:Lewis:TrialCondition'));
end
change(v.inputs(1),1000);press(v.calculate);
assert(isempty(v.caseHost.UserData)&&contains(v.caseStatus.Text,'46 m/s'));
change(v.inputs(1),120000);press(v.calculate);assert(isempty(v.caseHost.UserData)&&contains(v.caseStatus.Text,'90–110 kPa'));
change(v.inputs(2),120);press(v.calculate);assert(isempty(v.caseHost.UserData)&&contains(v.caseStatus.Text,'280–310 K'));
assert(isequaln(archive,v.surface.html.UserData.payload)&&isequaln(before,fig.UserData.workflowState));
press(fig.UserData.workflowButtons(1));assert(~isgraphics(v.tabs));
% 独立导出冻结的默认输入。旧无量纲公开测站数据脚本应可逐字节重建。
[mesh,stations,report]=runLewisCylinderCondition(101325,293.15,303.15);
mesh=attachLewisTunnelTrial(mesh);
page=exportAerothermalWeb(mesh,fullfile(out,'workbench'));assert(isfile(page));
writetable(stations,fullfile(out,'stations.csv'));
assert(isequaln(report.wallCoverage,coverage));
legacy=runLewisCylinderBenchmark();original=buildAerothermalViewerData(legacy);
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
old=fileread(fullfile(root,'outputs','v278_validation','gui01','workbench','case-data.js'));
assert(strcmp(strtrim(old),['window.ROCKET_CASE=' jsonencode(original,'ConvertInfAndNaN',true) ';']));
tw=fzero(@zeroDelta,[293.15 295]);[~,~,balanced]=runLewisCylinderCondition(101325,293.15,tw);
assert(abs(balanced.wallCoverage.calculatedHeatPowerW)<1e-10 && balanced.wallCoverage.jointComputedFaceCount==220);
receipt=struct('rendered',rendered,'faceEvents',[1 idx],'mouseDriven',false,'condition',mesh.lewisCondition, ...
    'hot',coverage,'cold',cold.report.wallCoverage,'zeroLimitPowerW',balanced.wallCoverage.calculatedHeatPowerW, ...
    'zeroLimitWallK',tw,'staleViewRemoved',true,'invalidFlowRejected',true,'oldBenchmarkUnchanged',true);
fid=fopen(fullfile(out,'workflow.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(receipt));
end
function v=zeroDelta(tw)
c=lewisCylinderCondition(101325,293.15,tw);v=c.adiabaticWallK-tw;
end
function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end
function change(field,value)
field.Value=value;field.ValueChangedFcn(field,[]);drawnow;
end
