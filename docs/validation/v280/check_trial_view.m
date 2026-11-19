function check_trial_view(out)
% 模块一：从工作台的公开圆柱入口计算同一工况，核对试算切换和边界变化。
fig=rocketAeroHeatWorkbench();fig.Visible='off';
clean=onCleanup(@()delete(fig)); %#ok<NASGU>
button=fig.UserData.controls.lewisCylinder;button.ButtonPushedFcn(button,[]);
v=fig.UserData.lewisCylinderReview;v.tabs.SelectedTab=v.caseTab;
v.calculate.ButtonPushedFcn(v.calculate,[]);drawnow;
active=v.caseHost.UserData;assert(~isempty(active),v.caseStatus.Text);
original=active.web.html.UserData.payload;
assert(strcmp(v.trialScheme.Enable,'on') && strcmp(v.trialScheme.Value,'fixed'));
v.trialButton.ButtonPushedFcn(v.trialButton,[]);drawnow;
baseline=active.web.html.UserData.payload;
assert(contains(field(baseline,'trialHeatFlux').scope,'120 次'));
v.trialScheme.Value='adaptive';v.trialScheme.ValueChangedFcn(v.trialScheme,[]);drawnow;
candidate=active.web.html.UserData.payload;
assert(isequaln(candidate.wallCoverage,original.wallCoverage) && ...
    candidate.wallCoverage.jointComputedFaceCount==220 && isnan(candidate.wallCoverage.wholeSurfaceHeatPowerW));
assert(isequaln(field(original,'heatFlux'),field(candidate,'heatFlux')));
assert(isequaln(field(original,'pressure'),field(candidate,'pressure')));

% 模块二：按两个归档的源帧和 2160 个原面值，分别核对实际工作台与独立导出。
base=fileparts(mfilename('fullpath'));
for k=1:2
    if k==1,page=baseline;folder=fullfile(base,'..','v279','evidence');
    else,page=candidate;folder=fullfile(base,'evidence');end
    r=jsondecode(fileread(fullfile(folder,'trial.json')));
    rows=readtable(fullfile(folder,'trial_faces.csv'));
    p=field(page,'trialPressure');q=field(page,'trialHeatFlux');
    assert(contains(p.scope,r.sourceFrameSHA256)&&contains(q.scope,r.sourceFrameSHA256));
    assert(contains(q.scope,sprintf('%d 次',r.innerIterations)) && strcmp(q.computedStateLabel,'未收敛试算'));
    assert(max(abs(p.values*1000-rows.trialPressurePa))<1e-8 && max(abs(q.values*1000-rows.trialHeatFluxWm2))<1e-9);
end
mesh=attachLewisTunnelTrial(active.mesh,'adaptive');
assert(isequaln(buildAerothermalViewerData(mesh),candidate));
exportAerothermalWeb(mesh,fullfile(out,'workbench'));
v.trialScheme.Value='fixed';v.trialScheme.ValueChangedFcn(v.trialScheme,[]);drawnow;
assert(isequaln(active.web.html.UserData.payload,baseline));
v.inputs(3).Value=283.15;v.inputs(3).ValueChangedFcn(v.inputs(3),[]);drawnow;
assert(isempty(v.caseHost.UserData)&&strcmp(v.trialScheme.Enable,'off')&&~isgraphics(active.web.tab));
v.calculate.ButtonPushedFcn(v.calculate,[]);drawnow;
cold=v.caseHost.UserData;assert(~isempty(cold)&&strcmp(v.trialScheme.Enable,'off'));
for s={'fixed','adaptive'}
    try
        attachLewisTunnelTrial(cold.mesh,s{1});error('冷壁工况误配归档');
    catch ex
        assert(strcmp(ex.identifier,'RAH:Lewis:TrialCondition'));
    end
end
fid=fopen(fullfile(out,'workflow.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(struct('checkedCases',{{'fixed','adaptive','fixed','cold'}}, ...
    'checkedFacesPerRun',2160,'jointEngineeringFaces',candidate.wallCoverage.jointComputedFaceCount, ...
    'partialPowerW',candidate.wallCoverage.calculatedHeatPowerW,'wholePowerW',candidate.wallCoverage.wholeSurfaceHeatPowerW, ...
    'sourceFrameSHA256',mesh.lewisTrial.sourceFrameSHA256,'mouseDriven',false)));
end

function f=field(page,key)
% 按变量唯一标识读取，兼容列表顺序的变化。
f=page.fields(strcmp({page.fields.key},key));assert(isscalar(f));
end
