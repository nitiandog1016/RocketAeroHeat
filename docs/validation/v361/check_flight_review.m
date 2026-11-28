function check_flight_review(out)
%CHECK_FLIGHT_REVIEW 实际操作同一侧栏，核对五个档案、当前边界和风洞归档。
%% 1. 冻结低速风洞尾迹工况。
app=rocketAeroHeatWorkbench(createCaseConfig('simpleLauncher'));
cleanup=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';
press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
selector=app.UserData.controls.launcherState;
selector.Value='wake';selector.ValueChangedFcn(selector,[]);
press(app.UserData.controls.confirm);
before=app.UserData.workflowState;
assert(before.currentStep==4 && isempty(before.result) && ~before.steps(5).available);
press(app.UserData.controls.openFlightRecords);
assert(~isfield(app.UserData.controls,'run') && ~isfield(app.UserData.controls,'readLauncher'));

%% 2. 逐档确认未返回、已有局部计算和给定压力的显示，浏览不产生求解任务。
mach=[.30 .85 1 1.2 2];
assert(isequal(app.UserData.controls.flightRecord.ItemsData,mach));
reference=jsondecode(fileread(fullfile(before.cfg.projectRoot,'docs','validation','v358','evidence','cases.json')));
review=cell(size(mach));
for i=1:numel(mach)
    field=app.UserData.controls.flightRecord;field.Value=mach(i);field.ValueChangedFcn(field,[]);
    c=app.UserData.controls; values=c.flightCoverage.Data;entry=reference.cases{i};note=strjoin(c.flightNote.Value,newline);
    assert(contains(c.flightCondition.Text,sprintf('%.2f m/s',entry.speedMps)));
    assert(contains(c.flightCondition.Text,'9216') && contains(note,'功率未知'));
    if i<=3
        assert(all(strcmp(values(:,2:3),'—'),'all') && contains(note,'未返回面场'));
        assert(contains(note,'球体对照')==(i==1));
    else
        assert(strcmp(values{3,2},num2str(entry.coverage.jointFaces)) && ...
            strcmp(values{3,3},sprintf('%.2f%%',100*entry.coverage.jointAreaM2/entry.coverage.areaM2)));
        assert(strcmp(values{4,2},num2str(entry.coverage.prescribedPressureFaces)));
        assert(contains(note,entry.heatMethod));
        for scopeLine=reshape(string(entry.pressureScope),1,[])
            assert(contains(note,scopeLine));
        end
    end
    nowState=app.UserData.workflowState;
    assert(isequaln(before.cfg,nowState.cfg) && isequaln(before.freestream,nowState.freestream));
    assert(before.revision==nowState.revision && before.solveId==nowState.solveId && ...
        isempty(nowState.result) && ~nowState.steps(5).available);
    review{i}=struct('mach',mach(i),'condition',c.flightCondition.Text,'visibleValues',{values},'note',note);
end

%% 3. 返回后仍可从同一风洞条件取得原有的完整结果，再切页时恢复日常入口。
press(app.UserData.controls.backToTunnel);
current=app.UserData.workflowState;
assert(isequaln(current.cfg,before.cfg) && ~isfield(app.UserData.controls,'flightRecord'));
press(app.UserData.controls.readLauncher);
current=app.UserData.workflowState;
assert(current.currentStep==5 && size(current.result.F,1)==9216 && current.steps(5).available);
assert(isequaln(current.freestream,before.freestream));
press(app.UserData.workflowButtons(4));press(app.UserData.controls.openFlightRecords);
press(app.UserData.workflowButtons(3));press(app.UserData.workflowButtons(4));
assert(isfield(app.UserData.controls,'readLauncher') && ~isfield(app.UserData.controls,'flightRecord'));
result=struct('cases',{review},'geometrySHA256',before.modelSource.sha256, ...
    'initialRevision',before.revision,'archivedWindTunnelFaces',size(current.result.F,1), ...
    'archivedWindTunnelElapsedS',current.freestream.elapsedPhysicalTimeS,'sameBoundary',true);
writelines(string(jsonencode(result,'PrettyPrint',true)),fullfile(out,'checked.json'),'Encoding','UTF-8');
end

function press(button)
button.ButtonPushedFcn(button,[]);
drawnow;
end
