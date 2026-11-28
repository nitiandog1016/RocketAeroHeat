function check_pressure_review(out)
%CHECK_PRESSURE_REVIEW 在工作台的正常入口核对档案数值、原面拾取和结果页资格。
%% 1. 用同形风洞模型进入第 4 步，再打开飞行研究档案。
cfg=createCaseConfig('simpleLauncher');
study=jsondecode(fileread(fullfile(cfg.projectRoot,'docs','validation','v364','research','surface.json')));
app=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
press(app.UserData.controls.confirm);before=app.UserData.workflowState;
press(app.UserData.controls.openFlightRecords);html=app.UserData.geometryView.html;
note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
assert(contains(note(),'同形研究外流') && contains(note(),'伯努利绝对静压'));
%% 2. 按真实拾取事件逐区核对 14D / 22D 的两个面值；切换工况时排除旧场。
faces=[1 1537 2017 6912 6913 8833 9216];
for face=faces
    waitRender(html);
    html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',face-1);
    html.DataChangedFcn(html,[]);drawnow;waitRender(html);
    assert(html.UserData.rendered.inspectedFace==face-1);
    k=study.faceGroup(face);
    p=study.freePressurePa+study.dynamicPressurePa*[study.groupCp14D(k) study.groupCp22D(k)];
    expected=sprintf('原面 %d · 研究静压（绝对）14D %.1f Pa，22D %.1f Pa。',face,p(1),p(2));
    assert(contains(note(),expected) && contains(note(),'没有壁面换热方程'));
end
c=app.UserData.controls.flightRecord;
for mach=[.85 1 1.2 2 .3]
    c.Value=mach;c.ValueChangedFcn(c,[]);drawnow;
    assert(contains(note(),'同形研究外流')==(mach==.3));
end
s=app.UserData.workflowState;
assert(isempty(s.result) && ~s.steps(5).available && s.solveId==before.solveId && isequaln(s.freestream,before.freestream));
proof=struct('pickedFaces',faces,'domainsD',[14 22],'checkedMach',[.3 .85 1 1.2 2], ...
    'methodAndHeatScopeShown',true,'finalResultAvailable',s.steps(5).available);
writelines(string(jsonencode(proof)),fullfile(out,'picked.json'),'Encoding','UTF-8');
end

function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end

function waitRender(html)
% 等待本次绘制的原面回执，避免把上一帧当作拾取结果。
timer=tic;
while toc(timer)<15
    drawnow;v=html.UserData;
    if ~isempty(v.error),error('%s',v.error.message);end
    if ~isempty(v.rendered) && v.rendered.revision==v.revision
        assert(v.rendered.webglError==0 && v.rendered.visiblePixelSamples>0);return
    end
    pause(.05);
end
error('当前几何绘制未返回。');
end
