function check_time_review(out)
%CHECK_TIME_REVIEW 在实际工作台核对时间步身份、空间诊断、原面拾取和最终结果资格。
cfg=createCaseConfig('simpleLauncher');app=rocketAeroHeatWorkbench(cfg);
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));press(app.UserData.controls.confirm);
before=app.UserData.workflowState;press(app.UserData.controls.openFlightRecords);
note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
source=fullfile(cfg.projectRoot,'docs','validation','v366','research','status.json');
study=jsondecode(fileread(source));
assert(contains(note(),'50 µs 整步及首个 25 µs 半步均未收敛，同刻热流差未知'));
shares=sprintf('%.1f%%、%.1f%%',100*study.fullEnergyRimShare,100*study.halfEnergyRimShare);
assert(contains(note(),shares));
%% 点选尾段、圆柱、平底；保留既有面值，残差描述仅在对应范围出现。
html=app.UserData.geometryView.html;
for face=[1 2017 9216]
    waitForFrame(html);
    html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',face-1);
    html.DataChangedFcn(html,[]);drawnow;waitForFrame(html);
    assert(html.UserData.rendered.inspectedFace==face-1 && contains(note(),'同刻热流差未知'));
    assert(contains(note(),sprintf('原面 %d · 4 / 2 µm 末帧热流',face)));
    assert(contains(note(),shares)==(face~=2017));
end
%% 其他马赫数不复用该步长资格；实际计算和既有边界不受档案浏览影响。
c=app.UserData.controls.flightRecord;
for mach=[.85 1 1.2 2 .3]
    c.Value=mach;c.ValueChangedFcn(c,[]);drawnow;
    assert(contains(note(),'同刻热流差未知')==(mach==.3));
end
s=app.UserData.workflowState;
assert(s.solveId==before.solveId && isequaln(s.freestream,before.freestream) && isempty(s.result) && ~s.steps(5).available);
proof=struct('picked',[1 2017 9216],'checkedMach',[.3 .85 1 1.2 2],'pairedTimeDifferenceAvailable',false, ...
    'sourceSHA256',cm011FileSHA256(source),'legacyProbeValuesKept',true,'finalPageAvailable',false);
writelines(string(jsonencode(proof)),fullfile(out,'time.json'),'Encoding','UTF-8');
end

function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end

function waitForFrame(html)
timer=tic;
while toc(timer)<15
    drawnow;v=html.UserData;
    if ~isempty(v.error),error('%s',v.error.message);end
    if ~isempty(v.rendered) && v.rendered.revision==v.revision
        assert(v.rendered.webglError==0 && v.rendered.visiblePixelSamples>0);return
    end
    pause(.05);
end
error('未收到当前帧的绘制回执。');
end
