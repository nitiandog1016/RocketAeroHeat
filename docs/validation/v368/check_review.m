function check_review(out)
%CHECK_REVIEW 核对正常侧栏里的重建对照与对应工况、区域和结果资格。
cfg=createCaseConfig('simpleLauncher');
app=rocketAeroHeatWorkbench(cfg);
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';
press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
press(app.UserData.controls.confirm);press(app.UserData.controls.openFlightRecords);
html=app.UserData.geometryView.html;
note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
src=fullfile(cfg.projectRoot,'docs','validation','v368','evidence','status.json');
g=jsondecode(fileread(src));
label=sprintf('邻近四列的能量残差范数为原方法的 %.2f 倍', ...
    g.methods.gg.fourColumnRhsNorm/g.methods.weighted.fourColumnRhsNorm);
assert(~contains(note(),label));
% 三个既有原面按真实的拾取消息进入；圆柱不应显示尾部方程读数。
for face=[9216 2017 1]
    waitFrame(html);
    html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',face-1);
    html.DataChangedFcn(html,[]);drawnow;waitFrame(html);
    assert(html.UserData.rendered.inspectedFace==face-1);
    assert(contains(note(),label)==(face~=2017));
    assert(contains(note(),sprintf('原面 %d · 4 / 2 µm 末帧热流',face)));
end
% 切换飞行档案后读数随工况撤销，全面积验收继续控制结果页。
c=app.UserData.controls.flightRecord;c.Value=.85;c.ValueChangedFcn(c,[]);drawnow;
s=app.UserData.workflowState;
assert(~contains(note(),label) && ~s.steps(5).available && isempty(s.result));
proof=struct('sourceSHA256',cm011FileSHA256(src),'pickedFaces',[9216 2017 1], ...
    'caseSwitchChecked',true,'finalResultAvailable',false);
writelines(string(jsonencode(proof)),fullfile(out,'gradient.json'),'Encoding','UTF-8');
end

function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end

function waitFrame(html)
% 在异步三维视图确认本次绘制后，再发送下一次面元拾取。
t=tic;
while toc(t)<15
    drawnow;s=html.UserData;
    if ~isempty(s.error),error('%s',s.error.message);end
    if ~isempty(s.rendered) && s.rendered.revision==s.revision
        assert(s.rendered.webglError==0 && s.rendered.visiblePixelSamples>0);
        return
    end
    pause(.05);
end
error('未收到本次三维绘制的回执。');
end
