function check_assembly_review(out)
%CHECK_ASSEMBLY_REVIEW 核对局部能量账只在对应工况和尾部范围显示。
cfg=createCaseConfig('simpleLauncher');app=rocketAeroHeatWorkbench(cfg);
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));press(app.UserData.controls.confirm);
before=app.UserData.workflowState;press(app.UserData.controls.openFlightRecords);
note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
src=fullfile(cfg.projectRoot,'docs','validation','v367','evidence','status.json');
r=jsondecode(fileread(src));expected=sprintf('流体节点 %d（X=%.2f mm）的中心对流项 %.3g',r.peakPoint,1000*r.modelXM,r.peakCentral);
assert(~contains(note(),expected));html=app.UserData.geometryView.html;
%% 点选平底、圆柱和尾段，区分流体控制体与所点的原三角面。
for face=[9216 2017 1]
    waitFrame(html);
    html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',face-1);
    html.DataChangedFcn(html,[]);drawnow;waitFrame(html);
    msg=note();assert(html.UserData.rendered.inspectedFace==face-1);
    assert(contains(msg,expected)==(face~=2017));
    assert(contains(msg,sprintf('净能量残差 %.3g',r.peakPhysicalResidual))==(face~=2017));
    assert(contains(msg,'两种步长仍未收敛')==(face~=2017));
    assert(contains(msg,sprintf('原面 %d · 4 / 2 µm 末帧热流',face)));
end
%% 检查研究资料不改变工况边界、计算状态和结果资格。
c=app.UserData.controls.flightRecord;c.Value=.85;c.ValueChangedFcn(c,[]);drawnow;
assert(~contains(note(),'局部组装'));
s=app.UserData.workflowState;
assert(isequaln(s.freestream,before.freestream) && s.solveId==before.solveId && ~s.steps(5).available && isempty(s.result));
proof=struct('pickedFaces',[9216 2017 1],'sourceSHA256',cm011FileSHA256(src),'fluidPoint',r.peakPoint, ...
    'caseSwitchChecked',true,'finalResultAvailable',false);
writelines(string(jsonencode(proof)),fullfile(out,'assembly.json'),'Encoding','UTF-8');
end

function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end

function waitFrame(html)
t=tic;
while toc(t)<15
    drawnow;v=html.UserData;
    if ~isempty(v.error),error('%s',v.error.message);end
    if ~isempty(v.rendered) && v.rendered.revision==v.revision
        assert(v.rendered.webglError==0 && v.rendered.visiblePixelSamples>0);return
    end
    pause(.05);
end
error('未收到当前原面视图回执。');
end
