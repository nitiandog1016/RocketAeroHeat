function check_short_step(out)
%CHECK_SHORT_STEP 从正常工作台核对六区原面的研究压力、热流、温度及资格。
cfg=createCaseConfig('simpleLauncher');
app=rocketAeroHeatWorkbench(cfg);
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';
press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
press(app.UserData.controls.confirm);press(app.UserData.controls.openFlightRecords);
folder=fullfile(cfg.projectRoot,'docs','validation');
path=fullfile(folder,'v369','evidence','surface.json');
t=jsondecode(fileread(path));
groups=jsondecode(fileread(fullfile(folder,'v365','research','surface.json')));
html=app.UserData.geometryView.html;
note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
stamp=sprintf('自数值初场起 %.2f µs，步内收敛',t.elapsedFromNumericalSeedUs);
assert(contains(note(),stamp) && contains(note(),sprintf('全表面 9216 面 · 入壁功率 %.3f W',t.regions.all.researchIntoWallW)));
% 每个区域使用不同原面，逐项核对用户实际看到的数值、单位和当前三维选区。
faces=[26 1602 2017 6951 8900 9216];
keys={'afterbody','contraction','cylinder','cone','sphericalNose','base'};
for j=1:numel(faces)
    face=faces(j);waitFrame(html);
    html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',face-1);
    html.DataChangedFcn(html,[]);drawnow;waitFrame(html);
    assert(html.UserData.rendered.inspectedFace==face-1);
    v=t.groupValues(groups.faceGroup(face),:);r=t.regions.(keys{j});
    msg=note();
    assert(contains(msg,stamp) && contains(msg,sprintf('入壁功率 %.3f W',r.researchIntoWallW)));
    assert(contains(msg,sprintf('原面 %d · 绝对壁压 %.1f Pa，入壁热流 %.1f W/m²',face,v(1),v(2))));
    assert(contains(msg,sprintf('给定壁温 %.2f K，壁节点流体 %.2f K，首层流体 %.2f K',t.wallGivenK,v(3),v(4))));
end
% 切换工况应撤下此帧；返回后仍保留其研究身份。
control=app.UserData.controls.flightRecord;
for mach=[.85 1 1.2 2 .3]
    control.Value=mach;control.ValueChangedFcn(control,[]);drawnow;
    assert(contains(note(),stamp)==(mach==.3));
end
s=app.UserData.workflowState;
assert(isempty(s.result) && ~s.steps(5).available);
proof=struct('sourceSHA256',cm011FileSHA256(path),'checkedFaces',faces,'checkedMach',[.3 .85 1 1.2 2], ...
    'elapsedFromNumericalSeedUs',t.elapsedFromNumericalSeedUs,'finalResultAvailable',false);
writelines(string(jsonencode(proof)),fullfile(out,'shortStep.json'),'Encoding','UTF-8');
end

function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end

function waitFrame(html)
% 等当前版本的三维绘制完成，再发出下一条拾取消息。
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
