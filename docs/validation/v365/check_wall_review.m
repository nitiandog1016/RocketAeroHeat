function check_wall_review(out)
%CHECK_WALL_REVIEW 从普通五步侧栏核对跨区面元的两档壁面值和档案切换。
cfg=createCaseConfig('simpleLauncher');
src=fullfile(cfg.projectRoot,'docs','validation','v365','research','surface.json');
wall=jsondecode(fileread(src));app=rocketAeroHeatWorkbench(cfg);
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));press(app.UserData.controls.confirm);
before=app.UserData.workflowState;press(app.UserData.controls.openFlightRecords);
html=app.UserData.geometryView.html;note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
assert(contains(note(),'等温壁层流试算') && contains(note(),'转捩待核实'));
%% 逐面核对同一份压力、热流、壁温和顺流剪切的字段与单位。
faces=[1 1537 2017 4500 6912 6913 8833 9216];
for face=faces
    waitForFrame(html);
    html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',face-1);
    html.DataChangedFcn(html,[]);drawnow;waitForFrame(html);
    assert(html.UserData.rendered.inspectedFace==face-1);
    k=wall.faceGroup(face);a=wall.normal4(k,:);b=wall.normal2(k,:);
    assert(contains(note(),sprintf('原面 %d · 4 / 2 µm 末帧热流 %.1f / %.1f W/m²',face,a(2),b(2))));
    assert(contains(note(),sprintf('绝对壁压 %.1f / %.1f Pa',a(1),b(1))));
    assert(contains(note(),sprintf('沿来流壁剪切 %.4g Pa',b(3))) && contains(note(),sprintf('首层沿来流速度 %.4g m/s',b(4))));
    assert(contains(note(),sprintf('给定壁温 300.00 K，壁节点流体 %.3f K，首层流体 %.3f K',b(5),b(6))));
end
%% 切换工况后不得残留这批热壁值，研究档案始终不占用最终结果页。
c=app.UserData.controls.flightRecord;
for mach=[.85 1 1.2 2 .3]
    c.Value=mach;c.ValueChangedFcn(c,[]);drawnow;
    assert(contains(note(),'等温壁层流试算')==(mach==.3));
end
s=app.UserData.workflowState;
assert(isempty(s.result) && ~s.steps(5).available && s.solveId==before.solveId && isequaln(s.freestream,before.freestream));
proof=struct('faces',faces,'sourceSHA256',cm011FileSHA256(src),'fieldsReadBack',{{'p','q','shear','U1','TwGiven','TwallFluid','T1'}}, ...
    'checkedMach',[.3 .85 1 1.2 2],'finalPageAvailable',false);
writelines(string(jsonencode(proof)),fullfile(out,'wall.json'),'Encoding','UTF-8');
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
