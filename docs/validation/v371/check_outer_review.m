function check_outer_review(out)
%CHECK_OUTER_REVIEW 在五步侧栏核对两网同刻云图、同色标、原面探针及研究资格。
cfg=createCaseConfig('simpleLauncher');
app=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';
press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
press(app.UserData.controls.confirm);press(app.UserData.controls.openFlightRecords);
folder=fullfile(cfg.projectRoot,'docs','validation');
path=fullfile(folder,'v371','evidence','comparison.json');
d=jsondecode(fileread(path));t=jsondecode(fileread(fullfile(folder,'v369','evidence','surface.json')));
wall=jsondecode(fileread(fullfile(folder,'v365','research','surface.json')));
html=app.UserData.geometryView.html;field=app.UserData.controls.flightField;
note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
waitFrame(html);assert(strcmp(field.Value,'geometry') && isempty(html.UserData.rendered.fieldRange));
modes={'fineHeat','coarseHeat','finePressure','coarsePressure','fineHeat','finePressure'};
faces=[1 1602 2017 6951 8900 9216];
for j=1:numel(modes)
    % 先切换场量，再从真实绘制版本拾取；两个网格共享同一物理量的全表面色标。
    mode=modes{j};field.Value=mode;field.ValueChangedFcn(field,[]);drawnow;waitFrame(html);
    col=1;if endsWith(mode,'Heat'),col=2;end
    v=t.groupValues;if startsWith(mode,'fine'),v=d.fineGroupValues;end
    limits=[min([d.fineGroupValues(:,col);t.groupValues(:,col)]) max([d.fineGroupValues(:,col);t.groupValues(:,col)])];
    payload=html.UserData.payload;actual=payload.faceField;
    assert(strcmp(actual.context,'single-step-research') && ~actual.productionQualified && ~payload.resetCamera);
    assert(isequal(actual.values,v(wall.faceGroup,col)) && isequal(actual.range,limits));
    assert(html.UserData.rendered.faces==9216 && max(abs(html.UserData.rendered.fieldRange(:)-limits(:)))<1e-9);
    face=faces(j);html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',face-1);
    html.DataChangedFcn(html,[]);drawnow;waitFrame(html);
    row=v(wall.faceGroup(face),:);
    assert(contains(note(),sprintf('原面 %d · 绝对壁压 %.5f Pa，入壁热流 %.5f W/m²',face,row(1),row(2))));
    assert(html.UserData.rendered.inspectedFace==face-1 && contains(note(),'外层加密对照'));
end
% 研究云图不跨工况传递；返回风洞结果和发布资格仍由正常流程控制。
mach=app.UserData.controls.flightRecord;
for m=[.85 1 1.2 2 .3]
    mach.Value=m;mach.ValueChangedFcn(mach,[]);drawnow;waitFrame(html);
    assert(strcmp(field.Value,'geometry') && isempty(html.UserData.rendered.fieldRange));
    assert(strcmp(field.Enable,'on')==(m==.3) && contains(note(),'外层加密对照')==(m==.3));
end
s=app.UserData.workflowState;assert(isempty(s.result) && ~s.steps(5).available);
proof=struct('sourceSHA256',cm011FileSHA256(path),'modes',{modes},'checkedFaces',faces,'matchedValuesPerMode',9216, ...
    'commonRangesChecked',true,'caseSwitchChecked',true,'finalResultAvailable',false);
writelines(string(jsonencode(proof)),fullfile(out,'outer.json'),'Encoding','UTF-8');
end

function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end

function waitFrame(html)
% 绘制回执须属于当前数据版本，且网页确实画出了可见像素。
t=tic;
while toc(t)<15
    drawnow;s=html.UserData;
    if ~isempty(s.error),error('%s',s.error.message);end
    if ~isempty(s.rendered) && s.rendered.revision==s.revision
        assert(s.rendered.webglError==0 && s.rendered.visiblePixelSamples>0);return
    end
    pause(.05);
end
error('未收到本次三维绘制的回执。');
end
