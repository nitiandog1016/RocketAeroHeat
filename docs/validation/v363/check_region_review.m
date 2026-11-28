function check_region_review(out)
%CHECK_REGION_REVIEW 用原面表独立检查侧栏、拾取事件及实际 WebGL 绘制回执。
%% 1. 从普通风洞工况入口打开档案，保留已确认状态用于返回核对。
cfg=createCaseConfig('simpleLauncher');app=rocketAeroHeatWorkbench(cfg);
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
setChoice(app.UserData.controls.launcherState,'wake');press(app.UserData.controls.confirm);
before=app.UserData.workflowState;press(app.UserData.controls.openFlightRecords);
html=app.UserData.geometryView.html;
geometry=jsondecode(fileread(fullfile(cfg.projectRoot,'cases','simple_launcher','geometry.json')));
archive=jsondecode(fileread(fullfile(cfg.projectRoot,'docs','validation','v358','evidence','cases.json')));
wind=jsondecode(fileread(fullfile(cfg.projectRoot,'docs','validation','v363','wind_domain.json')));
frame=waitFrame(html);
assert(frame.selectedFaceCount==0 && app.UserData.controls.flightRecord.Value==.30);

%% 2. 从列表选择六区，再通过既有网页拾取消息进入不同区域；核对两端编号。
reports=cell(numel(geometry.regions),1);
for i=1:numel(geometry.regions)
    r=geometry.regions(i);
    setChoice(app.UserData.controls.flightRegion,r.name);
    frame=waitFrame(html);
    faceIDs=cell2mat(html.UserData.payload.selectedFaceIDs);
    assert(isequal(faceIDs(:)',r.firstFace-1:r.lastFace-1));
    assert(frame.selectedFaceCount==numel(faceIDs) && frame.inspectedFace==-1);
    c=app.UserData.controls;
    assert(contains(c.flightCondition.Text,sprintf('%.6f m²',r.areaM2)));
    assert(strcmp(c.flightCoverage.ColumnName{3},'区内面积比'));
    assert(all(strcmp(c.flightCoverage.Data(:,2:3),'—'),'all'));
    note=strjoin(c.flightNote.Value,newline);domain=wind.regions.(r.name);
    assert(contains(note,'本区流动状态未判定') && contains(note,'同形稳态风洞（0.05 m/s）'));
    assert(contains(note,sprintf('风洞动压的 %.4f%%',domain.pressureRmsPercentDynamic)) && ...
        contains(note,sprintf('功率差 %.4f%%',domain.powerDifferencePercent)));
    next=geometry.regions(mod(i,numel(geometry.regions))+1);
    face=next.lastFace;sendProbe(html,html.UserData.revision,face-1);
    hit=waitFrame(html);c=app.UserData.controls;
    assert(strcmp(c.flightRegion.Value,next.name) && hit.inspectedFace==face-1);
    assert(hit.selectedFaceCount==next.lastFace-next.firstFace+1);
    assert(contains(strjoin(c.flightNote.Value,newline),sprintf('选中原面 %d',face)));
    reports{i}=struct('region',r.name,'facesDrawn',frame.selectedFaceCount, ...
        'pickedFace',face,'pickedRegion',next.name,'renderedProbeIndex',hit.inspectedFace, ...
        'visiblePixelSamples',hit.visiblePixelSamples);
end

%% 3. 同一圆柱段在五个档案中的面积分母和字段身份逐项取自各自回执。
setChoice(app.UserData.controls.flightRegion,'cylinder');
oldRevision=html.UserData.revision;
records=cell(numel(archive.cases),1);
for i=1:numel(archive.cases)
    rec=archive.cases{i};setChoice(app.UserData.controls.flightRecord,rec.mach);
    frame=waitFrame(html);c=app.UserData.controls;t=c.flightCoverage.Data;
    assert(frame.selectedFaceCount==4896 && frame.inspectedFace==-1);
    if rec.solverReturned
        r=rec.regions.cylinder;
        keys={'computedPressure','heat','joint','prescribedPressure'};
        for row=1:numel(keys)
            key=keys{row};
            assert(strcmp(t{row,2},sprintf('%d',r.([key 'Faces']))));
            assert(strcmp(t{row,3},sprintf('%.2f%%',100*r.([key 'AreaM2'])/r.areaM2)));
        end
    else
        assert(all(strcmp(t(:,2:3),'—'),'all'));
    end
    assert(contains(strjoin(c.flightNote.Value,newline),'本区流动状态')==(i==1));
    records{i}=struct('mach',rec.mach,'pressure',t{1,3},'heat',t{2,3},'joint',t{3,3},'prescribed',t{4,3});
end
rev=html.UserData.revision;sendProbe(html,oldRevision,0);sendProbe(html,rev,9216);
assert(html.UserData.revision==rev && strcmp(app.UserData.controls.flightRegion.Value,'cylinder'));
% 同一区内两张相邻原面的热流状态不同；不能用 99.02% 的区域汇总替代面状态。
for face=[2064 2065 6912 6913 8736 8737 9169]
    sendProbe(html,html.UserData.revision,face-1);frame=waitFrame(html);
    note=strjoin(app.UserData.controls.flightNote.Value,newline);
    if face<=6912 || face>=9169,expectedP='给定';else,expectedP='工程计算';end
    if face>=2065 && face<=8736,expectedQ='工程计算';else,expectedQ='缺项';end
    assert(frame.inspectedFace==face-1 && contains(note,['该面壁压 · ' expectedP]) && contains(note,['热流 · ' expectedQ]));
end
setChoice(app.UserData.controls.flightRegion,'all');frame=waitFrame(html);
assert(frame.selectedFaceCount==0 && strcmp(app.UserData.controls.flightCoverage.ColumnName{3},'全表面占比'));
assert(strcmp(app.UserData.controls.flightCoverage.Data{3,3},'14.75%'));

%% 4. 浏览后检查冻结边界、结果资格及原有读档路径。
s=app.UserData.workflowState;
assert(isequaln(s.cfg,before.cfg) && isequaln(s.freestream,before.freestream) && s.solveId==before.solveId);
assert(s.revision==before.revision && isempty(s.result) && ~s.steps(5).available);
press(app.UserData.controls.backToTunnel);frame=waitFrame(html);
payload=html.UserData.payload;
assert(frame.selectedFaceCount==0 && (~isfield(payload,'enableFaceProbe') || ~payload.enableFaceProbe));
press(app.UserData.controls.readLauncher);s=app.UserData.workflowState;
assert(s.currentStep==5 && s.steps(5).available && size(s.result.F,1)==9216 && isequaln(s.freestream,before.freestream));
proof=struct('regionChecks',{reports},'cylinderChecks',{records},'archivedTunnelFaces',size(s.result.F,1), ...
    'sameWindBoundary',true,'rejectedStaleAndOutOfRange',true,'mach2IndividualFaceChecks',[2064 2065 6912 6913 8736 8737 9169]);
writelines(string(jsonencode(proof,'PrettyPrint',true)),fullfile(out,'checked.json'),'Encoding','UTF-8');
end

function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end

function setChoice(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end

function sendProbe(html,revision,zeroIndex)
% 发送与网页命中三角形相同的零基消息，验证 MATLAB 桥接和双向显示。
html.Data=struct('kind','geometry-face-probe','revision',revision,'faceIndex',zeroIndex);
html.DataChangedFcn(html,[]);drawnow;
end

function frame=waitFrame(html)
timer=tic;
while toc(timer)<15
    drawnow;view=html.UserData;
    if ~isempty(view.error),error('%s',view.error.message);end
    if ~isempty(view.rendered) && view.rendered.revision==view.revision
        frame=view.rendered;
        assert(frame.faces==9216 && frame.webglError==0 && frame.visiblePixelSamples>0);
        assert(~isfield(view.payload,'faceField'));return
    end
    pause(.05);
end
error('未收到当前几何版本的 WebGL 回执。');
end
