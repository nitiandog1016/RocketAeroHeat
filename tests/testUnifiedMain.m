function testUnifiedMain()
%TESTUNIFIEDMAIN 校核当前默认入口；完整求解链路另由 testWorkflowWorkbench 覆盖。
% 模块一：保存测试将使用的工作区名称，创建自有对照窗口与无关变量。
% main 不得清除已有工作区或关闭其他窗口；测试结束只清理自己的窗口。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
names={'cfg','rawMeshOriginal','rawMesh','quality','mesh','freestream', ...
    'flight','projectRoot','appFigure','rocketAeroHeatTestUnrelated'};
previous=cell(size(names));existed=false(size(names));
for k=1:numel(names)
    existed(k)=evalin('base',sprintf('exist(''%s'',''var'')==1',names{k}));
    if existed(k),previous{k}=evalin('base',names{k});end
end
unrelated=uifigure('Name','RocketAeroHeat 入口测试自有对照窗口', ...
    'Visible','off');
marker=struct('token','preserve-unrelated-data','values',[17 29 41]);
unrelated.UserData=marker;
assignin('base','rocketAeroHeatTestUnrelated',marker);
existingFigures=findall(groot,'Type','figure');
cleanup=onCleanup(@()restoreTestState(unrelated,names,existed,previous)); %#ok<NASGU>

%% 模块二：真实执行默认 main，必须立即返回到未导入、未求解的五步骤工作台。
started=tic;evalin('base','main');entrySeconds=toc(started);
app=evalin('base','appFigure');
appCleanup=onCleanup(@()deleteIfValid(app)); %#ok<NASGU>
assert(isscalar(app)&&isgraphics(app)&&~any(app==existingFigures));
assert(~strcmp(app.WaitStatus,'waiting'), ...
    '默认入口不应使用 uiwait 阻塞五步骤工作台。');
state=app.UserData.workflowState;
assert(state.currentStep==1&&strcmp(state.status,'awaiting-model'));
assert(isempty(state.rawMesh)&&isempty(state.orientedMesh)&&isempty(state.mesh));
assert(isempty(state.result)&&isempty(state.freestream)&&~state.busy);
assert(state.revision==0&&state.solveId==0&&~state.boundaryConfirmed);
assert(~any([state.steps.completed])&&~any([state.steps(2:5).available]));
assert(evalin('base','isempty(rawMeshOriginal)&&isempty(rawMesh)&&isempty(mesh)&&isempty(freestream)'));
checkPreserved();checkOneAddedWindow();

%% 模块三：核对五个真实按钮的顺序、回调和共同宿主，初始前置条件必须明确。
labels={'输入模型','生成网格','添加边界条件','计算','查看结果'};
buttons=app.UserData.workflowButtons;
assert(numel(buttons)==5&&numel(findall(app,'Type','figure'))==1);
for k=1:5
    assert(isa(buttons(k),'matlab.ui.control.Button'));
    assert(startsWith(buttons(k).Text,sprintf('%d  %s',k,labels{k})));
    assert(isa(buttons(k).ButtonPushedFcn,'function_handle'));
    assert(ancestor(buttons(k),'figure')==app);
    assert(ancestor(buttons(k),'uipanel')==app.UserData.workflowSidebar);
    if k==1,assert(strcmp(buttons(k).Enable,'on'));
    else,assert(strcmp(buttons(k).Enable,'off'));end
end
buttons(1).ButtonPushedFcn(buttons(1),[]);drawnow;
assert(app.UserData.workflowState.currentStep==1&&isempty(app.UserData.workflowState.result));
assert(isgraphics(app.UserData.controls.load)&&isgraphics(app.UserData.controls.confirm));
assert(ancestor(app.UserData.controls.load,'figure')==app);
assert(strcmp(app.UserData.resultHost.Visible,'off'));
checkPreserved();checkOneAddedWindow();

%% 模块四：等待真实网页初始化回执；空模型不能显示旧几何或伪造的物理场。
view=app.UserData.geometryView;
assert(ancestor(view.html,'figure')==app);
webStarted=tic;
while toc(webStarted)<40
    drawnow;browser=view.html.UserData;
    if ~isempty(browser.error)
        error('RocketAeroHeat:TestMain:Browser','%s',browser.error.message);
    end
    if browser.ready&&~isempty(browser.rendered)&& ...
            browser.rendered.revision==browser.revision,break;end
    pause(.1);
end
browser=view.html.UserData;
assert(browser.ready&&~isempty(browser.rendered)&& ...
    browser.rendered.revision==browser.revision, ...
    '默认入口的几何网页没有收到当前版本初始化回执。');
assert(browser.rendered.webglError==0&&browser.rendered.width>0&&browser.rendered.height>0);
assert(browser.rendered.vertices==0&&browser.rendered.faces==0&&browser.rendered.visiblePixelSamples==0);
assert(isempty(browser.payload.vertices)&&isempty(browser.payload.faces)&& ...
    ~isfield(browser.payload,'fields'));
state=app.UserData.workflowState;
assert(isempty(state.rawMesh)&&isempty(state.result)&&state.solveId==0);
checkPreserved();checkOneAddedWindow();
fprintf(['Unified default main entry passed: nonblocking startup %.3f s, ' ...
    '5 actual ordered sidebar buttons, one new host, no automatic model/solve, ' ...
    'preserved unrelated workspace/windows, empty WebGL initialized (GL=0).\n'],entrySeconds);

    function checkPreserved()
        assert(all(isgraphics(existingFigures)), ...
            'main 不得关闭启动前已经存在的窗口。');
        assert(isgraphics(unrelated)&&isequaln(unrelated.UserData,marker));
        assert(isequaln(evalin('base','rocketAeroHeatTestUnrelated'),marker), ...
            'main 不得清除已有的无关变量。');
    end
    function checkOneAddedWindow()
        current=findall(groot,'Type','figure');
        assert(numel(current)==numel(existingFigures)+1&&any(current==app), ...
            '默认入口及切页只能增加同一个工作台窗口。');
    end
end

function restoreTestState(unrelated,names,existed,previous)
% 模块五：使用独立清理函数，失败退出时也不引用已结束的共享工作区。
deleteIfValid(unrelated);
for j=1:numel(names)
    if existed(j),assignin('base',names{j},previous{j});
    else,evalin('base',['clear ' names{j}]);end
end
end
function deleteIfValid(fig)
if isgraphics(fig),delete(fig);end
end
