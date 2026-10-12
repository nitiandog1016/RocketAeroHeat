function testUnifiedGUI(out)
% 验证同一窗口贯穿网格、飞行、压力及返回修改；调用真实控件回调。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig(); raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);');
cfg.output.verbose=false;
[raw,~]=alignGeometryAxis(raw,'pca');
[fig,host]=createAeroHeatWindow();
cleanup=onCleanup(@()safeDelete(fig)); %#ok<NASGU>
initialFigures=numel(findall(groot,'Type','figure'));
phase='mesh'; failure=[]; started=tic;
driver=timer('ExecutionMode','fixedSpacing','Period',0.5,'TimerFcn',@drive);
timerCleanup=onCleanup(@()delete(driver)); %#ok<NASGU>
start(driver);
[cfg,~,mesh,~,action]=meshGUI(cfg,raw,host);
assert(isempty(failure),messageOf(failure)); assert(strcmp(action,'continue'));
assert(isgraphics(fig) && numel(findall(groot,'Type','figure'))==initialFigures);
evalc('mesh=fixNormalDirection(mesh,cfg); mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
phase='flight'; setAeroHeatStage(fig,2);
[flight,fs,action]=createFlightCondition(cfg,mesh,host);
assert(isempty(failure),messageOf(failure)); assert(strcmp(action,'continue'));
assert(flight.alpha==20); cfg.flight=flight;
mesh=modifiedNewtonianPressure(mesh,fs,cfg);
same=pressureGUI(mesh,cfg,host); assert(same==fig);
assert(numel(findall(groot,'Type','figure'))==initialFigures);
if nargin>0, drawnow; exportapp(fig,fullfile(out,'unified_pressure.png')); end
phase='flightAgain';
button=fig.UserData.changeFlight; button.ButtonPushedFcn(button,[]);
assert(isempty(failure),messageOf(failure));
updated=evalin('base','mesh');
assert(updated.pressure.freestream.alpha==30);
assert(isgraphics(fig) && numel(findall(groot,'Type','figure'))==initialFigures);
phase='cancel';
button=fig.UserData.changeFlight; button.ButtonPushedFcn(button,[]);
assert(isgraphics(fig));
updated=evalin('base','mesh'); assert(updated.pressure.freestream.alpha==30);
phase='close';
button=fig.UserData.changeFlight; button.ButtonPushedFcn(button,[]);
assert(~isgraphics(fig));
assert(isempty(failure),messageOf(failure));
stop(driver);
fprintf('Unified GUI passed: one window; mesh -> flight -> pressure; recalculate, cancel and close.\n');

    function drive(~,~)
        if ~isgraphics(fig), return; end
        if toc(started)>180
            failure=MException('Test:Timeout','GUI 超时'); delete(fig); return;
        end
        if ~strcmp(fig.WaitStatus,'waiting'), return; end
        try
            if strcmp(phase,'mesh')
                invoke('确认：主轴沿 X 且鼻尖朝 +X');
                invoke('重新生成表面网格');
                invoke('确认网格并继续主程序');
                phase='between';
            elseif strcmp(phase,'flight') || strcmp(phase,'flightAgain')
                label=activeMatches(fig,'攻角 alpha / deg');
                assert(isscalar(label)); fields=label.Parent.Children;
                for k=1:numel(fields)
                    field=fields(k);
                    if isa(field,'matlab.ui.control.NumericEditField') && isequal(field.Layout.Row,4)
                        if strcmp(phase,'flight'), field.Value=20; else, field.Value=30; end
                        field.ValueChangedFcn(field,[]);
                    end
                end
                invoke('确认飞行条件并继续'); phase='between';
            elseif strcmp(phase,'cancel')
                invoke('取消主程序'); phase='between';
            elseif strcmp(phase,'close')
                fig.CloseRequestFcn(fig,[]);
            end
        catch ME
            failure=ME; safeDelete(fig);
        end
    end
    function invoke(text)
        button=activeMatches(fig,text); assert(isscalar(button));
        button.ButtonPushedFcn(button,[]);
    end
end
function objects=activeMatches(fig,text)
objects=findall(fig,'Text',text);
keep=true(size(objects));
for i=1:numel(objects)
    parent=objects(i);
    while ~isequal(parent,fig)
        if isprop(parent,'Visible') && strcmp(parent.Visible,'off'), keep(i)=false; break; end
        parent=parent.Parent;
    end
end
objects=objects(keep);
end
function safeDelete(fig)
if isgraphics(fig), delete(fig); end
end
function text=messageOf(failure)
text='GUI 测试失败';
if ~isempty(failure), text=getReport(failure,'extended','hyperlinks','off'); end
end
