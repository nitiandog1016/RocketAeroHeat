function testLegacyUnifiedMain()
% 从 main 入口驱动真实确认按钮，核对整个流程只使用同一个窗口。
addpath(fileparts(fileparts(mfilename('fullpath'))));
failure=[]; seen=[]; stage=1; started=tic;
driver=timer('ExecutionMode','fixedSpacing','Period',0.5,'TimerFcn',@drive);
cleanup=onCleanup(@()delete(driver)); %#ok<NASGU>
start(driver);
evalin('base','mainLegacy');
stop(driver);
if ~isempty(failure), rethrow(failure); end
assert(stage==3 && isgraphics(seen));
assert(evalin('base','isfield(mesh,''pressure'')'));
assert(evalin('base','isfield(mesh,''externalFlow'') && isfield(mesh,''surfaceStreamlines'') && isfield(mesh,''heating'')'));
assert(evalin('base','mesh.axisymmetricMarch.applied && mesh.boundaryLayerEdge.summary.axialMarchFaces>10000'));
assert(evalin('base','contains(mesh.boundaryLayerEdge.method,''内置工程外流'')'));
assert(evalin('base','mesh.stagnationHeating.available && mesh.stagnationHeating.radius>0'));
assert(evalin('base','mesh.heating.transition.pathUsed'));
ui=seen.UserData;
ui.heatChoice.Value=9;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequaln(ui.heatSurface.FaceVertexCData, ...
    ui.boundaryLayerEdge.shockAngleDeg));
assert(numel(findall(groot,'Type','figure'))==1);
delete(seen);
fprintf('Unified main entry passed: one figure throughout the complete main pipeline.\n');
    function drive(~,~)
        figs=findall(groot,'Type','figure');
        if isempty(figs), return; end
        fig=figs(1);
        if toc(started)>150, failure=MException('Test:Timeout','main 超时'); delete(fig); return; end
        if ~strcmp(fig.WaitStatus,'waiting'), return; end
        try
            if isempty(seen), seen=fig; else, assert(seen==fig); end
            assert(numel(figs)==1);
            if stage==1
                invoke(fig,'自动识别主轴 → X');
                invoke(fig,'确认：主轴沿 X 且鼻尖朝 +X');
                invoke(fig,'重新生成表面网格');
                invoke(fig,'确认网格并继续主程序'); stage=2;
            elseif stage==2
                invoke(fig,'确认飞行条件并继续'); stage=3;
            end
        catch ME
            failure=ME; delete(fig);
        end
    end
end
function invoke(fig,text)
button=findall(fig,'Text',text);
keep=true(size(button));
for i=1:numel(button)
    parent=button(i);
    while ~isequal(parent,fig)
        if isprop(parent,'Visible') && strcmp(parent.Visible,'off'), keep(i)=false; break; end
        parent=parent.Parent;
    end
end
button=button(keep); assert(isscalar(button));
button.ButtonPushedFcn(button,[]);
end
