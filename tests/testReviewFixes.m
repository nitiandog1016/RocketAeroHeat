function testReviewFixes()
% Run with: addpath('tests'); testReviewFixes
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
cfg = []; raw = []; curved = []; m = [];
evalc('cfg = createCaseConfig(); raw = loadSTLMesh(cfg); raw = calcFaceGeometry(raw);');
evalc('curved = analyzeSurfaceCurvature(raw,cfg);');
assert(all(isfinite(curved.targetSize)));
[raw,~] = alignGeometryAxis(raw,'pca');
evalc('m = autoRemeshSurface(raw,cfg); m = fixNormalDirection(m,cfg);');
assert(m.normalInfo.boundaryEdges == 0);
assert(m.normalInfo.nonManifoldEdges == 0);

limited = cfg;
limited.remesh.maxStations = 3;
limited.remesh.initialAxialSpacing = 0.5;
limited.remesh.minCircumPoints = 9;
limited.remesh.maxCircumPoints = 9;
evalc('m = autoRemeshSurface(raw,limited);');
assert(m.remeshInfo.numberSections == 3);
assert(m.remeshInfo.minCircumPoints == 9);
assert(m.remeshInfo.maxCircumPoints == 9);
assert(m.remeshInfo.boundaryEdges == 0);
assert(m.remeshInfo.nonManifoldEdges == 0);
fprintf('Numerical regression checks passed.\n');

% Exercise the real nested callbacks while meshGUI waits for input.
failure = [];
completed = false;
driver = timer('StartDelay',1,'ExecutionMode','fixedSpacing', ...
    'Period',1,'TimerFcn',@exerciseGUI);
cleanup = onCleanup(@()delete(driver)); %#ok<NASGU>
start(driver);
[~,~,returnedMesh,returnedQuality,action] = meshGUI(cfg,raw);
if ~isempty(failure)
    rethrow(failure);
end
assert(completed);
assert(strcmp(action,'cancel'));
assert(isempty(returnedMesh) && isempty(returnedQuality));
fprintf('GUI invalidation checks passed.\n');

    function exerciseGUI(~,~)
        fig = [];
        try
            confirm = findall(groot,'Text','确认：主轴沿 X 且鼻尖朝 +X');
            if isempty(confirm)
                return;
            end
            assert(isscalar(confirm));
            fig = ancestor(confirm,'figure');
            % Initial rendering also processes timers, before uiwait starts.
            if ~strcmp(fig.WaitStatus,'waiting')
                return;
            end
            stop(driver);
            invoke(confirm);
            generate = findall(fig,'Text','重新生成表面网格');
            proceed = findall(fig,'Text','确认网格并继续主程序');
            reset = findall(fig,'Text','恢复默认网格参数');
            grid = generate.Parent;
            controls = grid.Children;
            checked = 0;
            for k = 1:numel(controls)
                control = controls(k);
                if ~(isa(control,'matlab.ui.control.NumericEditField') || ...
                        isa(control,'matlab.ui.control.Slider'))
                    continue;
                end
                assert(~isempty(control.ValueChangedFcn));
                evalc('invoke(generate);');
                assert(strcmp(proceed.Enable,'on'));
                % Programmatic property changes do not fire ValueChangedFcn.
                control.Value = control.Value + 1;
                control.ValueChangedFcn(control,[]);
                assert(strcmp(proceed.Enable,'off'));
                checked = checked + 1;
                fprintf('GUI control %d passed.\n',checked);
            end
            assert(checked == 10);
            evalc('invoke(generate);');
            assert(strcmp(proceed.Enable,'on'));
            invoke(reset);
            assert(strcmp(proceed.Enable,'off'));
            completed = true;
        catch ME
            failure = ME;
        end
        if ~isempty(fig) && isvalid(fig)
            fig.CloseRequestFcn(fig,[]);
        else
            close all force;
        end
    end

    function invoke(button)
        button.ButtonPushedFcn(button,[]);
    end
end
