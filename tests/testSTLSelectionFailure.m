function testSTLSelectionFailure()
% A failed file selection must leave the old configuration and geometry intact.
addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg = []; raw = [];
evalc('cfg=createCaseConfig(); raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);');
mockFolder = tempname;
mkdir(mockFolder);
file = fopen(fullfile(mockFolder,'uigetfile.m'),'w');
assert(file >= 0);
fprintf(file,['function [name,folder]=uigetfile(varargin)\n' ...
    'name=''missing_review_test.stl'';\n' ...
    'folder=fileparts(mfilename(''fullpath''));\nend\n']);
fclose(file);
cleanupMock = onCleanup(@()removeMock(mockFolder)); %#ok<NASGU>
addpath(mockFolder,'-begin');
failure = [];
completed = false;
driver = timer('StartDelay',1,'ExecutionMode','fixedSpacing','Period',1, ...
    'TimerFcn',@exerciseGUI);
cleanupTimer = onCleanup(@()delete(driver)); %#ok<NASGU>
start(driver);
[returnedCfg,returnedRaw,~,~,action] = meshGUI(cfg,raw);
if ~isempty(failure), rethrow(failure); end
assert(completed && strcmp(action,'cancel'));
assert(strcmp(returnedCfg.geometry.stlFile,cfg.geometry.stlFile));
assert(strcmp(returnedCfg.geometry.stlName,cfg.geometry.stlName));
assert(isequal(returnedRaw.V,raw.V) && isequal(returnedRaw.F,raw.F));
fprintf('Failed STL selection regression passed.\n');

    function exerciseGUI(~,~)
        button = findall(groot,'Text','选择新的 STL 文件');
        if isempty(button), return; end
        fig = ancestor(button,'figure');
        if ~strcmp(fig.WaitStatus,'waiting'), return; end
        stop(driver);
        try
            evalc('button.ButtonPushedFcn(button,[]);');
            pathField = findall(fig,'Value',cfg.geometry.stlFile);
            assert(~isempty(pathField));
            completed = true;
        catch ME
            failure = ME;
        end
        fig.CloseRequestFcn(fig,[]);
    end
end

function removeMock(folder)
rmpath(folder);
% Delete only the two temporary items created by this test.
delete(fullfile(folder,'uigetfile.m'));
rmdir(folder);
end
