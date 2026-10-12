function runAllTests()
%RUNALLTESTS 独立执行全部 test*.m，并汇总每项结果。
% 模块作用：每个测试失败后仍继续下一项；最后统一报告失败清单。
% 使用函数工作区，避免 main.m 在 MATLAB base 工作区的变量操作
% 干扰测试调度器自身。

%% 1. 发现测试并加入项目路径
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fileparts(mfilename('fullpath')));
files=dir(fullfile(root,'tests','test*.m'));
failures={};

%% 2. 逐项执行，保留可复现的通过/失败名称
for index=1:numel(files)
    name=erase(files(index).name,'.m');
    try
        feval(name);
        fprintf('PASS %s\n',name);
    catch problem
        fprintf(2,'FAIL %s: %s\n',name,problem.message);
        % 保留实际失败行和调用栈，避免只记录“断言失败”而无法定位。
        fprintf(2,'%s\n',getReport(problem,'extended','hyperlinks','off'));
        failures{end+1}=name; %#ok<AGROW>
    end
end

%% 3. 汇总，不让部分通过冒充完整通过
fprintf('RESULT %d/%d passed\n',numel(files)-numel(failures),numel(files));
if ~isempty(failures)
    error('RocketAeroHeat:TestSuite:Failure', ...
        '失败测试：%s',strjoin(failures,'、'));
end
end
