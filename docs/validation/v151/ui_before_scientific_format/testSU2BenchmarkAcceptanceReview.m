function report=testSU2BenchmarkAcceptanceReview(out,source)
%TESTSU2BENCHMARKACCEPTANCEREVIEW 从主入口核对真实完整场、指标、曲线和旧记录。
% 模块一：绑定局部1199及独立初场，26项门槛和未获物理资格分别保存。
% 模块二：审核导出测站/单位/方法、伪造通过拒绝，以及返回原输入。
assert(~isfolder(out));mkdir(out);fig=rocketAeroHeatWorkbench();guard=onCleanup(@()close(fig));
before=fig.UserData.workflowState;
fig.UserData.actions.readBenchmarkReview(source);drawnow;
view=fig.UserData.benchmarkReview;review=view.review;
assert(review.exitCode==0 && review.lastRecordedIteration==1199 && review.savedLocalIteration==1199);
assert(isempty(review.sourceLocalIteration) && review.hasSavedNewField);
assert(~review.validationPassed && ~review.physicsQualified && review.convergenceAcceptanceImplemented);
assert(contains(view.summary.Text,'独立相容初场'));
assert(size(view.gateTable.Data,1)==26 && sum(strcmp(view.gateTable.Data(:,4),'未通过'))==21);
assert(~isempty(view.comparisonTab));view.tabs.SelectedTab=view.comparisonTab;drawnow;
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.engineJob));
file=fullfile(out,'benchmark_review.json');view.exportTo(file);
assert(strcmp(fileread(file),fileread(source)));
data=readtable([file '.stations.csv'],'TextType','string');
assert(height(data)==51 && sum(data.Quantity=="pressure")==24 && sum(data.Quantity=="heat")==27);
assert(all(data.SavedIteration==1199) && all(~data.ValidationPassed));
assert(all(data.Method==string(review.model)));
bad=jsondecode(fileread(source));bad.acceptance.gates(3).passed=true;
badfile=fullfile(out,'bad_gate.json');write(badfile,bad);rejected=false;
try,readSU2BenchmarkReview(badfile);catch,rejected=true;end
assert(rejected);
legacy=readSU2BenchmarkReview(fullfile(pwd,'outputs','v150_validation','review01','report.json'));
assert(legacy.exitCode==1 && ~legacy.hasSavedNewField && legacy.sourceLocalIteration==1);
write(fullfile(out,'ready_for_capture.json'),struct('ready',true));
t=tic;
while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<120,drawnow;pause(.2);end
captured=isfile(fullfile(out,'capture_done.json'));
fig.UserData.actions.selectStep(1);drawnow;
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.benchmarkReview));
report=struct('mainWorkbenchEntry',true,'independentInitialIdentity',true,'currentFieldIteration',1199, ...
    'gateCount',26,'failedGateCount',21,'stationRows',height(data),'falseGateRejected',rejected, ...
    'legacyV150Accepted',true,'inputStatePreserved',true,'exportExact',true, ...
    'nativeScreenshotObserved',captured,'solverRuns',0,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
