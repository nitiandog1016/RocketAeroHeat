function report=testSU2BenchmarkReview(out,source)
%TESTSU2BENCHMARKREVIEW 从主五步骤入口查看真实失败，验证身份、导出和返回。
assert(~isfolder(out));mkdir(out);fig=rocketAeroHeatWorkbench();guard=onCleanup(@()close(fig));
before=fig.UserData.workflowState;assert(isfield(fig.UserData.controls,'readBenchmarkReview'));
fig.UserData.actions.readBenchmarkReview(source);drawnow;
view=fig.UserData.benchmarkReview;review=view.review;
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.engineJob));
assert(review.exitCode==1 && review.lastRecordedIteration==325 && review.sourceLocalIteration==1);
assert(~review.hasSavedNewField && isempty(review.savedLocalIteration) && ~review.physicsQualified);
assert(contains(view.heading.Text,'SST-Vm') && contains(view.heading.Text,'失败'));
assert(contains(view.summary.Text,'W/m') && contains(view.summary.Text,'尚无新保存场'));
assert(isfile(view.picture.ImageSource) && strcmp(view.picture.ScaleMethod,'fit'));
exportapp(fig,fullfile(out,'workbench.png'));
view.exportTo(fullfile(out,'benchmark_review.json'));
assert(strcmp(fileread(source),fileread(fullfile(out,'benchmark_review.json'))));
% 伪造最后历史步为新保存场，应在进入查看器前拒绝。
bad=jsondecode(fileread(source));bad.hasSavedNewField=true;bad.savedLocalIteration=325;
file=fullfile(out,'bad_saved_frame.json');write(file,bad);rejected=false;
try,readSU2BenchmarkReview(file);catch,rejected=true;end
assert(rejected);
fig.UserData.actions.selectStep(1);drawnow;
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.benchmarkReview));
assert(strcmp(fig.UserData.resultHost.Visible,'off'));
report=struct('mainWorkbenchEntry',true,'actualFailureShown',true,'inputStatePreserved',true, ...
    'fakeSavedFrameRejected',rejected,'exportExact',true,'returnedToModelStep',true, ...
    'lastRecordedIteration',325,'sourceLocalIteration',1,'savedNewField',false, ...
    'solverRuns',0,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
