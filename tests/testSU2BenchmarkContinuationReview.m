function report=testSU2BenchmarkContinuationReview(out,source)
%TESTSU2BENCHMARKCONTINUATIONREVIEW 核对实际公开续算及原工作台导出的完整来源。
% 模块一，主五步骤载入真实审核，检查源步号、本段与累计次数。
% 模块二，逐站比对导出与方法，拒绝错误累计身份并核对旧初场记录。
assert(~isfolder(out));mkdir(out);fig=rocketAeroHeatWorkbench();guard=onCleanup(@()close(fig));
before=fig.UserData.workflowState;
fig.UserData.actions.readBenchmarkReview(source);drawnow;
view=fig.UserData.benchmarkReview;review=view.review;
assert(review.sourceLocalIteration==1199 && review.continuation.completedUpdatesBeforeRun==1200);
assert(review.completedUpdatesInSeries==1200+review.lastRecordedIteration+1);
assert(contains(view.summary.Text,sprintf('系列 %d / 12000',review.completedUpdatesInSeries)));
assert(~review.physicsQualified && review.newQualifiedArea==0);
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.engineJob));
assert(size(view.gateTable.Data,1)==26);
file=fullfile(out,'benchmark_review.json');view.exportTo(file);
assert(strcmp(fileread(file),fileread(source)));
if ~isempty(review.currentFieldAssessment)
    data=readtable([file '.stations.csv'],'TextType','string');
    assert(height(data)==51 && all(data.SourceLocalIteration==1199));
    assert(all(data.CompletedUpdatesInSeries==review.completedUpdatesInSeries));
    assert(all(data.SavedIteration==review.savedLocalIteration));
    assert(all(data.StabilitySampling==string(review.temporalStability.sampling)));
    assert(all(data.Method==string(review.model)));
    names={'pressure','heat'};
    for j=1:2
        actual=data(data.Quantity==string(names{j}),:);expected=review.currentFieldAssessment.stations.(names{j});
        assert(max(abs(actual.Predicted-expected.predicted)./max(abs(expected.predicted),1))<1e-12);
    end
end
bad=jsondecode(fileread(source));bad.completedUpdatesInSeries=bad.completedUpdatesInSeries+1;
badfile=fullfile(out,'bad_lineage.json');write(badfile,bad);rejected=false;
try,readSU2BenchmarkReview(badfile);catch,rejected=true;end
assert(rejected);
legacy=readSU2BenchmarkReview(fullfile(pwd,'outputs','v151_validation','review01','report.json'));
assert(isempty(legacy.sourceLocalIteration) && legacy.savedLocalIteration==1199);
fig.UserData.actions.selectStep(1);drawnow;
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.benchmarkReview));
report=struct('passed',true,'sourceLocalIteration',1199, ...
    'savedLocalIteration',review.savedLocalIteration,'completedUpdatesInSeries',review.completedUpdatesInSeries, ...
    'gateCount',26,'wrongLineageRejected',rejected,'legacyIndependentStatePassed',true, ...
    'newMouseOrScreenshotCheck',false);
write(fullfile(out,'report.json'),report);
end
function write(file,value)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(value),'char');
end
