function report=testEilmerExplicitEntry()
% 模块一：通过现有MATLAB入口读取显式研究结果，不重跑任何CFD。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
base=fullfile(root,'outputs','v183_validation');out=fullfile(base,'matlab01');
assert(~isfolder(out));mkdir(out);timer=tic;
opt=struct('caseDirectory',fullfile(base,'small01'), ...
    'recordDirectory',fullfile(out,'inspection receipt'), ...
    'runReceipt',fullfile(base,'small01_recovery','report.json'),'seconds',20);
result=runEilmerCase('inspect',opt);
assert(result.passed && result.solverCompleted && ~result.physicsQualified);
assert(result.nativeCompletion.finalStep==6);
assert(strcmp(result.nativeCompletion.timeMeaning,'local-pseudo-time-not-flight-time'));
% 模块二：有符号热流、原面及给定壁温通过原来的read通路导出。
opt.recordDirectory=fullfile(out,'read receipt');
opt.loadsFile=fullfile(base,'small01','lmrsim','loads','0001','blk-0000-bndry-0-wall.dat');
wall=runEilmerCase('read',opt);
assert(wall.passed && wall.wallFaces==7 && wall.newQualifiedAreaM2==0);
assert(wall.negativeHeatFaces==5 && wall.positiveHeatFaces==2);
assert(strcmp(wall.numericsProfile,'rah-explicit-local-v1'));
assert(strcmp(wall.timeMeaning,'local-pseudo-time-not-flight-time'));
% 模块三：旧失败回执仍按失败展示，不能被修正后的解析器自动变成成功。
opt.recordDirectory=fullfile(out,'original failed receipt');
opt.runReceipt=fullfile(base,'small01_run','report.json');
failed=runEilmerCase('inspect',opt);
assert(~failed.passed && ~failed.solverCompleted && ~failed.completedWallDataAvailable);
report=struct('passed',true,'inspectAndReadPassed',true,'originalFailurePreserved',true, ...
    'newCFDRuns',0,'newQualifiedAreaM2',0,'elapsedSeconds',toc(timer));
f=fopen(fullfile(out,'report.json'),'w','n','UTF-8');closer=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(report));
end
