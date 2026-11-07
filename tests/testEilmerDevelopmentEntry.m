function report=testEilmerDevelopmentEntry()
% 模块一：MATLAB真实传递发展配置/含空格路径，只准备小域，不重跑CFD。
% 模块二：读取真实内存失败，明确返回failed和无壁面数据，不能被exit0误导。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
base=fullfile(root,'outputs','v180_validation');out=fullfile(base,'matlab01');
assert(~isfolder(out));mkdir(out);timer=tic;
condition=jsondecode(fileread(fullfile(base,'small-condition.json')));
opt=struct('caseDirectory',fullfile(out,'prepared with spaces'), ...
    'recordDirectory',fullfile(out,'prepare receipt'), ...
    'sourceCase',fullfile(root,'outputs','v178_validation','native03'), ...
    'sourcePreparedReceipt',fullfile(root,'outputs','v178_validation','native03_prepare','report.json'), ...
    'sourceGeometryAudit',fullfile(root,'outputs','v178_validation','native03_geometry.json'), ...
    'condition',condition,'seconds',30,'outputMiB',128);
prepared=runEilmerCase('development',opt);
assert(prepared.passed && prepared.uniformInitialFieldBitwiseIdentical);
opt=struct('caseDirectory',fullfile(base,'real01'),'recordDirectory',fullfile(out,'failed inspection'), ...
    'runReceipt',fullfile(base,'real01_run','report.json'),'seconds',20);
failed=runEilmerCase('inspect',opt);
assert(~failed.passed && ~failed.solverCompleted && ~failed.completedWallDataAvailable);
assert(failed.launcherExitCode==0 && contains(failed.failureReason,'OutOfMemoryError'));
assert(failed.newQualifiedAreaM2==0);
report=struct('passed',true,'developmentPreparationPassed',true,'realNativeFailureIdentified',true, ...
    'newCFDRuns',0,'spaceInPathPassed',true,'elapsedSeconds',toc(timer));
f=fopen(fullfile(out,'report.json'),'w','n','UTF-8');closer=onCleanup(@()fclose(f));fprintf(f,'%s\n',jsonencode(report));
end
