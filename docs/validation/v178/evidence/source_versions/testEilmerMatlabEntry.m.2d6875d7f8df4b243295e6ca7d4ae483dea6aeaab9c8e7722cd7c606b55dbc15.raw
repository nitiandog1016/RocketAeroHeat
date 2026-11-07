function report = testEilmerMatlabEntry()
% 验证MATLAB入口的实际WSL参数传递、原结果回读及覆盖/错误来源拒绝，不重跑CFD。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root);
base=fullfile(root,'outputs','v178_validation');
out=fullfile(base,'matlab01');
assert(~isfolder(out),'已有测试不重复覆盖'); mkdir(out);
opt=struct('caseDirectory',fullfile(base,'native03'), ...
    'recordDirectory',fullfile(out,'read with spaces'), ...
    'loadsFile',fullfile(base,'native03','lmrsim','loads','0000','blk-0000-bndry-0-wall.dat'), ...
    'runReceipt',fullfile(base,'interface01','observed_run.json'),'seconds',30);
timer=tic;
r=runEilmerCase('read',opt);
assert(r.passed && r.wallFaces==7 && r.nativeStep==4 && r.newQualifiedAreaM2==0);
data=readtable(fullfile(opt.recordDirectory,'wall.csv'));
assert(isequal(data.originalFaceId,(1:7)'));
assert(any(data.heatFluxIntoWall_W_m2<0) && any(data.heatFluxIntoWall_W_m2>0));
assert(all(data.givenWallTemperature_K==300) && all(data.wallFaceFluidTemperature_K==300));
assert(abs(sum(data.area_m2.*data.heatFluxIntoWall_W_m2)-r.heatPowerIntoWallW)<2e-10);
rejected=0;
try, runEilmerCase('read',opt); catch e, assert(strcmp(e.identifier,'RocketAeroHeat:EilmerOverwrite'));rejected=rejected+1;end
bad=opt;bad.recordDirectory=fullfile(out,'failed_source');
bad.runReceipt=fullfile(base,'interface01','failed_run.json');
try,runEilmerCase('read',bad);catch e,assert(strcmp(e.identifier,'RocketAeroHeat:EilmerFailed'));rejected=rejected+1;end
assert(rejected==2);
% 原生准备与运行已由Python小域验证；MATLAB此处只做新的桥接行为，不重复计算。
report=struct('passed',true,'originalFaces',7,'rejected',rejected,'newNativeRuns',0, ...
    'elapsedSeconds',toc(timer),'newQualifiedAreaM2',0,'spaceInPathPassed',true);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
end
