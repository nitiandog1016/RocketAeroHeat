function report=testRocketDirectionPreparationBinding(folder)
%TESTROCKETDIRECTIONPREPARATIONBINDING 拒绝把其他结果组或被改源文件的方法贴到当前工作台。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);actual=fullfile(root,'outputs','v83_validation','peakweightns01','peak_direction_ns_summary.json');
s=jsondecode(fileread(actual));rejections=0;
for j=1:2
    bad=s;if j==1,bad.candidateCase='krylov22';else,bad.candidateFileSHA256=repmat('0',1,64);end
    file=fullfile(folder,sprintf('bad_source_%d.json',j));fid=fopen(file,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(bad));fclose(fid);
    rejected=false;destination=fullfile(folder,sprintf('must_not_export_%d',j));
    try,exportRocketWallVelocityStudy(destination,fullfile(root,'outputs','v83_validation'),'krylov24',[],file);
    catch ex,rejected=strcmp(ex.identifier,'WallVelocityView:DirectionSource');end
    assert(rejected&&~exist(destination,'dir'));rejections=rejections+1;
end
report=struct('status','passed-direction-preparation-wrong-case-and-source-SHA-rejections','rejections',rejections, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'direction_preparation_binding_evidence.mat'),'report','-v7');fprintf('%s\n',jsonencode(report));
end
