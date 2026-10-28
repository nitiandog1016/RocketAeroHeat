function report=testRocketLinearStudyBinding(folder,summaryFile)
%TESTROCKETLINEARSTUDYBINDING 拒绝错端点/错字节/伪合格研究，且不创建结果目录。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);assert(~exist(folder,'dir'));mkdir(folder);source=jsondecode(fileread(summaryFile));rejected=false(3,1);
for j=1:3
    bad=source;
    if j==1,bad.sourceCase='krylov15';elseif j==2,bad.sourceEvidenceSHA256=repmat('0',1,64);else,bad.productionQualified=true;end
    input=fullfile(folder,sprintf('bad%d.json',j));fid=fopen(input,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(bad));fclose(fid);
    target=fullfile(folder,sprintf('must_not_exist%d',j));
    try,exportRocketWallVelocityStudy(target,fullfile(root,'outputs','v83_validation'),'krylov17',input);
    catch ex,rejected(j)=strcmp(ex.identifier,'WallVelocityView:LinearStudy');end
    assert(rejected(j)&&~exist(target,'dir'));
end
report=struct('status','passed-three-source-or-qualification-rejections-before-export','rejected',rejected, ...
    'newActualNSEvaluations',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
save(fullfile(folder,'linear_study_binding_evidence.mat'),'report','summaryFile','-v7');
end
