function az=readShortSupportAzimuth(cases,job,s,stations,folder)
% 模块一：只给匹配源帧、面网和条件的 14D 两半步结果加载曲线，新任务使用自己的反力包。
az=[];local=false;
if ~isempty(job)
    receipt=jsondecode(fileread(fullfile(job,'completed.json')));
    local=isfield(receipt,'stationAzimuthSHA256');
end
if local
    file=fullfile(job,'export','stationAzimuth.zip');key=receipt.stationAzimuthSHA256;
else
    if ~strcmp(s.frameSHA256,'0132e9b266ef344fcbf13fd9a7d3dd994fa56ef3192caa3fbf7c8fc1da5f002e'),return;end
    file=fullfile(cases,'stationAzimuth65.zip');key='1c4167c5f2d314419b2db464fe13e449e9ba27bc58305bb9617caa82b47c5cc3';
end
assert(strcmp(cm011FileSHA256(file),key),'RAH:Short:Azimuth','周向曲线包与回执不符。');
% 使用独立临时目录，避免与主面表的同名摘要冲突。
scratch=tempname(folder);mkdir(scratch);clean=onCleanup(@()discard(scratch)); %#ok<NASGU>
unzip(file,scratch);a=jsondecode(fileread(fullfile(scratch,'summary.json')));
assert(strcmp(a.schema,'short65-azimuth-v1') && strcmp(a.frameSHA256,s.frameSHA256) && ...
    strcmp(a.meshSHA256,s.meshSHA256) && strcmp(a.caseSHA256,s.caseSHA256) && strcmp(a.bodySHA256,s.bodySHA256) && ...
    strcmp(a.reactionSHA256,s.releaseSHA256) && strcmp(a.wallBudgetSHA256,s.energy.files.wall_budget_npz) && ...
    strcmp(a.wallSHA256,s.files.wall_csv) && ...
    a.elapsedS==s.elapsedSeconds && a.dtS==s.dtSeconds,'RAH:Short:Azimuth','周向曲线与所选时刻的壁面结果不符。');
% 模块二：未解区保持 NaN。保留匹配到的标准三角面号，导入网的编号在映射完成后填写。
t=readtable(fullfile(scratch,'points.csv'),'TextType','string');valid=t.status=="sphere";
assert(isequal(unique(t.angleDeg),stations.angleDeg) && all(isfinite(t.nu(valid))) && all(isnan(t.nu(~valid))));
assert(all(isfinite(t.qInWm2(valid))) && all(isnan(t.qInWm2(~valid))) && all(t.faceID(valid)>0 & t.faceID(valid)<=s.faces));
zero=t(t.phiDeg==0,:);assert(isequal(zero.angleDeg,stations.angleDeg) && max(abs(zero.nu-stations.computedNu))<1e-8);
t.referenceFaceID=t.faceID;t.fileFaceID=t.faceID;t=removevars(t,'faceID');
az=struct('summary',a,'points',t,'archiveSHA256',key,'stlSHA256','');
end
function discard(folder)
for name={'points.csv','summary.json'}
    f=fullfile(folder,name{1});if isfile(f),delete(f);end
end
rmdir(folder);
end
