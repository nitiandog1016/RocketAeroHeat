function check_source_cache(out)
% 模块用途：在独立目录从固定作者提交首次获取、离线回读，并检验损坏文件的拒收。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
app=fullfile(out,'app');caseDir=fullfile(app,'cases','rfz_rec');mkdir(caseDir);
for file=["readRFZReferenceSamples.m","cm011FileSHA256.m"],copyfile(fullfile(root,file),app);end
for file=["deg190_reference.json","deg190_boundary.json"],copyfile(fullfile(root,'cases','rfz_rec',file),caseDir);end
prior=pwd;restore=onCleanup(@()cd(prior));cd(app); %#ok<NASGU>
assert(strcmp(which('readRFZReferenceSamples'),fullfile(app,'readRFZReferenceSamples.m')));
first=readRFZReferenceSamples();second=readRFZReferenceSamples();
assert(height(first)==72706 && isequal(first,second));
source=jsondecode(fileread(fullfile(caseDir,'deg190_reference.json')));
cache=fullfile(app,'outputs','rfz_reference','published_samples',source.sourceCommit);
files=dir(fullfile(cache,'*_data.npy'));assert(numel(files)==6 && all([files.bytes]==581776));
bad=fullfile(out,'bad');mkdir(bad);
f=fopen(fullfile(bad,'X_data.npy'),'wb');assert(f>=0);fwrite(f,uint8('wrong source'));fclose(f);
rejected=false;
try,readRFZReferenceSamples(bad);catch ex,rejected=strcmp(ex.identifier,'RAH:RFZ:Reference');end
assert(rejected);
rec=struct('files',numel(files),'rawBytes',sum([files.bytes]),'samples',height(first), ...
    'coldAndCachedIdentical',true,'corruptSourceRejected',rejected);
f=fopen(fullfile(out,'cache.json'),'w');assert(f>=0);fprintf(f,'%s\n',jsonencode(rec,'PrettyPrint',true));fclose(f);
end
