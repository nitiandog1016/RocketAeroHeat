function samples=readRFZReferenceSamples(sourceDirectory)
% 模块一：固定作者的工况、提交和六个 NPY 摘要，保留原始样点顺序及符号。
root=fileparts(mfilename('fullpath'));
caseDir=fullfile(root,'cases','rfz_rec');
boundaryFile=fullfile(caseDir,'deg190_boundary.json');
assert(strcmp(cm011FileSHA256(boundaryFile),'b6647eb6dd7cb68e8b0db272240bdae0a9f4f686c660f63cf1eb62e1758de073'), ...
    'RAH:RFZ:Reference','参考字段清单与固定版本不符。');
boundary=jsondecode(fileread(boundaryFile));
source=jsondecode(fileread(fullfile(caseDir,'deg190_reference.json')));
assert(strcmp(boundary.dataset,'EUCASS_2023/deg190') && strcmp(boundary.dataset,source.dataset) && ...
    strcmp(cm011FileSHA256(fullfile(caseDir,'deg190_reference.json')),boundary.caseSHA256), ...
    'RAH:RFZ:Reference','工况身份不匹配。');
hashes=boundary.fieldSHA256;
hashes.X='e511731392daafa3e97379a6a40c3dc307ce26705519fc1a843d5327eca4cc8d';
hashes.Y='62f415a45076b7604d9feac76886ede6a3073678c606a3266139e4b219b59a02';
hashes.Z='4e11f63efb2a4495a263a71ab049bce3a4f2415ef12966ca6df25d079cd01c38';
relative=['Aerodynamics/Stage_1/' source.dataset];
local=fullfile(root,'outputs','rfz_reference');
archive=fullfile(local,'author_data_20261124','files',relative);
cache=fullfile(local,'published_samples',source.sourceCommit);
keys={'X','Y','Z','p','q','T'};
n=boundary.referencePoints;values=zeros(n,numel(keys));
for i=1:numel(keys)
    key=keys{i};file=[key '_data.npy'];
    if nargin>0
        path=fullfile(sourceDirectory,file);
    elseif isfile(fullfile(archive,file))
        path=fullfile(archive,file);
    else
        path=fullfile(cache,file);
        if ~isfile(path)
            url=['https://raw.githubusercontent.com/rfz-model/RFZ_data/' source.sourceCommit '/' relative '/' file];
            fetchSample(url,path,hashes.(key),n);
        end
    end
    values(:,i)=readVector(path,hashes.(key),n);
end
samples=array2table([(1:n)' values],'VariableNames', ...
    {'sampleID','authorX_m','authorY_m','authorZ_m','authorPressure_Pa','authorHeatFlux_Wm2','authorWallTemperature_K'});
end

function v=readVector(path,expected,n)
% 模块二：在同一份字节上校验来源并解析固定的一维小端双精度数组。
[fid,msg]=fopen(path,'rb');assert(fid>=0,'RAH:RFZ:Reference','读取失败 %s，%s',path,msg);
closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
b=fread(fid,128+8*n+1,'*uint8');
hash=java.security.MessageDigest.getInstance('SHA-256');hash.update(typecast(b,'int8'));
actual=lower(reshape(dec2hex(typecast(hash.digest(),'uint8'),2).',1,[]));
assert(strcmp(actual,expected) && numel(b)==128+8*n && isequal(b(1:10)',uint8([147 double('NUMPY') 1 0 118 0])), ...
    'RAH:RFZ:Reference','样点文件或字段顺序与作者版本不符：%s',path);
[~,~,byteOrder]=computer;
v=typecast(b(129:end),'double');if byteOrder=='B',v=swapbytes(v);end
assert(all(isfinite(v)),'RAH:RFZ:Reference','样点包含非有限值。');
end

function fetchSample(url,path,expected,n)
% 模块三：首次打开按需取得缺少的原文件。校验通过后发布到离线缓存。
folder=fileparts(path);if ~isfolder(folder),mkdir(folder);end
partial=[tempname(folder) '.npy'];
cleanup=onCleanup(@()removePartial(partial)); %#ok<NASGU>
websave(partial,url,weboptions('Timeout',45,'ContentType','binary'));
readVector(partial,expected,n);
if ~isfile(path)
    [ok,msg]=movefile(partial,path);
    assert(ok || isfile(path),'RAH:RFZ:Reference','无法保存参考数据：%s',msg);
end
end

function removePartial(path)
if isfile(path),delete(path);end
end
