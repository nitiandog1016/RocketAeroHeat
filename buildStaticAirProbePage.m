function page=buildStaticAirProbePage(evidence)
%BUILDSTATICAIRPROBEPAGE 把已审核的静止空气节点和来源封装为可离线打开的同一页面。
% 模块一：从固定清单读取并校验同一份原始字节，再解码为页面数据。
root=fileparts(mfilename('fullpath'));
if nargin==0,evidence=fullfile(root,'docs','validation','v261','evidence');end
manifest=verifiedText(fullfile(evidence,'manifest.json'), ...
    'd2efc5a87430bf8b44d86d056f28b4580d66bc0b39b92d814f8482c557cb2f5f');
hashes=jsondecode(manifest);
data=struct;
for group={'default','stable'}
    id=group{1};sample=struct;
    for entry={'case.cfg','initialization.json','initial.csv','terminal.csv'}
        name=entry{1};rel=[id '/' name];
        key=matlab.lang.makeValidName(rel);
        digest=hashes.(key);
        item=struct('path',rel,'sha256',digest, ...
            'text',verifiedText(fullfile(evidence,id,name),digest));
        [~,stem]=fileparts(name);sample.(stem)=item;
    end
    data.(id)=sample;
end
% 模块二：同一份 HTML 供主窗口和独立文件使用，所有资源均内嵌。
template=fileread(fullfile(root,'web','static-air-probe.html'));
source=fileread(fullfile(root,'web','static-air-probe.js'));
json=strrep(jsonencode(data),'<','\u003c');
page=strrep(template,'<!--PROBE_SCRIPT-->', ['<script>' source '</script>']);
page=strrep(page,'<!--ARCHIVE_JSON-->',json);
end

function text=verifiedText(path,expected)
fid=fopen(path,'rb');assert(fid>=0,'StaticAir:Read','无法读取来源文件 %s',path);
guard=onCleanup(@()fclose(fid));
bytes=fread(fid,inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
assert(strcmp(digest,expected),'StaticAir:SourceChanged','来源文件与归档记录不符：%s',path);
text=native2unicode(bytes','UTF-8');
end
