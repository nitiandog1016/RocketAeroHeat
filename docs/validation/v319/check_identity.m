function check_identity(out)
% 模块一：检查米制文本、显式旋转、整体平移及每面独立顶点，逐面比对独立真值。
[ref,~]=thermistorWorkbenchModel('leads');
shuffled=readmatrix(fullfile(fileparts(out),'fixture','reference_ids.txt'));
f=ref.F(shuffled,[2 3 1]);n=size(f,1);
mesh=struct('F',reshape(1:3*n,3,[])','V',ref.V(reshape(f',[],1),:));
offset=[.002 -.001 .003];mesh.V=mesh.V+offset;
m=mapThermistorSTL(mesh,ref);
assert(isequal(m.referenceFaceID,shuffled) && norm(m.translationM-offset)<1e-15);
r=[1 0 0;0 0 -1;0 1 0];sourceV=mesh.V*r;
path=fullfile(out,'translated-m.stl');fid=fopen(path,'w');assert(fid>=0);
fprintf(fid,'solid probe\n');
fmt=['facet normal 0 0 0\nouter loop\n' repmat('vertex %.17g %.17g %.17g\n',1,3) 'endloop\nendfacet\n'];
fprintf(fid,fmt,sourceV');fprintf(fid,'endsolid probe\n');fclose(fid);
cfg=struct('geometry',struct('stlFile',path,'scale',1,'rotationMatrix',r,'rotationCenter',[0 0 0]));
raw=loadSTLMesh(cfg,true);raw.V=raw.V*r';
source=struct('file',path,'unit','m','sha256',cm011FileSHA256(path),'scale',1);
a=bindThermistorSTL(raw,ref,source,cfg.geometry);assert(isequal(a.referenceFaceID,shuffled));
% 模块二：故意构造错误几何及来源，逐项检查拒绝原因，保留失败面号。
bad=mesh;bad.F(249,:)=bad.F(249,[1 3 2]);
proof.reversed=reject(@()mapThermistorSTL(bad,ref),'RAH:Bead:NormalArea','249');
bad=mesh;bad.F(46,:)=bad.F(45,:);
proof.duplicate=reject(@()mapThermistorSTL(bad,ref),'RAH:Bead:Duplicate','46');
bad=mesh;bad.F(95,:)=[];
proof.missing=reject(@()mapThermistorSTL(bad,ref),'RAH:Bead:Missing',num2str(shuffled(95)));
bad=mesh;bad.V(bad.F(19,1),3)=bad.V(bad.F(19,1),3)+3e-8;
proof.shape=reject(@()mapThermistorSTL(bad,ref),'RAH:Bead:Vertex','19');
bad=mesh;bad.V=mesh.V*1000;
proof.units=reject(@()mapThermistorSTL(bad,ref),'RAH:Bead:Vertex','');
bad=raw;bad.F([2 3],:)=bad.F([3 2],:);
proof.fileOrder=reject(@()bindThermistorSTL(bad,ref,source,cfg.geometry),'RAH:Bead:FileFace','2');
fid=fopen(path,'a');fprintf(fid,'\n');fclose(fid);
proof.fileChanged=reject(@()bindThermistorSTL(raw,ref,source,cfg.geometry),'RAH:Bead:Source','');
proof.maxTextErrorM=a.maxVertexErrorM;proof.textFaces=n;proof.splitVertices=size(mesh.V,1);
fid=fopen(fullfile(out,'identity.json'),'w','n','UTF-8');fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function msg=reject(call,identifier,detail)
try
    call();
catch ex
    assert(strcmp(ex.identifier,identifier) && contains(ex.message,detail),'%s %s',ex.identifier,ex.message);
    msg=ex.message;return;
end
error('负例未被拒绝：%s',identifier);
end
