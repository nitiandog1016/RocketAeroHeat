function check_file_count(out)
% 用同一珠线接头的两个登记三角形单独检查文件面数，六个输入面全部是子面。
[body,~]=thermistorWorkbenchModel('leads');faces=body.F([183 184],:);
[vi,~,fi]=unique(faces(:));ref=struct('V',body.V(vi,:),'F',reshape(fi,2,3));
stl=ref;stl.F=zeros(6,3);
for j=1:2
    f=ref.F(j,:);new=size(stl.V,1)+1;stl.V(new,:)=mean(ref.V(f,:),1);
    stl.F(3*j-2:3*j,:)=[f' f([2 3 1])' repmat(new,3,1)];
end
path=fullfile(out,'two-parent-patch.stl');stlwrite(triangulation(stl.F,stl.V),path,'text');
g=struct('stlFile',path,'scale',1,'rotationMatrix',eye(3),'rotationCenter',[0 0 0]);
m=loadSTLMesh(struct('geometry',g),true);
source=struct('unit','m','scale',1,'file',path,'sha256',cm011FileSHA256(path));
map=bindThermistorSTL(m,ref,source,g);assert(isequal(map.referenceFaceID,[1;1;1;2;2;2]));
% 在完整的二进制样件上再核实一次文件面号、面积及父面。
p=fullfile(fileparts(out),'fixture');g.scale=.001;g.stlFile=fullfile(p,'probe-partition-mm.stl');
m=loadSTLMesh(struct('geometry',g),true);source.scale=g.scale;source.file=g.stlFile;source.sha256=cm011FileSHA256(source.file);
b=bindThermistorSTL(m,body,source,g);truth=readmatrix(fullfile(p,'face_truth.csv'));
assert(isequal(b.referenceFaceID,truth(:,1)));
fid=fopen(fullfile(out,'file.json'),'w');fprintf(fid,'%s',jsonencode(struct('patchFaces',numel(map.childIDs),'fullFaces',size(truth,1),'source',source.sha256)));fclose(fid);
end
