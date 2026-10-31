function report=testSU2SSTInterface(out)
%TESTSU2SSTINTERFACE 内部非零壁热平板的真实SST接口任务，不作公开热流验证。
% 模块一：生成二维小网格，显式全湍流假设；固定预算100步/120秒/2线程。
% 模块二：原生真实终态、SST列名/模型来源与k能量扣除审核，未收敛仍标研究。
assert(~isfolder(out));mkdir(out);folder=fullfile(out,'plate');
buildSU2FlatPlateCase(folder,struct('nx',40,'ny',32,'iterations',100,'outputEvery',50));
python='D:/anacondapython/python.exe';root=fileparts(fileparts(mfilename('fullpath')));
[code,text]=system(sprintf('"%s" -X utf8 "%s" --folder "%s" --provenance "controlled full-turbulent interface assumption; not flight measurement"', ...
    python,fullfile(root,'scripts','prepareSU2SSTCase.py'),folder));assert(code==0,text);
receipt=runEmbeddedSU2Case(folder,struct('pythonExecutable',python,'threads',2,'maximumSeconds',120));
data=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
% 二维长度分担与原生节点合成面积一致；仅构造收敛审核所需的真实负荷。
verts=data.vertices;area=zeros(size(verts,1),1);start=1;
for j=1:numel(data.offsets)
    ids=data.connectivity(start:data.offsets(j))+1;assert(numel(ids)==2);
    area(ids)=area(ids)+norm(verts(ids(1),:)-verts(ids(2),:))/2;start=data.offsets(j)+1;
end
r=struct('receipt',receipt,'reconstructedEngineHeatPowerRaw',area'*data.heatFluxRaw);
a=auditSU2ResultConvergence(folder,r);assert(a.checkpointIdentityConfirmed,a.reason);
assert(a.checkpointIteration==99 && numel(a.checkpointResiduals)==4);
model=describeSU2Runtime(receipt);assert(strcmp(model.mode,'sst-2003m') && contains(model.physicsCaption,'SST'));
report=struct('audit',a,'model',model,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
