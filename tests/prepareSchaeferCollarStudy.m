function report=prepareSchaeferCollarStudy(out)
%PREPARESCHAEFERCOLLARSTUDY 保持原壁面和首24层，仅延长近壁棱柱区，准备同初态配对。
% 模块一：读取已完成2999步及原混合网格，生成唯一32层候选并审核完整拓扑。
% 模块二：双方使用同一物性/绝热边界；新节点只在旧真实四面体内插值守恒初猜。
% 模块三：仅准备新目录，不启动600步数值比较，不改变旧数据或放宽门槛。
assert(~isfolder(out));mkdir(out);[surface,truth]=buildSchaeferFergusonBody();
oldFolder='outputs/hybrid_mesh/a75f5d87-b861-4478-b972-3095018283e6';
oldEngine='outputs/embedded_workbench/96284313-167b-4afc-a198-39d8a5cc2385';
old=readSU2HybridVolume(oldFolder,surface);input=load(old.meshReceipt.input,'optionsJSON');opt=jsondecode(input.optionsJSON);
opt.layerHeights=[0 cumsum(2e-6*1.25.^(0:31))];opt.pythonExecutable='D:/anacondapython/python.exe';opt.maximumSeconds=120;
newFolder=fullfile(out,'mesh32');candidate=buildSU2HybridVolume(surface,newFolder,opt);
nv=size(surface.V,1);assert(isequal(old.vertices(1:25*nv,:),candidate.vertices(1:25*nv,:)));
saved=load(fullfile(oldEngine,'input_snapshot.mat'),'manifest');flight=saved.manifest.flight;
flight.iterations=600;flight.outputEvery=600;flight.namedCheckpoints=false;
items={old,candidate};names={'control','candidate'};
for j=1:2
    folder=fullfile(out,names{j});manifest=writeSU2RocketCase(items{j},folder,flight);
    text=fileread(fullfile(folder,'case.cfg'));text=strrep(text,'RESTART_SOL= NO','RESTART_SOL= YES');
    text=[text newline 'READ_BINARY_RESTART= NO' newline 'SOLUTION_FILENAME= initial' newline];
    writeText(fullfile(folder,'case.cfg'),text);
    projection=buildSU2WallProjection(items{j},surface);save(fullfile(folder,'projection.mat'),'projection','manifest','-v7');
end
source=char(java.io.File(fullfile(oldEngine,'restart_flow_002999.csv')).getCanonicalPath());
copyfile(source,fullfile(out,'control','initial.csv')); % 最终去重另行核对，绝不原位改写初始场。
root=fileparts(fileparts(mfilename('fullpath')));python='D:/anacondapython/python.exe';
cmd=sprintf('"%s" "%s" --source-mesh "%s" --target-mesh "%s" --restart "%s" --output "%s"', ...
    python,fullfile(root,'scripts','transferSU2HybridState.py'),fullfile(root,oldFolder,'hybrid_volume.mat'), ...
    fullfile(root,newFolder,'hybrid_volume.mat'),source,fullfile(root,out,'candidate','initial.csv'));
[code,message]=system(cmd);assert(code==0,'初场迁移未通过：%s',message);
report=struct('oldMeshFolder',oldFolder,'newMeshFolder',newFolder,'oldEngineFolder',oldEngine, ...
    'originalFaces',size(surface.F,1),'oldCells',numel(old.cells),'candidateCells',numel(candidate.cells), ...
    'oldThickness',sum(2e-6*1.25.^(0:23)),'newThickness',opt.layerHeights(end), ...
    'first24LayersExact',true,'oldMeshUnchanged',true,'truth',truth,'physicsQualified',false,'newQualifiedArea',0);
writeText(fullfile(out,'preparation.json'),jsonencode(report));
fprintf('Original wall and 24 layers preserved; one extended-collar candidate prepared, no solver started.\n');
end
function writeText(file,text)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',text);
end
