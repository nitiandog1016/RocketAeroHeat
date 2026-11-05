function report = testSU2IdealAirInitialExport(output)
%TESTSU2IDEALAIRINITIALEXPORT 从正常MATLAB网格/配置入口准备新物性初场的两个输入。
% 模块一。复用v145已通过接口的小型混合网格及完整真实网格，保持原节点和原壁。
% 模块二。写出新的独立计算目录，工况沿用原记录，计算留给显式后台任务。
assert(~isfolder(output));mkdir(output);
source = jsondecode(fileread('outputs/v145_validation/main01/engine_folder.json'));
original = jsondecode(fileread(fullfile(source.folder,'case_manifest.json')));
folders = {'outputs/v145_validation/near_wall_fixture','outputs/v145_validation/mesh'};
names = {'small','real'};
records = cell(1,2);
for index = 1:2
    data = load(fullfile(folders{index},'mesh_input.mat'));
    surface = struct('V',data.V,'F',data.F);
    volume = readSU2HybridVolume(folders{index},surface);
    flight = original.flight;
    flight.iterations = 2;flight.outputEvery = 1;flight.namedCheckpoints = true;
    if index == 1
        flight.temperature = 226.5;flight.pressure = 1171.8;
    end
    folder = fullfile(output,names{index});
    info = writeSU2RocketCase(volume,folder,flight);
    assert(isequal(volume.sourceVertices,surface.V) && isequal(volume.sourceFaces,surface.F));
    if index == 2
        assert(info.nodes==335726 && info.cells==705490 && info.originalSurfaceFaceCount==19752);
        assert(abs(info.wallArea-original.wallArea)<1e-14);
        assert(strcmp(sha(fullfile(folder,'mesh.su2')),sha(fullfile(source.folder,'mesh.su2'))));
    else
        assert(info.nodes==38 && info.cells==92 && info.originalSurfaceFaceCount==8);
        assert(strcmp(sha(fullfile(folder,'mesh.su2')),sha('outputs/v145_validation/fixture02/case/mesh.su2')));
    end
    records{index} = struct('name',names{index},'nodes',info.nodes,'cells',info.cells, ...
        'originalFaces',info.originalSurfaceFaceCount,'areaM2',info.wallArea,'meshSHA256',sha(fullfile(folder,'mesh.su2')));
end
report = struct('passed',true,'cases',{records},'nativeRuns',0,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(output,'export_report.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end

function digest=sha(file)
% 模块三。逐块核对原网格身份，避免把全部大文件一次读入内存。
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));
hash=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),hash.update(fread(fid,1024*1024,'*uint8'));end
digest=lower(reshape(dec2hex(typecast(hash.digest(),'uint8'),2)',1,[]));
end
