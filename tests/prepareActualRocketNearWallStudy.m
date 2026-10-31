function report=prepareActualRocketNearWallStudy(out)
%PREPAREACTUALROCKETNEARWALLSTUDY 核验固定外域候选并导出配对算例，不启动引擎。
% 模块一：独立读入两网格，检查体积、内面双向共享、完整原墙面与远场。
% 模块二：统一原精确工况/200步/100间隔；初态迁移另行审核后才能运行。
assert(~isfolder(out));mkdir(out);d=load('outputs/v118_validation/actual01/surface.mat');
flight=d.flight;flight.iterations=200;flight.outputEvery=100;flight.namedCheckpoints=true;
folders={'outputs/v118_validation/actual01/mesh','outputs/v120_validation/mesh32'};
names={'control','candidate'};report=struct;
for j=1:2
    v=readSU2HybridVolume(folders{j},d.s);folder=fullfile(out,names{j});
    manifest=writeSU2RocketCase(v,folder,flight);projection=buildSU2WallProjection(v,d.s);
    save(fullfile(folder,'projection.mat'),'projection','manifest','-v7');
    assert(size(d.s.F,1)==19752 && abs(v.bodyArea-1.9646292565715553)<1e-12);
    report.(names{j})=struct('nodes',size(v.vertices,1),'cells',numel(v.cells),'topology',v.topology, ...
        'areaM2',v.bodyArea,'fluidVolumeM3',v.fluidVolume,'meshConstruction',v.meshConstruction);
end
assert(abs(report.control.fluidVolumeM3-report.candidate.fluidVolumeM3)<1e-9);
report.physicsQualified=false;report.newQualifiedArea=0;report.gmshRun=false;report.engineRun=false;
fid=fopen(fullfile(out,'preparation.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Both full-domain meshes independently checked and cases written; no engine started.\n');
end
