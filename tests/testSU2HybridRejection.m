function report=testSU2HybridRejection(out)
%TESTSU2HYBRIDREJECTION 复用实际网格检查篡改/换面拒绝和原投影与旧裁剪路径一致。
% 本测试不调用Gmsh或SU2；仅在独立小目录复制回执及小报告，原科学文件不改写。
assert(~isfolder(out));mkdir(out);folder='outputs/v113_validation/actual02';
d=load('outputs/v113_validation/actual_input.mat','V','F');s=struct('V',d.V,'F',d.F);
volume=readSU2HybridVolume(folder,s);bad=s;bad.F(1,[1 2])=bad.F(1,[2 1]);
rejected=false;try,readSU2HybridVolume(folder,bad);catch,rejected=true;end;assert(rejected);
% 使用回执复制与原MAT硬链接不修改原对象；篡改的小报告必须被SHA阻止。
copyfile(fullfile(folder,'mesh_process.json'),fullfile(out,'mesh_process.json'));
copyfile(fullfile(folder,'mesh_report.json'),fullfile(out,'mesh_report.json'));
java.nio.file.Files.createLink(java.io.File(fullfile(out,'hybrid_volume.mat')).toPath(), ...
    java.io.File(fullfile(folder,'hybrid_volume.mat')).toPath());
fid=fopen(fullfile(out,'mesh_report.json'),'a');fprintf(fid,' ');fclose(fid);
rejected=false;try,readSU2HybridVolume(out,s);catch ex,rejected=contains(ex.message,'审核记录SHA');end;assert(rejected);
% 非平面四面体表面的准确三角权重与既有裁剪实现交叉核对，不造物理场。
v=[0 0 0;1 0 0;0 1 0;0 0 1];f=[1 3 2;1 2 4;2 3 4;3 1 4];
small=struct('V',v,'F',f,'param',struct('faceBandID',(1:4)'));
sample=struct('vertices',v,'sourceVertices',v,'sourceFaces',f,'bodyWall',{mat2cell(f,ones(4,1),3)}, ...
    'wallSourceBand',(1:4)','wallSourceFaceID',(1:4)');
a=buildSU2WallProjection(sample,small);sample=rmfield(sample,'wallSourceFaceID');b=buildSU2WallProjection(sample,small);
assert(norm(a.weights-b.weights,inf)<2e-14 && norm(a.sourceNodeAreaVector-b.sourceNodeAreaVector,inf)<2e-14);
report=struct('actualMeshReadWithReportHash',true,'changedOriginalFaceRejected',true, ...
    'tamperedReportRejected',true,'exactTriangleVersusClippingPassed',true,'solverInvoked',false,'mesherInvoked',false);
fid=fopen(fullfile(out,'checks.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
