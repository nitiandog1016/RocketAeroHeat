function [mesh,info]=inspectSU2MeshCandidate(folder,surface)
%INSPECTSU2MESHCANDIDATE 只读定位已完成但质量未过的网格，不创建可计算体网格。
% 模块一：核对实际失败终态、原表面、输出SHA及完整拓扑审核证据。
% 模块二：独立重算四面体SICN，保留原体编号；绝不伪造成功回执或物理资格。
receipt=jsondecode(fileread(fullfile(folder,'mesh_process.json')));
assert(strcmp(receipt.status,'terminal') && receipt.exitCode==1 && ~receipt.timedOut && receipt.sourceStable && receipt.inputStable, ...
    '仅支持源稳定且已终态的失败候选。');
for name={'hybrid_volume.mat','mesh_report.json','geometry_audit.json'}
    file=fullfile(folder,name{1});key=strrep(name{1},'.','_');
    assert(isfield(receipt.outputSHA256,key) && strcmp(sha(file),receipt.outputSHA256.(key)),'候选文件身份不符。');
end
report=jsondecode(fileread(fullfile(folder,'mesh_report.json')));geometry=jsondecode(fileread(fullfile(folder,'geometry_audit.json')));
assert(geometry.passed && ~report.qualityGatePassed && report.tetsBelowSICN01>0,'本入口仅查看完整拓扑通过而质量失败的网格。');
d=load(fullfile(folder,'hybrid_volume.mat'));assert(isequal(d.sourceVertices,surface.V) && isequal(double(d.sourceFaces),surface.F),'候选不属于当前原表面。');
assert(isequal(d.vertices(1:size(surface.V,1),:),surface.V));
T=double(d.tetrahedra);offset=size(d.wedges,1)+size(d.pyramids,1);
assert(size(T,1)==geometry.tetrahedra && size(d.vertices,1)==report.nodes);
diagnostic=struct('vertices',d.vertices,'cells',{mat2cell(T,ones(size(T,1),1),4)}, ...
    'cellTypes',10*ones(size(T,1),1),'sourceVertices',surface.V,'sourceFaces',surface.F,'meshConstruction',report);
[mesh,info]=buildSU2MeshQualityPreview(diagnostic);
assert(info.lowQualityCells==report.tetsBelowSICN01 && abs(info.minimumSICN-report.minimumTetSICN)<1e-12);
% 诊断临时数组只含四面体，恢复棱柱/金字塔在前的真实体编号。
info.sourceCellIDs=info.sourceCellIDs+offset;mask=info.displaySourceCellIDs>0;info.displaySourceCellIDs(mask)=info.displaySourceCellIDs(mask)+offset;
info.acceptedForCalculation=false;info.candidateExitCode=receipt.exitCode;info.candidateDirectory=folder;
info.exteriorAlgorithm='Delaunay';if isfield(report,'exteriorAlgorithm'),info.exteriorAlgorithm=report.exteriorAlgorithm;end
info.volumeCounts=struct('nodes',report.nodes,'prisms',report.prisms,'pyramids',report.pyramids,'tetrahedra',report.tetrahedra);
info.fullTopologyEvidence=geometry;info.hybridVolumeSHA256=receipt.outputSHA256.hybrid_volume_mat;
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
