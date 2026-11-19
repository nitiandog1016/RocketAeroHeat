function candidate=prepareCM011WorkbenchMesh(oriented,source,halfSectors)
%PREPARECM011WORKBENCHMESH 由已确认的 SI 原面生成半域体网；当前面元编号原样保留。
if nargin<3,halfSectors=12;end
assert(isscalar(halfSectors)&&ismember(halfSectors,[12 24]),'RAH:CM011:Mesh','半周单元数须为 12 或 24。');
root=fileparts(mfilename('fullpath'));work=fullfile(root,'outputs',['cm011_import_' char(java.util.UUID.randomUUID())]);
assert(mkdir(work),'RAH:CM011:Output','无法创建网格目录。');
V=oriented.V;F=oriented.F;sourceArrays=fullfile(work,'oriented.mat');save(sourceArrays,'V','F','-v7');
folder=fullfile(work,'mesh');
runCM011Python({fullfile(root,'scripts','cm011Imported.py'),'prepare',folder,sourceArrays, ...
    source.file,source.sha256,sprintf('%.17g',source.scale),'--half-sectors',num2str(halfSectors)},fullfile(work,'prepare.log'));
geometry=jsondecode(fileread(fullfile(folder,'geometry.json')));
assert(numel(geometry.rayRadians)==halfSectors+1,'RAH:CM011:Mesh','实际周向网格与选择不符。');
candidate=oriented;
candidate.cm011Prepared=struct('directory',folder,'manifestSHA256',cm011FileSHA256(fullfile(folder,'prepared.json')), ...
    'sourceSHA256',source.sha256,'halfSectors',halfSectors,'nodes',geometry.nodes,'cells',geometry.hexahedra+geometry.prisms, ...
    'maxVertexDistanceM',geometry.importedSurface.maxVertexDistanceM);
% 新剖分分别显示 CAD 轮廓的三角面抽样偏差和实际求解壁面的投影误差。
if isfield(geometry.importedSurface,'maxFacetSampleDistanceM')
    candidate.cm011Prepared.maxFacetSampleDistanceM=geometry.importedSurface.maxFacetSampleDistanceM;
    candidate.cm011Prepared.maxVertexToSolverWallM=geometry.importedSurface.maxVertexToSolverWallM;
end
end
