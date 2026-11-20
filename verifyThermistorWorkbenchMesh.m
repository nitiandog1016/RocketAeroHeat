function verifyThermistorWorkbenchMesh(mesh,cells)
%VERIFYTHERMISTORWORKBENCHMESH 在使用结果前分别核对全珠原面和已生成的求解体网。
assert(isscalar(cells) && ismember(cells,[2 4]) && isfield(mesh,'thermistorPrepared'), ...
    'RAH:Bead:Mesh','请先生成所选热敏珠网格。');
p=mesh.thermistorPrepared;
[canonical,source]=thermistorWorkbenchModel();
assert(isequal(mesh.V,canonical.V) && isequal(mesh.F,canonical.F) && strcmp(p.bodySHA256,source.bodySHA256), ...
    'RAH:Bead:Geometry','原面几何或原编号与名义试件不符。');
root=fileparts(mfilename('fullpath'));
contract=jsondecode(fileread(fullfile(root,'cases','thermistor_b','grids.json')));
if cells==2,fixed=contract.grids.outer;else,fixed=contract.grids.outerAngular4;end
manifest=fullfile(p.directory,'prepared.json');g=jsondecode(fileread(manifest));
assert(cells==p.angularCells && g.angularCells==cells && strcmp(cm011FileSHA256(manifest),p.manifestSHA256) && ...
    strcmp(p.meshSHA256,fixed.mesh_su2) && strcmp(cm011FileSHA256(fullfile(p.directory,'mesh.su2')),fixed.mesh_su2) && ...
    strcmp(cm011FileSHA256(fullfile(p.directory,'nodes.npz')),fixed.nodes_npz), ...
    'RAH:Bead:Mesh','所选体网、节点或预处理回执已改变，请重新生成。');
end
