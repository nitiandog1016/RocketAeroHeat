function mesh=prepareThermistorWorkbenchMesh(oriented,angularCells)
%PREPARETHERMISTORWORKBENCHMESH 生成独立扇区体网，保留全珠原面和原生壁面预览。
assert(isscalar(angularCells) && ismember(angularCells,[2 4]),'RAH:Bead:Angle','请选择已支持的周向网格。');
[canonical,source]=thermistorWorkbenchModel();
assert(isequal(canonical.V,oriented.V) && isequal(canonical.F,oriented.F), ...
    'RAH:Bead:Geometry','当前几何与登记的 B 型珠体不符。');
root=fileparts(mfilename('fullpath'));
work=fullfile(root,'outputs',['thermistor_grid_' char(java.util.UUID.randomUUID())]);mkdir(work);
folder=fullfile(work,'mesh');
runCM011Python({fullfile(root,'scripts','runThermistorCase.py'),folder,'--mesh-only','--angular-cells',num2str(angularCells)}, ...
    fullfile(work,'prepare.log'));
manifest=fullfile(folder,'prepared.json');g=jsondecode(fileread(manifest));
mesh=canonical;
mesh.thermistorPrepared=struct('directory',folder,'manifestSHA256',cm011FileSHA256(manifest), ...
    'meshSHA256',g.meshSHA256,'bodySHA256',source.bodySHA256,'angularCells',g.angularCells, ...
    'nodes',g.nodes,'cells',g.hexahedra+g.prisms,'wallFaces',g.boundaryFaces.wall,'preview',g.preview);
verifyThermistorWorkbenchMesh(mesh,angularCells);
end
