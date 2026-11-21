function mesh=prepareThermistorWorkbenchMesh(oriented,angularCells,windowOpen)
%PREPARETHERMISTORWORKBENCHMESH 准备所选扇区或全周体网，保留原面及真实壁节点预览。
if nargin<3,windowOpen=[];end
assert(isscalar(angularCells) && ismember(angularCells,[2 4 48]),'RAH:Bead:Angle','请选择已支持的周向网格。');
[canonical,source]=thermistorWorkbenchModel();
assert(isequal(canonical.V,oriented.V) && isequal(canonical.F,oriented.F), ...
    'RAH:Bead:Geometry','当前几何与登记的 B 型珠体不符。');
root=fileparts(mfilename('fullpath'));
work=fullfile(root,'outputs',['thermistor_grid_' char(java.util.UUID.randomUUID())]);mkdir(work);
if angularCells==48
    % 全周网可从几何重建；耗时的首次准备期间继续处理工作台事件。
    manifest=fullfile(work,'prepared.json');
    pulse=[];
    if ~isempty(windowOpen),pulse=@refreshWindow;end
    runCM011Python({fullfile(root,'scripts','runThermistorStartupJob.py'),'prepare',manifest},fullfile(work,'prepare.log'),pulse,true);
    mesh=canonical;p=jsondecode(fileread(manifest));
    p.manifest=manifest;p.manifestSHA256=cm011FileSHA256(manifest);mesh.thermistorPrepared=p;
    verifyThermistorWorkbenchMesh(mesh,angularCells);return;
end
folder=fullfile(work,'mesh');
runCM011Python({fullfile(root,'scripts','runThermistorCase.py'),folder,'--mesh-only','--angular-cells',num2str(angularCells)}, ...
    fullfile(work,'prepare.log'));
manifest=fullfile(folder,'prepared.json');g=jsondecode(fileread(manifest));
mesh=canonical;
mesh.thermistorPrepared=struct('directory',folder,'manifestSHA256',cm011FileSHA256(manifest), ...
    'meshSHA256',g.meshSHA256,'bodySHA256',source.bodySHA256,'angularCells',g.angularCells, ...
    'nodes',g.nodes,'cells',g.hexahedra+g.prisms,'wallFaces',g.boundaryFaces.wall,'preview',g.preview);
verifyThermistorWorkbenchMesh(mesh,angularCells);
    function refreshWindow()
        drawnow limitrate;
        if ~windowOpen(),error('RAH:Bead:Closed','工作台已关闭，本次网格准备仍在独立目录完成。');end
    end
end
