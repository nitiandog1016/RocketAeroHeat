function mesh=prepareThermistorLeadMesh(oriented,windowOpen)
% 模块一：用选中的带线原面生成对应三维体网，返回实际壁节点组成的扇区预览。
if nargin<2,windowOpen=[];end
[reference,source]=thermistorWorkbenchModel('leads');
assert(isequal(reference.V,oriented.V) && isequal(reference.F,oriented.F),'RAH:Lead:Geometry','带线试件几何或面序不符。');
root=fileparts(mfilename('fullpath'));work=fullfile(root,'outputs',['thermistor_leads_' char(java.util.UUID.randomUUID())]);
mkdir(work);folder=fullfile(work,'prepared');
runCM011Python({fullfile(root,'scripts','thermistorLeadWorkbench.py'),'prepare',folder}, ...
    fullfile(work,'prepare.log'),@pulse,true);
path=fullfile(folder,'prepared.json');p=jsondecode(fileread(path));
assert(strcmp(source.bodySHA256,p.bodySHA256) && strcmp(p.directory,folder),'RAH:Lead:Source','准备记录与试件不符。');
p.manifest=path;p.manifestSHA256=cm011FileSHA256(path);
mesh=oriented;mesh.thermistorPrepared=p;
    function pulse()
        drawnow limitrate;
        if ~isempty(windowOpen) && ~windowOpen(),error('RAH:Lead:Closed','网格已在独立目录准备，工作台已关闭。');end
    end
end
