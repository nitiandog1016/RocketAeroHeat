function mesh=prepareThermistorLeadMesh(oriented,windowOpen,variant,refineLeads,stlSource,geometry)
% 模块一：为选定的带线或光滑中段加密试件生成体网，并返回真实壁节点预览。
if nargin<2,windowOpen=[];end
if nargin<3,variant='leads';end
if nargin<4,refineLeads=false;end
variant=validatestring(variant,{'leads','midbody','sidewind'});options={};script='thermistorLeadWorkbench.py';
if strcmp(variant,'midbody'),options={'--midbody'};end
model=variant;
if strcmp(variant,'sidewind')
    script='thermistorSidewindWorkbench.py';model='leads';
    if isequal(refineLeads,4),script='thermistorMeasuredMesh.py';end
elseif refineLeads
    assert(strcmp(variant,'leads'));[~,options]=thermistorLeadGridChoice(refineLeads);
end
[reference,source]=thermistorWorkbenchModel(model);
if nargin>=5
    % 普通 STL 经逐面核对后从登记几何重建体网，保留实际导入壁面。
    assert(strcmp(variant,'sidewind') && isequal(refineLeads,4),'RAH:Bead:Grid','导入带线 B 珠支持 96 档侧吹网。');
    oriented.thermistorImport=bindThermistorSTL(oriented,reference,stlSource,geometry);
    options={'--rebuild'};
else
    assert(isequal(reference.V,oriented.V) && isequal(reference.F,oriented.F),'RAH:Lead:Geometry','试件几何或面序不符。');
end
root=fileparts(mfilename('fullpath'));work=fullfile(root,'outputs',['thermistor_' variant '_' char(java.util.UUID.randomUUID())]);
mkdir(work);folder=fullfile(work,'prepared');
runCM011Python({fullfile(root,'scripts',script),'prepare',folder,options{:}}, ...
    fullfile(work,'prepare.log'),@pulse,true);
path=fullfile(folder,'prepared.json');p=jsondecode(fileread(path));
assert(strcmp(source.bodySHA256,p.bodySHA256) && strcmp(p.directory,folder),'RAH:Lead:Source','准备记录与试件不符。');
if nargin>=5,assert(~p.reusedGrid);end
p.manifest=path;p.manifestSHA256=cm011FileSHA256(path);
mesh=oriented;mesh.thermistorPrepared=p;
    function pulse()
        drawnow limitrate;
        if ~isempty(windowOpen) && ~windowOpen(),error('RAH:Lead:Closed','网格已在独立目录准备，工作台已关闭。');end
    end
end
