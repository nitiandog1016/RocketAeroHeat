function candidate=prepareCM011WorkbenchMesh(oriented,source)
%PREPARECM011WORKBENCHMESH 由已确认的 SI 原面生成半域体网；当前面元编号原样保留。
root=fileparts(mfilename('fullpath'));work=fullfile(root,'outputs',['cm011_import_' char(java.util.UUID.randomUUID())]);
assert(mkdir(work),'RAH:CM011:Output','无法创建网格目录。');
V=oriented.V;F=oriented.F;sourceArrays=fullfile(work,'oriented.mat');save(sourceArrays,'V','F','-v7');
folder=fullfile(work,'mesh');
runCM011Python({fullfile(root,'scripts','cm011Imported.py'),'prepare',folder,sourceArrays, ...
    source.file,source.sha256,sprintf('%.17g',source.scale)},fullfile(work,'prepare.log'));
geometry=jsondecode(fileread(fullfile(folder,'geometry.json')));
candidate=oriented;
candidate.cm011Prepared=struct('directory',folder,'manifestSHA256',cm011FileSHA256(fullfile(folder,'prepared.json')), ...
    'sourceSHA256',source.sha256,'nodes',geometry.nodes,'cells',geometry.hexahedra+geometry.prisms, ...
    'maxVertexDistanceM',geometry.importedSurface.maxVertexDistanceM);
end
