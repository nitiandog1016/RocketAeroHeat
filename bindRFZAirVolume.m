function air = bindRFZAirVolume(ticket,raw,mesh,cfg)
% 模块一：核对落盘资产、标准坐标和当前工况。体网文件中的编号均从 1 开始。
summaryFile=fullfile(ticket.directory,'air.json');
assetFile=fullfile(ticket.directory,'air.mat');
assert(strcmp(cm011FileSHA256(summaryFile),ticket.summarySHA256),'RAH:RFZ:Air','体网回执已变更。');
s=jsondecode(fileread(summaryFile));
input=cfg.rfzSurfaceInput;
assert(strcmp(s.schema,'rfz-cad-air-v1') && strcmp(s.fileSHA256,cm011FileSHA256(assetFile)) && ...
    strcmp(s.stlSHA256,input.stlSHA256) && strcmp(s.cadTableSHA256,input.cadTableSHA256) && ...
    strcmp(s.canonicalMeshSHA256,phaseAwareMeshIdentity(raw)),'RAH:RFZ:Air','体网与原始几何不匹配。');
d=load(assetFile,'vertices','tetrahedra','wall','outer','wallOwner','cellSICN');
assert(isequal(double(d.wall),raw.F,mesh.F) && isequal(d.vertices(1:size(raw.V,1),:),raw.V), ...
    'RAH:RFZ:Air','原面顺序或顶点与已建空气域不符。');
rotation=cfg.geometry.rotationMatrix;
center=cfg.geometry.rotationCenter(:)';
expected=raw.V;
if ~isequal(rotation,eye(3)),expected=(rotation*(expected-center)')'+center;end
assert(max(abs(expected-mesh.V),[],'all')<1e-12,'RAH:RFZ:Air','当前网格坐标与给定刚体朝向不符。');
[rebuilt,wall]=buildRFZSurfaceInput(mesh,cfg.geometry,cfg.rfzCavityMode,cfg.rfzCavityK);
assert(isequal(rebuilt,input) && isequal(cfg.heating.wallTemperatureByFace,wall.givenWallK), ...
    'RAH:RFZ:Air','给定壁温或工况身份已变更。');

% 模块二：原面号逐行索引空气单元；形状质量随单元保存，温度仍取当前确认输入。
owner=d.wallOwner;
assert(isequal(size(owner),[input.faces 1]) && isinteger(owner) && ...
    numel(d.cellSICN)==s.tetrahedra && size(d.tetrahedra,1)==s.tetrahedra && ...
    all(owner>=1 & owner<=s.tetrahedra) && all(isfinite(d.cellSICN) & d.cellSICN>0), ...
    'RAH:RFZ:Air','体单元编号或形状质量无效。');
air=struct('input',input,'owner',owner,'cellSICN',d.cellSICN,'summary',s, ...
    'thermal',buildWallNodeConstraints(d,wall));
end
