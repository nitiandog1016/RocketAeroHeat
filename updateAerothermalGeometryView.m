function updateAerothermalGeometryView(view,mesh,caption,options)
%UPDATEAEROTHERMALGEOMETRYVIEW 将实际 V/F 同步到无物理场的几何网页。
% 模块一：核验输入和选项；仅传送几何，不制造压力、热流或温度字段。
% options.showWire：是否显示实际三角形边，默认 false。
% options.resetCamera：更新时适配模型，默认 true；false 保留用户视角。
if ~isstruct(view)||~isfield(view,'html')||~isvalid(view.html),return;end
if nargin<2,mesh=[];end
if nargin<3||isempty(caption),caption='几何预览';end
if nargin<4||isempty(options),options=struct;end
if ~isstruct(options)||~isscalar(options)
    error('RocketAeroHeat:GeometryView:Options','几何显示选项必须是标量结构体。');
end
showWire=readFlag(options,'showWire',false);
resetCamera=readFlag(options,'resetCamera',true);
enableFaceProbe=readFlag(options,'enableFaceProbe',false);
if ~(ischar(caption)&&isrow(caption))&&~(isstring(caption)&&isscalar(caption))
    error('RocketAeroHeat:GeometryView:Caption','几何标题必须是单行文本。');
end
V=zeros(0,3);F=zeros(0,3);
if ~isempty(mesh)
    if ~isstruct(mesh)||~isscalar(mesh)||~isfield(mesh,'V')||~isfield(mesh,'F')
        error('RocketAeroHeat:GeometryView:Mesh','预览网格必须包含 V 和 F。');
    end
    V=mesh.V;F=mesh.F;
    if ~isnumeric(V)||~isreal(V)||size(V,2)~=3||any(~isfinite(V(:)))|| ...
            ~isnumeric(F)||~isreal(F)||size(F,2)~=3||any(~isfinite(F(:)))|| ...
            any(F(:)~=fix(F(:)))||any(F(:)<1)||any(F(:)>size(V,1))|| ...
            xor(isempty(V),isempty(F))
        error('RocketAeroHeat:GeometryView:Mesh','几何顶点或三角面编号无效。');
    end
    if ~isempty(V)&&max(max(V,[],1)-min(V,[],1))<=0
        error('RocketAeroHeat:GeometryView:Extent','几何的实际尺寸为零。');
    end
end
% 模块二：每次更新都建立新版本并清除旧绘制回执。坐标统一为 SI。
state=view.html.UserData;revision=state.revision+1;
payload=struct('kind','aerothermal-geometry','schemaVersion',1, ...
    'vertices',double(V),'faces',double(F)-1,'caption',char(caption), ...
    'coordinateNote','SI · XYZ / m · 主轴 X · 鼻尖 +X · 真实几何比例', ...
    'showWire',showWire,'resetCamera',resetCamera,'revision',revision);
% 作者点云使用独立的点序、坐标和逐点色标。连续面元仍由 V/F 分支处理。
if isfield(options,'referencePoints')
    p=options.referencePoints;
    validateattributes(p.xyz,{'double'},{'real','finite','nonempty','ncols',3});
    validateattributes(p.values,{'double'},{'real','finite','column','numel',size(p.xyz,1)});
    validateattributes(p.range,{'double'},{'real','finite','numel',2});
    assert(isempty(mesh) && p.range(2)>p.range(1) && ~isfield(options,'givenWallK'), ...
        'RAH:Reference:Display','参考点云需要独立坐标及有效色标。');
    payload.vertices=p.xyz;
    payload.coordinateNote='作者机体坐标 / m · 发布半域 · 喷口人工封盖';
    payload.pointField=rmfield(p,'xyz');payload.pointField.context='published-samples';
    payload.pointField.productionQualified=false;
    payload.initialCameraAzimuthRad=2.15;payload.initialCameraElevationRad=.48;
end
% 纯几何探针沿用现有三角形拾取；此开关只传回面号和交点。
if enableFaceProbe,payload.enableFaceProbe=true;end
% 给定壁温来自确认的边界。整面等温时色标只显示这一温度。
if isfield(options,'givenWallK')
    t=options.givenWallK(:);
    assert(isnumeric(t) && isreal(t) && numel(t)==size(F,1) && ~isempty(t) && all(isfinite(t) & t>0), ...
        '逐面给定壁温必须与当前面数一致，以 K 为单位且严格为正。');
    payload.faceField=struct('context','given-wall','name','给定壁温','unit','K','values',t, ...
        'range',[min(t),max(t)],'productionQualified',false);
end
% 选区来自当前几何的一基面号；网页接收对应的零基索引。
if isfield(options,'selectedFaces') && ~isempty(options.selectedFaces)
    ids=options.selectedFaces;
    assert(isnumeric(ids) && isreal(ids) && isvector(ids) && all(isfinite(ids)) && ...
        all(ids==fix(ids)) && all(ids>=1 & ids<=size(F,1)),'几何选区含无效面号。');
    payload.selectedFaceIDs=num2cell(double(ids(:)')-1);
end
state.payload=payload;state.revision=revision;state.rendered=[];state.error=[];
view.html.UserData=state;view.html.Data=payload;
end

function value=readFlag(options,name,fallback)
% 模块三：逻辑选项只接受单一布尔或 0/1，避免字符串被误解释为真。
value=fallback;if ~isfield(options,name),return;end
raw=options.(name);
if ~(islogical(raw)||isnumeric(raw))||~isreal(raw)||~isscalar(raw)|| ...
        ~isfinite(raw)||~ismember(raw,[0 1])
    error('RocketAeroHeat:GeometryView:OptionValue','选项 %s 必须为 true/false。',name);
end
value=logical(raw);
end
