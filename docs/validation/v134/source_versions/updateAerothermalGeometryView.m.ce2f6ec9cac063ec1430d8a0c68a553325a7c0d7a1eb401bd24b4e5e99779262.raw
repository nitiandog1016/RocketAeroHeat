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
