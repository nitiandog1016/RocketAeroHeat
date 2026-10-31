function view=createAerothermalGeometryView(host,mesh)
%CREATEAEROTHERMALGEOMETRYVIEW 在统一工作台创建无物理场的离线几何预览。
% 模块一：复用调用者的宿主；几何导入和网格检查不打开新的图形窗口。
% host 可为 uigridlayout、uipanel、uitab 或 uifigure。网格布局中直接嵌入
% uihtml，其他宿主中创建无边距单元布局。mesh 可省略或为空。
if nargin<2,mesh=[];end
if ~isvalid(host),error('RocketAeroHeat:GeometryView:Host','几何预览宿主无效。');end
if isa(host,'matlab.ui.container.GridLayout')
    parent=host;
else
    parent=uigridlayout(host,[1 1]);parent.Padding=[0 0 0 0];
    parent.RowSpacing=0;parent.ColumnSpacing=0;
end
html=uihtml(parent,'HTMLSource',fullfile(fileparts(mfilename('fullpath')), ...
    'web','geometry.html'));
% 模块二：完整负载与网页回执分别保存；不把回执当成新的几何负载。
html.UserData=struct('payload',[],'ready',false,'rendered',[], ...
    'error',[],'revision',0);
html.DataChangedFcn=@receive;
view=struct('html',html,'host',host);
updateAerothermalGeometryView(view,mesh,'几何预览');

    function receive(src,~)
        message=src.Data;state=src.UserData;
        if ~isstruct(message)||~isfield(message,'kind'),return;end
        if strcmp(message.kind,'ready')
            state.ready=true;src.UserData=state;src.Data=state.payload;
        elseif strcmp(message.kind,'geometry-rendered')
            % 模块三：仅当前版本回执有效，防止旧网页帧覆盖新网格状态。
            if isfield(message,'revision')&&message.revision==state.revision
                state.rendered=message;state.error=[];src.UserData=state;
            end
        elseif strcmp(message.kind,'geometry-error')
            if isfield(message,'revision')&&message.revision==state.revision
                state.error=message;src.UserData=state;
            end
        end
    end
end
