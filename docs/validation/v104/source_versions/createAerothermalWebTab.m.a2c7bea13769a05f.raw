function ui=createAerothermalWebTab(tabs,mesh)
%CREATEAEROTHERMALWEBTAB 在统一窗口嵌入离线 WebGL 结果页。
% 模块一：建立网页组件；计算仍在 MATLAB，渲染和面元拾取在网页。
tab=uitab(tabs,'Title','三维结果工作台');
grid=uigridlayout(tab,[1 1]); grid.Padding=[0 0 0 0];
html=uihtml(grid,'HTMLSource',fullfile(fileparts(mfilename('fullpath')),'web','index.html'));
payload=buildAerothermalViewerData(mesh);
html.UserData=struct('payload',payload,'selectedFace',[],'ready',false,'rendered',[]);
html.DataChangedFcn=@receive;
html.Data=payload;
ui=struct('tab',tab,'html',html);
    function receive(src,~)
        % 模块二：只处理白名单事件，避免点击消息覆盖完整模型。
        message=src.Data; state=src.UserData;
        if ~isstruct(message) || ~isfield(message,'kind'),return;end
        if strcmp(message.kind,'ready')
            state.ready=true; src.UserData=state; src.Data=state.payload;
        elseif strcmp(message.kind,'rendered')
            state.rendered=message;src.UserData=state;
        elseif strcmp(message.kind,'selection') && isfield(message,'face')
            id=message.face;
            if isnumeric(id) && isscalar(id) && isfinite(id) && id==fix(id) && ...
                    id>=1 && id<=size(state.payload.faces,1)
                state.selectedFace=id; src.UserData=state;
            end
        end
    end
end
