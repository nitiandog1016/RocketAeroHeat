function updateAerothermalWebTab(ui,mesh)
%UPDATEAEROTHERMALWEBTAB 热响应或工况改变后同步完整结果，防止显示陈旧数据。
if ~isvalid(ui.html),return;end
payload=buildAerothermalViewerData(mesh);
state=ui.html.UserData; state.payload=payload;
ui.html.UserData=state; ui.html.Data=payload;
end
