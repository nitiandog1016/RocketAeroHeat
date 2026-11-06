function updateAerothermalWebTab(ui,mesh)
%UPDATEAEROTHERMALWEBTAB 热响应或工况改变后同步完整结果，防止显示陈旧数据。
if ~isvalid(ui.html),return;end
payload=buildAerothermalViewerData(mesh);
state=ui.html.UserData; state.payload=payload;
ui.html.UserData=state; ui.html.Data=payload;
% 新结果同时替换区域证据页，避免云图已更新而证据仍指向旧工况。
if isfield(ui,'evidenceView'),ui.evidenceView.update(payload);end
end
