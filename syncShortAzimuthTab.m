function syncShortAzimuthTab(tabs,mesh,evidenceTab,view)
% 模块一：在同一科研工作台维护与当前结果匹配的周向页。更换步长或模型后清除旧页。
tab=findobj(tabs,'Type','uitab','Tag','Short65Azimuth');
if ~isfield(mesh,'shortSupportAzimuth')
    if ~isempty(tab),delete(tab);end
    return;
end
az=mesh.shortSupportAzimuth;seal=[az.archiveSHA256 az.stlSHA256];
if ~isempty(tab) && strcmp(tab.UserData.seal,seal),return;end
if ~isempty(tab),delete(tab);end
tab=uitab(tabs,'Title','七站周向换热','Tag','Short65Azimuth');grid=uigridlayout(tab,[1 1]);grid.Padding=[0 0 0 0];
html=uihtml(grid,'HTMLSource',fullfile(fileparts(mfilename('fullpath')),'web','short65-azimuth.html'));
payload=shortAzimuthPayload(az);tab.UserData=struct('html',html,'seal',seal,'rendered',[]);
html.DataChangedFcn=@receive;html.Data=payload;
    function receive(~,~)
        % 模块二：曲线发来的面号只在本次确有球面读数时允许传入现有三维探针。
        m=html.Data;if ~isstruct(m) || ~isfield(m,'kind'),return;end
        switch m.kind
            case 'azimuth-ready'
                html.Data=payload;
            case 'azimuth-rendered'
                if strcmp(m.frame,payload.frame),tab.UserData.rendered=m;end
            case 'azimuth-locate'
                if strcmp(m.frame,payload.frame) && isscalar(m.face) && any(az.points.fileFaceID==m.face & az.points.status=="sphere")
                    tabs.SelectedTab=evidenceTab;view.selectField('heatFlux');view.inspectFace(m.face,[]);
                end
        end
    end
end
