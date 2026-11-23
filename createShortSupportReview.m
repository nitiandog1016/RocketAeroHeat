function view=createShortSupportReview(host,g)
% 模块一：在科研工作台显示带支管名义外形，接头、测站及面积共用同一套面号。
tabs=uitabgroup(host,'Tag','short-sphere-review','Visible','off','Units','normalized','Position',[0 0 1 1]);
try
    tab=uitab(tabs,'Title','双支管银球');layout=uigridlayout(tab,[4 1]);
    layout.RowHeight={32,30,'1x',36};layout.Padding=[12 8 12 8];
    s=g.summary;
    uilabel(layout,'Text',sprintf('球径 %.2f mm   |   支管直径 %.4f mm   |   两端名义截断位置 Y = ±%.2f mm', ...
        1e3*s.diameterM,1e3*s.tubeDiameterM,1e3*s.endYM),'FontWeight','bold');
    items=[{'两个接头'};cellstr(compose('测站   %g°',g.stations.angleDeg))];
    choice=uidropdown(layout,'Items',items,'ItemsData',0:height(g.stations),'Value',0);
    scene=createAerothermalGeometryView(layout);info=uilabel(layout,'WordWrap','on');
    choice.ValueChangedFcn=@(~,~)selectStation();
    selectStation();
    % 模块二：外露银球、管外壁与人工截断分行记账，已遮部分回查原球编号。
    audit=uitab(tabs,'Title','外露面积与原面号');box=uigridlayout(audit,[5 1]);
    box.RowHeight={32,148,32,'1x',32};box.Padding=[12 8 12 8];
    uilabel(box,'Text',sprintf('原银球 %.4f mm²   |   接头遮盖 %.4f mm²   |   银球外露 %.4f mm²', ...
        1e6*s.smoothGridM2,1e6*s.coveredM2,1e6*s.regionAreaM2(1)),'FontWeight','bold');
    group=table(["外露银球";"+Y 支管外壁";"−Y 支管外壁";"+Y 人工截断";"−Y 人工截断"], ...
        s.regionFaceCount,1e6*s.regionAreaM2,'VariableNames',{'region','count','area'});
    uitable(box,'Data',group,'ColumnName',{'几何分区','三角面数','面积 / mm²'},'RowName',[]);
    affected=g.parents(g.parents.coveredM2>0,:);
    uilabel(box,'Text',sprintf('接头涉及 %d 个原球面，其中 %d 个完全遮盖',height(affected),sum(affected.exposedM2==0)));
    ledger=uitable(box,'Data',[affected.parentID,1e6*affected{:,2:4}], ...
        'ColumnName',{'原球面号','原面积 / mm²','外露 / mm²','遮盖 / mm²'},'RowName',[]);
    uibutton(box,'Text','导出几何、逐面面积与七站定位','ButtonPushedFcn',@saveGeometry);
    tabs.SelectedTab=tab;
    view=struct('tabs',tabs,'geometryTab',tab,'areaTab',audit,'scene',scene,'choice',choice,'info',info, ...
        'ledger',ledger,'geometry',g,'exportTo',@(path)exportShortSupportedSphere(g,path));
catch ex
    delete(tabs);rethrow(ex);
end
    function selectStation()
        % 模块三：用原生面号标出所选测站，并以真实三维线段标出两侧接头。
        html=scene.html;state=html.UserData;payload=state.payload;
        payload.vertices=g.V;payload.faces=g.F-1;payload.caption='Short Test 65 · 双支管名义装置';
        payload.coordinateNote='SI · +X 迎风点 · ±Y 支管 · +Z 测线';
        payload.resetCamera=~isfield(payload,'initialCameraAzimuthRad');
        payload.initialCameraAzimuthRad=.56;payload.initialCameraElevationRad=.47;
        edges=[g.joinPositive;g.joinNegative];
        payload.sectionSegments=reshape(g.V(edges(:),:),size(edges,1),2,3);
        k=choice.Value;payload.inspectedFace=-1;
        if k==0
            info.Text='轮廓线标出球面与管外壁的两处接缝。支管外伸为名义值；双支管工况的壁压与热流待求。';
        else
            row=g.stations(k,:);payload.inspectedFace=row.faceID-1;
            info.Text=sprintf('测站 %g°   |   面号 %d   |   原球面号 %d   |   位置 [%.4f, %.4f, %.4f] mm', ...
                row.angleDeg,row.faceID,row.parentID,1e3*[row.xM row.yM row.zM]);
        end
        state.revision=state.revision+1;payload.revision=state.revision;
        state.payload=payload;state.rendered=[];state.error=[];html.UserData=state;html.Data=payload;
    end
    function saveGeometry(~,~)
        [name,folder]=uiputfile('*.zip','导出双支管银球几何','Short65-supported.zip');
        if isequal(name,0),return;end
        exportShortSupportedSphere(g,fullfile(folder,name));
    end
end
