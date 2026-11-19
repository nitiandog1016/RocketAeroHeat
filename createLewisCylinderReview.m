function view=createLewisCylinderReview(parent,mesh,stations,report)
%CREATELEWISCYLINDERREVIEW 在当前工作台共享原面云图、铜带读数和算法来源。
tabs=uitabgroup(parent,'Units','normalized','Position',[0 0 1 1], ...
    'Tag','lewis-cylinder-review','Visible','off');
try
    surface=createAerothermalWebTab(tabs,mesh);
    tab=uitab(tabs,'Title','风洞测站对照');
    grid=uigridlayout(tab,[4 1]);
    grid.RowHeight={36,'1x',190,56};grid.Padding=[16 10 16 10];
    uilabel(grid,'Text','NASA Lewis 光滑圆柱    |    Re_D = 177 000    |    低湍流空气', ...
        'FontSize',16,'FontWeight','bold');
    ax=uiaxes(grid);hold(ax,'on');
    % 实验图读范围与测量误差分列；实线连接的是六个有限宽度铜带的预测。
    angle=stations.angleDeg;
    plot(ax,angle,stations.predictedNu,'-o','Color',[.11 .34 .53],'DisplayName','Frössling 铜带平均');
    errorbar(ax,angle,stations.figureNu,stations.figureNu-stations.readLowNu, ...
        stations.readHighNu-stations.figureNu,'d','LineStyle','none','Color',[.61 .28 .11], ...
        'DisplayName','图9菱形及读图范围');
    xlabel(ax,'距迎流驻点的角度 / °');ylabel(ax,'Nu_D');
    ax.XTick=angle';xlim(ax,[-3 53]);ax.YGrid='on';legend(ax,'Location','best');
    rows=stations{:,{'angleDeg','lowerDeg','upperDeg','figureNu','predictedNu','relativeError','tableFitNu'}};
    rows(:,6)=100*rows(:,6);
    uitable(grid,'Data',round(rows,2),'ColumnName',{'中心 / °','起点 / °','终点 / °', ...
        '实测图读 Nu','带平均预测 Nu','相对差 / %','表 III 回归 Nu'},'RowName',[]);
    uilabel(grid,'WordWrap','on','Text',sprintf(['六站相对 RMS %.2f%%。图读范围反映扫描和定位。' ...
        '原报告给出全部数据平均误差 5.7%%、最大 7.8%%。绝对热流与壁压未知。'],100*report.stationRmsRelative));
    tabs.SelectedTab=tab;
    view=struct('tabs',tabs,'comparisonTab',tab,'surface',surface,'stations',stations,'report',report);
catch ex
    delete(tabs);rethrow(ex);
end
end
