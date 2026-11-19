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
    % 单独的给定工况接受用户热边界，实测试验及其原始读数继续保存在相邻页。
    caseTab=uitab(tabs,'Title','给定边界工况');
    panel=uigridlayout(caseTab,[3 1]);panel.RowHeight={40,40,'1x'};
    bar=uigridlayout(panel,[1 7]);bar.Padding=[0 0 0 0];bar.ColumnWidth={'fit','1x','fit','1x','fit','1x',110};
    labels={'上游静压 / Pa','来流静温 / K','圆柱壁温 / K'};defaults=[101325 293.15 303.15];
    inputs=gobjects(1,3);
    for j=1:3
        uilabel(bar,'Text',labels{j});
        inputs(j)=uieditfield(bar,'numeric','Value',defaults(j),'Limits',[eps inf], ...
            'ValueChangedFcn',@invalidate);
    end
    calculate=uibutton(bar,'Text','计算工况','ButtonPushedFcn',@calculateCase);
    feedback=uigridlayout(panel,[1 2]);feedback.Padding=[0 0 0 0];feedback.ColumnWidth={'1x',146};
    status=uilabel(feedback,'Text','输入研究工况的绝对压力和温度后计算。','WordWrap','on');
    trialButton=uibutton(feedback,'Text','载入未收敛三维试算','Enable','off','ButtonPushedFcn',@showTrial);
    host=uipanel(panel,'BorderType','none');
    tabs.SelectedTab=tab;
    view=struct('tabs',tabs,'comparisonTab',tab,'surface',surface,'stations',stations,'report',report, ...
        'caseTab',caseTab,'inputs',inputs,'calculate',calculate,'trialButton',trialButton,'caseHost',host,'caseStatus',status);
catch ex
    delete(tabs);rethrow(ex);
end
    function invalidate(~,~)
        delete(host.Children);host.UserData=[];trialButton.Enable='off';status.Text='边界已改变，请重新计算。';
    end
    function calculateCase(~,~)
        invalidate([],[]);
        try
            [result,~,summary]=runLewisCylinderCondition(inputs(1).Value,inputs(2).Value,inputs(3).Value);
            inner=uitabgroup(host,'Units','normalized','Position',[0 0 1 1]);
            web=createAerothermalWebTab(inner,result);
            host.UserData=struct('web',web,'report',summary,'mesh',result);
            if isequal([inputs.Value],defaults),trialButton.Enable='on';end
            status.Text=summary.wallCoverage.coverageMessage;
        catch ex
            delete(host.Children);host.UserData=[];status.Text=ex.message;
        end
    end
    function showTrial(~,~)
        try
            active=host.UserData;meshWithTrial=attachLewisTunnelTrial(active.mesh);
            updateAerothermalWebTab(active.web,meshWithTrial);
            trialButton.Enable='off';
            status.Text='变量列表已加入「未收敛试算」壁压和热流。工程覆盖及部分功率仍按原工程场计。';
        catch ex
            status.Text=ex.message;
        end
    end
end
