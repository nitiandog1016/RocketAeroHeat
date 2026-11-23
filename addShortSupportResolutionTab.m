function page=addShortSupportResolutionTab(tabs,pair)
% 模块一：用相同原面分区比较两套近壁网的功率和热点，并显示实测七站。
tab=uitab(tabs,'Title','近壁网格敏感性');box=uigridlayout(tab,[5 1]);
box.RowHeight={34,150,'1x',52,38};box.Padding=[12 8 12 8];
heading=uilabel(box,'FontWeight','bold','Text',sprintf('五层 %.6f W    十层 %.6f W    整件变化 %+.3f%%', ...
    pair.coarseW,pair.fineW,100*(pair.fineW/pair.coarseW-1)));
r=pair.regions;region=["外露银球";"+Y 管壁";"−Y 管壁";"+Y 湿端盖";"−Y 湿端盖"];
t=table(region,[r.coarseW]',[r.fineW]',100*[r.relativePowerChange]',[r.coarsePeakWm2]',[r.finePeakWm2]', ...
    'VariableNames',{'region','coarseW','fineW','changePct','coarsePeak','finePeak'});
ledger=uitable(box,'Data',t,'ColumnName',{'分区','五层 (W)','十层 (W)','变化 (%)','五层峰值 (W/m²)','十层峰值 (W/m²)'},'RowName',[]);
ax=uiaxes(box);hold(ax,'on');plot(ax,pair.anglesDeg,pair.coarseNu,'--o','Color',[.47 .49 .48],'DisplayName','五层网');
line=plot(ax,pair.anglesDeg,pair.fineNu,'-o','Color',[.14 .34 .50],'DisplayName','十层网');
plot(ax,pair.anglesDeg,pair.measuredNu,'d','Color',[.65 .31 .12],'DisplayName','Short 表 III');
xlabel(ax,'迎风点沿 +Z 测线的极角 (°)');ylabel(ax,'当地表温物性的 Nu');ax.YGrid='on';legend(ax,'Location','best');
note=uilabel(box,'WordWrap','on','Text',sprintf(['原面热流绝对差：面积加权中位数 %.2f W/m²，95%% 分位 %.2f W/m²。' ...
    '外露银球折合量热相对偏差：五层 %+.2f%%，十层 %+.2f%%。'], ...
    pair.heatAbsDifference.p50,pair.heatAbsDifference.p95,100*pair.coarsePowerError,100*pair.finePowerError));
uilabel(box,'WordWrap','on','Text','三维页主值来自十层网，变量菜单可选择五层网或同面差值；两网共用色标。');
page=struct('tab',tab,'table',ledger,'heading',heading,'note',note,'line',line);
end
