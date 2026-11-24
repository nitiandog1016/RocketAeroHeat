function page=addShortSupportSurfaceTab(tabs,stations,s)
% 模块一：在三维云图旁展示同一稳态的分区功率、热峰、守恒及实验精度。
tab=uitab(tabs,'Title','五区热量与七站');
box=uigridlayout(tab,[5 1]);box.RowHeight={34,148,'1x',58,34};box.Padding=[12 8 12 8];
if isfield(s,'energy')
    farW=s.energy.farW;curveLabel=sprintf('扰动后 %.6f ms',1000*s.elapsedSeconds);
else
    farW=s.thermal.early.scalarFarW;curveLabel='双支管稳态计算';
end
uilabel(box,'FontWeight','bold','Text',sprintf('%d / %d 面  ·  %.6f mm²  ·  壁放热 %.7f W  ·  远场输出 %.7f W', ...
    s.faces,s.faces,1e6*s.areaM2,s.outW,farW));
names=["外露银球";"+Y 支管外壁";"−Y 支管外壁";"+Y 湿端盖";"−Y 湿端盖"];r=s.regions;
data=table(names,[r.faces]',1e6*[r.areaM2]',[r.outW]',-[r.minFluxIn]',[r.peakOutFace]', ...
    'VariableNames',{'region','faces','area','outW','peakOut','face'});
regions=uitable(box,'Data',data,'RowName',[], ...
    'ColumnName',{'分区','原面数','面积 (mm²)','放热 (W)','最大向外热流 (W/m²)','原面号'});
ax=uiaxes(box);hold(ax,'on');
line=plot(ax,stations.angleDeg,stations.computedNu,'-o','Color',[.15 .32 .49],'DisplayName',['壁反力 · ' curveLabel]);
if ismember('gradientNu',stations.Properties.VariableNames)
    plot(ax,stations.angleDeg,stations.gradientNu,'--s','Color',[.20 .44 .35],'DisplayName','法向温度梯度 · 二次');
end
plot(ax,stations.angleDeg,stations.reportedNu,'d','Color',[.66 .32 .12],'DisplayName','Short 表 III');
xlabel(ax,'从迎风点沿 +Z 测线的极角 (°)');ylabel(ax,'当地表温物性的 Nu');
ax.YGrid='on';xlim(ax,[-4 180]);legend(ax,'Location','best');
note=uilabel(box,'WordWrap','on','Text',sprintf(['仅外露银球的平均热流折合整球 %.4f W，试验量热 %.4f W，偏差 %+.2f%%。' ...
    '七站 RMSE %.2f Nu。支管尺寸和温度采用名义假设；全表面联合实验验证面积 0 m²。'], ...
    s.silverEquivalentFullSphereW,s.observedFullSphereW,100*s.equivalentRelativeError,s.stationRMSE));
btn=uibutton(box,'Text','导出七站位置、读数与误差 CSV','ButtonPushedFcn',@saveStations);
page=struct('tab',tab,'table',regions,'note',note,'line',line,'exportButton',btn,'exportTo',@writeStations);
    function saveStations(~,~)
        [name,folder]=uiputfile('*.csv','导出七站','Short65-supported-stations.csv');
        if ~isequal(name,0),writeStations(fullfile(folder,name));end
    end
    function writeStations(path)
        writetable(stations,path,'Encoding','UTF-8');
    end
end
