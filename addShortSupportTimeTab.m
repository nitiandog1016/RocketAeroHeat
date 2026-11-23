function page=addShortSupportTimeTab(tabs,p)
% 模块一：三次实际物理步各取自己的前一时间层，分别给出壁放热、远场输出和储热率。
tab=uitab(tabs,'Title','同刻步长与储热');box=uigridlayout(tab,[5 1]);
box.RowHeight={32,106,32,'1x',52};box.Padding=[12 8 12 8];
uilabel(box,'FontWeight','bold','Text',sprintf('Short Test 65    扰动后 %.6f ms    BDF1    Re = %.0f',1e3*p.elapsedS,p.condition.re));
es=[p.fullEnergy,p.firstHalfEnergy,p.halfEnergy];
data=table(["整步";"第一半步";"第二半步"],1000*[es.elapsedS]',[es.wallW]',[es.farW]',1e6*[es.storageW]',100*[es.relativeBalance]', ...
    'VariableNames',{'step','time','wall','far','storage','balance'});
ledger=uitable(box,'Data',data,'ColumnName',{'物理步','终时刻 (ms)','壁放热 (W)','远场输出 (W)','空气储热率 (µW)','相对净差 (%)'},'RowName',[]);
% 模块二：测站以原生弦插值，面积 RMS 采用双支管原面分母，局部与整体误差分别显示。
uilabel(box,'Text',sprintf('反向面积    整步 %.4f mm²    两半步 %.4f mm²    热峰面号 %d / %d', ...
    1e6*p.fullReverse.areaM2,1e6*p.halfReverse.areaM2,p.fullPeak.faceID,p.halfPeak.faceID));
st=table(p.angles(:),p.experimentNu(:),p.fullNu(:),p.halfNu(:), ...
    'VariableNames',{'angle','measured','full','half'});
stations=uitable(box,'Data',st,'ColumnName',{'测站极角 (°)','Short 表 III · Nu','整步 · Nu','两半步 · Nu'},'RowName',[]);
a=p.stepComparison;
uilabel(box,'WordWrap','on','Text',sprintf('同刻面积 RMS 差：壁压 %.4g Pa，入壁热流 %.4g W/m²。七站实验供参照，测量与数值扰动没有同步时间。',a.pressureRMSPa,a.heatRMSWm2));
page=struct('tab',tab,'ledger',ledger,'stations',stations);
end
