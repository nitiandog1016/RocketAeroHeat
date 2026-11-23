function page=addShortSupportDirectionTab(tabs,d)
% 模块一：按试验射线的真实测站位置对照两网，区分测点插值与三角面平均。
tab=uitab(tabs,'Title','七站近壁流向');box=uigridlayout(tab,[4 1]);
box.RowHeight={36,'1x',52,52};box.Padding=[12 8 12 8];
uilabel(box,'Text','Short Test 65   |   沿 +Z 测线，自迎风点计角','FontWeight','bold');
c=d.coarse.stations;f=d.fine.stations;
t=table([f.angleDeg]',[f.faceID]',direction(c),direction(f),value(c,'shearPa'),value(f,'shearPa'), ...
    value(c,'speedMps'),value(f,'speedMps'),'VariableNames',{'angle','face','coarse','fine','tau5','tau10','u5','u10'});
ledger=uitable(box,'Data',t,'ColumnName',{'极角 (°)','原面号','五层流向','十层流向','五层剪切 (Pa)', ...
    '十层剪切 (Pa)','五层首层速度 (m/s)','十层首层速度 (m/s)'},'RowName',[]);
uilabel(box,'WordWrap','on','Text',sprintf('双指标同向反转的面积：五层 %.3f mm²，十层 %.3f mm²，两网共同 %.3f mm²。', ...
    1e6*d.coarse.total.reverse.areaM2,1e6*d.fine.total.reverse.areaM2,1e6*d.overlap.bothM2));
uilabel(box,'WordWrap','on','Text',['剪切与速度均为顺流方向的有符号分量。0° 的切向基退化。' ...
    '测站按相邻原生节点插值；三维探针读取所选原面的平均量。两网首层位置分别约 25.4 和 12.7 µm。']);
page=struct('tab',tab,'table',ledger);
end

function result=value(stations,key)
% JSON 的空值保留为 NaN，测站方向未定义时仍保持七行对齐。
result=nan(numel(stations),1);
for i=1:numel(stations)
    v=stations(i).(key);if ~isempty(v),result(i)=v;end
end
end

function text=direction(stations)
% 将双符号诊断和轴点退化转为便于核对的测点说明。
state=[stations.status]';text=repmat("未判定",size(state));
text(state>0)="顺向";text(state<0)="反向";
text(strcmp({stations.reason},'noDirection'))="方向退化";
end
