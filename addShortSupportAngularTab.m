function page=addShortSupportAngularTab(tabs,p)
% 模块一：相同十层网壁法向分布下，对照两个周向分辨率的实际功率、背风换热和热峰。
tab=uitab(tabs,'Title','64 与 128 扇区对照');box=uigridlayout(tab,[5 1]);
box.RowHeight={32,142,30,'1x',46};box.Padding=[12 8 12 8];
uilabel(box,'FontWeight','bold','Text',sprintf('整件放热    64 扇区 %.6f W    128 扇区 %.6f W    变化 %+.2f%%', ...
    p.coarseW,p.fineW,100*(p.fineW/p.coarseW-1)));
r=p.regions;names=["外露银球";"+Y 管壁";"−Y 管壁";"+Y 湿端盖";"−Y 湿端盖"];
data=table(names,[r.coarseW]',[r.fineW]',[r.coarsePeakWm2]',[r.finePeakWm2]', ...
    'VariableNames',{'region','power64','power128','peak64','peak128'});
regions=uitable(box,'Data',data,'ColumnName',{'分区','64 区 (W)','128 区 (W)','64 区峰值 (W/m²)','128 区峰值 (W/m²)'},'RowName',[]);
uilabel(box,'Text',sprintf('最大放热面号    %d → %d      反向面积    %.3f → %.3f mm²', ...
    p.coarsePeak.faceID,p.finePeak.faceID,1e6*p.baseDirection.total.reverse.areaM2,1e6*p.direction.total.reverse.areaM2));
% 模块二：测点采用两套原生弦各自的插值权重，三维探针仍取原面平均量。
st=p.direction.stations;old=p.baseDirection.stations;
result=table([st.angleDeg]',[st.measuredNu]',[old.computedNu]',[st.computedNu]', ...
    direction(old),direction(st),'VariableNames',{'angle','measured','nu64','nu128','flow64','flow128'});
stations=uitable(box,'Data',result,'ColumnName',{'极角 (°)','实测 Nu','64 区 Nu','128 区 Nu','64 区流向','128 区流向'},'RowName',[]);
uilabel(box,'WordWrap','on','Text',sprintf('外露球面折合量热偏差 %+.2f%% → %+.2f%%。测站按位置对照，温界取同一极角函数。', ...
    100*p.coarsePowerError,100*p.finePowerError));
page=struct('tab',tab,'regions',regions,'stations',stations);
end

function labels=direction(stations)
% 顺流按来流在局部切平面上的投影定义，驻点及异号测点单列。
code=[stations.status]';labels=repmat("未判定",size(code));
labels(code>0)="顺向";labels(code<0)="反向";
labels(strcmp({stations.reason},'noDirection'))="方向退化";
end
