function [mesh, stations, report] = runLewisCylinderBenchmark()
%RUNLEWISCYLINDERBENCHMARK 计算低速圆柱的周向换热，并与六条铜带的原始测量对照。
% 模块一：冻结 TM-83557 图9的一个工况，读入可追溯菱形；理论不使用实验列作输入。
root=fileparts(mfilename('fullpath'));
file=fullfile(root,'tests','data','lewis_smooth_cylinder','stations.csv');
stations=readtable(file);
assert(isequal(stations.angleDeg,(0:10:50)') && ...
    all(stations.readLowNu<=stations.figureNu & stations.figureNu<=stations.readHighNu), ...
    'RAH:Lewis:Source','测站的角度或图读范围无效。');
diameter=0.066; span=0.152; stripSpan=0.066; stripWidth=0.0051;
bandDeg=rad2deg(2*stripWidth/diameter);
stations.lowerDeg=stations.angleDeg-bandDeg/2;
stations.upperDeg=stations.angleDeg+bandDeg/2;
stations.pointNu=calcLewisCylinderNu(stations.angleDeg,0);
stations.predictedNu=calcLewisCylinderNu(stations.angleDeg,bandDeg);
stations.relativeError=stations.predictedNu./stations.figureNu-1;
stations.tableFitNu=stations.tableA.*177000.^stations.tableB;
% 模块二：X 为圆柱跨向，+Y 为迎流驻点，+Z 为正角。暴露面只有圆柱侧壁。
% 两平端与试验段侧壁接触；中央 66 mm 是铜带跨度，外侧没有跨向均匀性证据。
angles=(0:359)'; numAngles=numel(angles); axial=[-span -stripSpan stripSpan span]/2;
vertices=[kron(axial(:),ones(numAngles,1)), ...
    repmat((diameter/2)*[cosd(angles) sind(angles)],numel(axial),1)];
faces=zeros(6*numAngles,3); nu=nan(size(faces,1),1);
regions=repmat("靠近风洞侧壁",size(faces,1),1);
for band=1:3
    left=(band-1)*numAngles+(1:numAngles)';
    next=(band-1)*numAngles+mod((1:numAngles)',numAngles)+1;
    first=(band-1)*2*numAngles+(1:numAngles)';
    second=first+numAngles;
    faces(first,:)=[left next next+numAngles];
    faces(second,:)=[left next+numAngles left+numAngles];
    if band==2
        signedAngle=mod(angles+0.5+180,360)-180;
        % 一个等角面带共用该段的平均值；没有将六个观测值扩展到对称面。
        ringNu=calcLewisCylinderNu(signedAngle,1);
        nu([first;second])=[ringNu;ringNu];
        regions([first;second])="中央跨度 · 模型角域外";
        regions([first(isfinite(ringNu));second(isfinite(ringNu))])="中央跨度 · 迎流附着区";
    end
end
mesh=struct('V',vertices,'F',faces);
mesh.localGeometry.region=regions;
mesh.lewisCylinder=struct('nusselt',nu,'reynolds',177000,'diameterM',diameter, ...
    'spanM',span,'stripSpanM',stripSpan,'stripWidthM',stripWidth);
mesh.evidenceIdentity=struct('origin','research','productionQualified',false, ...
    'scope','NASA Lewis 光滑圆柱低湍流工况，中央跨度迎流附着区的带平均换热对照');
% 模块三：完整暴露面积、模型覆盖及实测铜带面积分别登记；绝对压力和热流未辨识。
edges1=vertices(faces(:,2),:)-vertices(faces(:,1),:);
edges2=vertices(faces(:,3),:)-vertices(faces(:,1),:);
area=vecnorm(cross(edges1,edges2,2),2,2)/2;
report=struct('sideAreaM2',sum(area),'analyticSideAreaM2',pi*diameter*span, ...
    'mountingEndAreaM2',2*pi*(diameter/2)^2,'nuAreaM2',sum(area(isfinite(nu))), ...
    'measuredStripAreaM2',height(stations)*stripSpan*stripWidth,'absolutePressureHeatAreaM2',0, ...
    'stationRmsRelative',sqrt(mean(stations.relativeError.^2)), ...
    'maxAbsoluteRelative',max(abs(stations.relativeError)),'bandWidthDeg',bandDeg);
end
