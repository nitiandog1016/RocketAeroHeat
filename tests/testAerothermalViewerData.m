function testAerothermalViewerData(out)
%TESTAEROTHERMALVIEWERDATA 检查面积覆盖、无解保护及不同温度的来源。
% 两个不等面积面元避免“按数量统计”错误碰巧通过。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
m.V=[0 0 0;2 0 0;0 1 0;0 0 2];m.F=[1 2 3;1 4 2];
m.localGeometry.region={'sidewall';'aftCap'};
m.pressure=struct('model','测试压力方法','staticPressure',[10000;20000], ...
    'Cp',[1;2],'isPrescribedBase',[false;true],'isPrescribedLeeward',[false;false]);
m.boundaryLayerEdge=struct('pressure',[10000;50000], ...
    'temperature',[900;999],'algebraicallyValid',[true;false], ...
    'faceMethodID',[1;0],'faceMethodNames',{{'无解','切楔'}});
m.heating=struct('heatFlux',[200000;12345],'eligibleFace',[true;false], ...
    'faceMethodID',[1;0],'faceMethodNames',{{'无解','层流'}}, ...
    'validity','测试适用性','recoveryTemperature',[1000;NaN], ...
    'wallTemperatureByFace',[300;300]);
m.heatingCoverage=struct('reasonCode',[0;1],'reasonNames',{{'已计算','尾盖无模型'}});
m.wallThermal=struct('time',[0 1],'temperature',[300 350;300 300], ...
    'validFace',[true;false],'method','测试厚度导热', ...
    'convectiveHeatFlux',[200000 180000;0 0]);
data=buildAerothermalViewerData(m);
f=@(key)data.fields(strcmp({data.fields.key},key));
p=f('pressure');assert(isequal(p.status,[1;2]));
assert(abs(p.computedAreaFraction-1/3)<1e-12);
assert(abs(p.prescribedAreaFraction-2/3)<1e-12);
q=f('heatFlux');assert(q.values(1)==200 && isnan(q.values(2)));
assert(q.status(2)==0 && abs(q.missingAreaFraction-2/3)<1e-12);
e=f('edgeTemperature');assert(isnan(e.values(2)));
w=f('wallSolved');assert(w.values(1)==350 && isnan(w.values(2)));
t=f('wallInput');assert(isequal(t.values,[300;300]) && all(t.status==2));
assert(isequal(data.faces,m.F-1));
assert(isequal(data.missingReason,{'已计算';'尾盖无模型'}));
for k=1:numel(data.fields)
    f0=data.fields(k);
    assert(abs(f0.computedAreaFraction+f0.prescribedAreaFraction+f0.missingAreaFraction-1)<1e-12);
end
% JSON 中未解值必须是 null，不能作为零送到网页。
encoded=jsonencode(data,'ConvertInfAndNaN',true);assert(contains(encoded,'null'));
if nargin>0
    exportAerothermalWeb(m,fullfile(out,'fixture'));
    % 单面元 JSON 会将数值列压成标量，网页必须支持此类最小数据集。
    one=struct('V',m.V(1:3,:),'F',[1 2 3]);
    one.pressure=struct('model','单面元测试','staticPressure',10000,'Cp',1);
    exportAerothermalWeb(one,fullfile(out,'single_face'));
end
fprintf('Viewer payload passed: unequal-area coverage, prescribed boundaries, invalid finite inputs, wall provenance.\n');
end
