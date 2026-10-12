function testAerothermalViewerHistory(out)
%TESTAEROTHERMALVIEWERHISTORY 用真实热响应和可审计小样本验证历史导出。
% 模块一：无辐射集中热容的独立递推，检验温度、单位和积分时间轴。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
m.V=[0 0 0;2 0 0;0 1 0;0 0 2];m.F=[1 2 3;1 4 2];
m.localGeometry=struct('area',[1;2],'region',{{'sidewall';'aftCap'}});
m.heating=struct('heatTransferCoefficient',[100;200],'recoveryTemperature',[1000;900], ...
    'eligibleFace',[true;false]);
cfg.flight.Twall=300;
mat=struct('density',1000,'specificHeat',1000,'thickness',.01, ...
    'thermalConductivity',10,'emissivity',0,'radiationTemperature',300);
t=[0 2 7];m=calcWallThermalResponse(m,t,mat,cfg);m=calcAerothermalExposure(m);
% 查看器的静态场来自另一个单点工况，不能被时间滑块改写。
m=rmfield(m,'heating');m.pressure=struct('model','单点压力校核','staticPressure',[10000;20000],'Cp',[1;2]);
d=buildAerothermalViewerData(m);f=@(key)d.fields(strcmp({d.fields.key},key));
w=f('wallSolved');h=w.history;
assert(isequal(h.time,t) && strcmp(h.layout,'time-major'));
assert(isequal(h.status,[2;0;1;0;1;0]));
assert(all(isnan(h.values(2:2:end))));
expected=(10000/2*300+100*1000)/(10000/2+100);
assert(abs(h.values(3)-expected)<1e-10);
q=f('wallFinalFlux');assert(abs(q.history.values(3)-100*(1000-expected)/1000)<1e-10);
e=f('cumulativeNetHeat');assert(abs(e.history.values(3)-2*q.history.values(3)/1000)<1e-12);
assert(contains(h.methodNames{h.methodIndex(1)+1},'给定'));
assert(contains(h.methodNames{h.methodIndex(3)+1},'集中热容'));
p=f('pressure');assert(isempty(p.history) && isequal(p.values,[10;20]));
b=f('wallBackFlux');assert(all(b.history.status(1:2:end)==2));
assert(contains(b.history.methodNames{b.history.methodIndex(1)+1},'绝热'));
% 模块二：逐时间方法必须改变；有效掩码优先于有限占位数值。
fake=m;fake.wallThermal.coupling=struct('faceMethodIDHistory',[1 2 1;0 0 0], ...
    'finalHeating',struct('faceMethodNames',{{'无解','层流','湍流'}}), ...
    'lostFace',[false;true],'firstFailureTime',[NaN;2]);
fake.wallThermal.temperature(2,:)=300;
d2=buildAerothermalViewerData(fake);q=d2.fields(strcmp({d2.fields.key},'wallFinalFlux')).history;
assert(contains(q.methodNames{q.methodIndex(1)+1},'层流'));
assert(contains(q.methodNames{q.methodIndex(3)+1},'湍流'));
assert(all(isnan(q.values(2:2:end))) && contains(q.faceNote{2},'首次失效 t=2'));
% 不匹配的累计热量和畸形时间轴应显式拒绝，不能对齐到错误时刻。
bad=m;bad.aerothermalExposure.time=[0 2 8];rejected(bad,'ExposureTime');
bad=m;bad.wallThermal.time=[0 2 2];rejected(bad,'HistoryTime');
bad=m;bad.wallThermal.temperature=300;rejected(bad,'HistorySize');
% 模块三：导出单面元、单时刻和真实局部反馈示例供浏览器复测。
if nargin>0
    exportAerothermalWeb(m,fullfile(out,'history_fixture'));
    exportAerothermalWeb(fake,fullfile(out,'history_methods'));
    one=struct('V',m.V(1:3,:),'F',[1 2 3]);
    one.wallThermal=struct('time',0,'temperature',300,'convectiveHeatFlux',70000, ...
        'validFace',true,'method','单时刻给定 h 校核');
    exportAerothermalWeb(one,fullfile(out,'history_single'));
    one.wallThermal.time=[0 1 3];one.wallThermal.temperature=[300 NaN 350];
    one.wallThermal.convectiveHeatFlux=[70000 NaN 65000];
    one.wallThermal.method='显示边界测试：人工设置中间站缺失，不是物理算例';
    exportAerothermalWeb(one,fullfile(out,'history_gap'));
    gap=buildAerothermalViewerHistory(one);
    assert(isequal(gap(1).history.status,[2;0;1]));
end
fprintf('Viewer history passed: real implicit thermal solution, dose, status, per-time provenance, invalid histories and singleton export.\n');
end
function rejected(m,suffix)
try,buildAerothermalViewerData(m);error('Test:NotRejected','应拒绝无效历史');
catch err,assert(strcmp(err.identifier,['RocketAeroHeat:Viewer:' suffix]),err.message);end
end
