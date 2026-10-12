function testAerothermalFlowHistory(out)
%TESTAEROTHERMALFLOWHISTORY 验证逐站来源、缺失站、独立有效性及自适应接受步。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
m=struct('V',[0 0 0;2 0 0;0 1 0;0 0 2],'F',[1 2 3;1 4 2]);
p=struct('model','第一站工程压力','staticPressure',[10000;20000], ...
    'Cp',[1;2],'isPrescribedBase',[false;true]);
fs=struct('Mach',4,'altitude',20000,'alpha',0,'beta',0,'unitFlowDirection',[-1 0 0]);
p.freestream=fs;m.pressure=p;
s1=captureAerothermalFlowSnapshot(m);
m.pressure.staticPressure=[30000;40000];m.pressure.Cp=[3;4];
m.pressure.model='末站工程压力';m.pressure.freestream.alpha=10;
s2=captureAerothermalFlowSnapshot(m);
m.wallThermal=struct('time',[0 1 2],'validFace',[false;false], ...
    'temperature',nan(2,3),'convectiveHeatFlux',nan(2,3),'method','测试热壁失效');
meta1=struct('time',0,'flight',fs,'flowSnapshot',s1);
meta2=struct('time',2,'flight',m.pressure.freestream,'flowSnapshot',s2);
m.wallThermal.coupling.boundaryMetadata={meta1,[],meta2};
d=buildAerothermalViewerData(m);f=d.fields(strcmp({d.fields.key},'pressure'));
assert(numel(f)==1 && isequaln(f.history.values,[10;20;NaN;NaN;30;40]));
assert(isequal(f.history.status,[1;2;0;0;1;2]));
assert(contains(f.history.methodNames{f.history.methodIndex(1)+1},'第一站'));
assert(contains(f.history.methodNames{f.history.methodIndex(5)+1},'末站'));
assert(contains(f.history.methodNames{f.history.methodIndex(6)+1},'给定'));
assert(contains(f.history.timeNote{3},'攻角 10'));
assert(all(d.fields(strcmp({d.fields.key},'wallSolved')).status==0));
bad=m;bad.wallThermal.coupling.boundaryMetadata{3}.time=1.9;reject(bad,'FlowTime');
bad=m;bad.wallThermal.coupling.boundaryMetadata{3}.flowSnapshot.condition(3)=9;reject(bad,'FlowCondition');
bad=m;bad.wallThermal.coupling.boundaryMetadata{3}.flowSnapshot.faceCount=1;reject(bad,'FlowSize');
bad=m;bad.flightHeatingHistory.couplingInput=struct('vertices',m.V,'faces',flipud(m.F));reject(bad,'FlowMesh');
legacy=m;legacy.wallThermal=rmfield(legacy.wallThermal,'coupling');
old=buildAerothermalViewerData(legacy);assert(isempty(old.fields(strcmp({old.fields.key},'pressure')).history));

%% 2. 没有壁温求解的原始飞行计算站也能浏览；工况不符必须拒绝
raw=rmfield(m,'wallThermal');
raw.flightHeatingHistory=struct('time',[0 2],'flowStations',{{s1,s2}}, ...
    'flight',struct('Mach',[4 4],'altitude',[20000 20000],'alpha',[0 10],'beta',[0 0]));
data=buildAerothermalViewerData(raw);p=data.fields(strcmp({data.fields.key},'pressure'));
assert(isequal(p.history.values,[10;20;30;40]));

%% 3. 人工刚性热边界强制自动二分，外流快照必须只来自实际接受时刻
one=struct('V',m.V(1:3,:),'F',[1 2 3]);
mat=struct('density',1000,'specificHeat',1000,'thickness',.01, ...
    'thermalConductivity',.5,'emissivity',.8,'radiationTemperature',300,'thermalNodes',6);
one.wallThermal=solveCoupledWallConduction(@boundary,boundary(300,0),[0 .1],mat,300);
assert(one.wallThermal.coupling.timeSubdivisions>0);
data=buildAerothermalViewerData(one);p=data.fields(strcmp({data.fields.key},'pressure'));
assert(isequal(p.history.time,one.wallThermal.time));
assert(max(abs(p.history.values(:)'-(10+one.wallThermal.time)))<1e-12);
if nargin>0
    exportAerothermalWeb(m,fullfile(out,'flow_gap'));
    exportAerothermalWeb(raw,fullfile(out,'flow_only'));
    exportAerothermalWeb(one,fullfile(out,'flow_adaptive'));
end
fprintf('Flow history passed: independent validity, prescribed pressure, station provenance, missing/stale records, flow-only view and adaptive accepted times.\n');
    function h=boundary(T,t)
        snapshot=s1;snapshot.faceCount=1;
        for j=1:numel(snapshot.fields)
            f=snapshot.fields(j);f.values=10+t;f.status=1;f.methodIndex=0;
            f.methodNames={'人工时间标识字段；不代表物理压力求解'};snapshot.fields(j)=f;
        end
        h=struct('heatTransferCoefficient',(100+1000/(1+exp((T-320)/2)))*(1+.1*t), ...
            'recoveryTemperature',1000+100*t,'eligibleFace',true, ...
            'boundaryMetadata',struct('time',t,'flight',fs,'flowSnapshot',snapshot));
    end
end
function reject(m,suffix)
try,buildAerothermalViewerData(m);error('Test:MissingError','应拒绝无效历史');
catch err,assert(strcmp(err.identifier,['RocketAeroHeat:Viewer:' suffix]),err.message);end
end
