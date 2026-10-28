function testEntropyLayerEdgeState()
%TESTENTROPYLAYEREDGESTATE 独立验证流管质量匹配根、热力守恒与失败合同。
% 模块作用：制造可解析复算的质量坐标根，不用被测模块产生参考解；
% NASA7 用独立温度未知量的激波求根和熵反算核对。制造输入不能视为
% 自算激波形状、边界层耦合、真实鼻部匹配或实验热流验证。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
fs=makeState(220,5000,4);
cfg=struct('heating',struct('airThermalModel','caloricallyPerfect'));
table=makeTable(fs,'caloricallyPerfect');
target=[.004 .013;.028 .055];
angle=interp1(table.massCoordinate,table.shockCrossingAngleDeg,target,'linear');
reference=caloricReference(fs,angle,fs.pressure);
station=makeStation(fs.pressure,.2,target./(reference.rho.*reference.U*.2));
originalStation=station; originalTable=table;

%% 1. 独立完全气体制造根：质量方程、穿越角、能量和熵分别核对
edge=recoverEntropyLayerEdgeState(fs,station,table,cfg);
assert(edge.available && all(edge.algebraicallyValid(:)));
assert(isequaln(station,originalStation) && isequaln(table,originalTable));
assert(isequal(size(edge.speed),size(target)));
assert(max(abs(edge.massCoordinate(:)-target(:)))<1e-9);
assert(max(abs(edge.shockCrossingAngleDeg(:)-angle(:)))<1e-6);
assert(max(abs(edge.temperature(:)./reference.T(:)-1))<1e-8);
assert(max(abs(edge.speed(:)./reference.U(:)-1))<1e-8);
assert(max(abs(edge.totalPressure(:)./reference.p0(:)-1))<1e-8);
matched=edge.density.*edge.speed.*edge.streamlineMetricInput.*edge.massThicknessInput;
assert(max(abs(edge.massCoordinate(:)-matched(:)))<1e-9);
assert(max(abs(edge.massResidualNormalized(:)))<=edge.relativeMassTolerance);
assert(all(edge.rootCountDetected(:)==1) && all(edge.numericalUniquenessChecked(:)));
assert(edge.scanSampleCount==1+(numel(table.massCoordinate)-1)*16);
assert(~edge.geometryVerified && ~edge.defaultFlowIntegrated);
assert(contains(edge.uniquenessMethod,'非数学唯一性证明'));
assert(strcmp(edge.massCoordinateUnit,'kg/s'));
assert(strcmp(edge.fieldMetadata.pressure.valueKind,'input'));
assert(strcmp(edge.fieldMetadata.speed.valueKind,'computed'));
assert(strcmp(edge.inputMetadata.massThickness.valueKind,'input'));
assert(strcmp(edge.inputMetadata.streamlineMetric.unit,'m'));
assert(contains(edge.entropySource,table.shockShapeSource));
assert(contains(edge.totalEnthalpySource,table.upstreamStateSource));
assert(all(edge.finalBracketLower(:)<=edge.massCoordinate(:)) && ...
    all(edge.finalBracketUpper(:)>=edge.massCoordinate(:)));
cp=fs.gamma*fs.R/(fs.gamma-1);
assert(max(abs(cp*edge.temperature(:)+edge.speed(:).^2/2-edge.upstreamTotalEnthalpy))<1e-5);
entropy=cp*log(edge.temperature/fs.temperature)-fs.R*log(edge.pressure/fs.pressure);
assert(max(abs(entropy(:)-edge.entropyIncrease(:)))<1e-6);

%% 2. 度量协变、数组行列/标量广播、来源逐站保留
scaledTable=table; scaledTable.massCoordinate=table.massCoordinate*7;
scaled=station; scaled.streamlineMetric=station.streamlineMetric*7;
covariant=recoverEntropyLayerEdgeState(fs,scaled,scaledTable,cfg);
assert(all(covariant.algebraicallyValid(:)));
assert(max(abs(covariant.massCoordinate(:)/7-edge.massCoordinate(:)))<1e-10);
assert(max(abs(covariant.speed(:)-edge.speed(:)))<1e-8);
column=makeStation(fs.pressure,.2,station.massThickness(:));
ec=recoverEntropyLayerEdgeState(fs,column,table,cfg);
assert(isequal(size(ec.temperature),[4 1]) && all(ec.algebraicallyValid));
row=makeStation(fs.pressure,.2,station.massThickness(:)');
er=recoverEntropyLayerEdgeState(fs,row,table,cfg);
assert(isequal(size(er.temperature),[1 4]) && all(er.algebraicallyValid));
assert(max(abs(er.speed'-ec.speed))<1e-8);
partial=row; partial.physicalPressureSource={'制造压力甲','','制造压力丙','制造压力丁'};
ep=recoverEntropyLayerEdgeState(fs,partial,table,cfg);
assert(ep.statusCode(2)==2 && isnan(ep.speed(2)) && all(ep.algebraicallyValid([1 3 4])));
assert(strcmp(ep.fieldMetadata.pressure.source{3},'制造压力丙'));

%% 3. 表端点、零h/零质量厚度/驻点极限；没有表内零坐标则拒绝
endTarget=[table.massCoordinate(1) table.massCoordinate(end)];
endRef=caloricReference(fs,table.shockCrossingAngleDeg([1 end]),fs.pressure);
ends=makeStation(fs.pressure,.2,endTarget./(.2*endRef.rho.*endRef.U));
ee=recoverEntropyLayerEdgeState(fs,ends,table,cfg);
assert(all(ee.algebraicallyValid) && isequal(ee.massCoordinate,endTarget));
zeroMetric=makeStation(fs.pressure,0,.01);
ez=recoverEntropyLayerEdgeState(fs,zeroMetric,table,cfg);
assert(ez.algebraicallyValid && ez.massCoordinate==0 && ez.zeroMassLimit);
assert(ez.speed>0 && ez.massResidual==0);
zeroThickness=makeStation(fs.pressure,.2,0);
et=recoverEntropyLayerEdgeState(fs,zeroThickness,table,cfg);
assert(et.algebraicallyValid && et.massCoordinate==0 && et.zeroMassLimit);
stagReference=caloricReference(fs,table.shockCrossingAngleDeg(1),fs.pressure);
stag=makeStation(stagReference.p0,.2,0);
es=recoverEntropyLayerEdgeState(fs,stag,table,cfg);
assert(es.algebraicallyValid && es.speed==0 && es.stagnationLimit && es.massCoordinate==0);
assert(es.temperature==es.upstreamTotalTemperature);
noZeroTable=table; noZeroTable.massCoordinate=table.massCoordinate+.01;
noZero=recoverEntropyLayerEdgeState(fs,zeroMetric,noZeroTable,cfg);
assert(noZero.statusCode==6); assertMissing(noZero);

%% 4. NASA7：独立温度根法复算激波、总熵与最终质量根
thermal=struct; tableT=makeTable(fs,'thermallyPerfect');
psiTarget=.028; gammaTarget=interp1(tableT.massCoordinate,tableT.shockCrossingAngleDeg,psiTarget);
independent=thermalReference(fs,gammaTarget,fs.pressure);
stationT=makeStation(fs.pressure,.2,psiTarget/(.2*independent.rho*independent.U));
thermalEdge=recoverEntropyLayerEdgeState(fs,stationT,tableT,thermal);
assert(thermalEdge.algebraicallyValid);
assert(abs(thermalEdge.massCoordinate-psiTarget)<1e-9);
assert(abs(thermalEdge.temperature-independent.T)<1e-5);
assert(abs(thermalEdge.speed-independent.U)<1e-5);
assert(abs(thermalEdge.totalPressure/independent.p0-1)<1e-8);
[~,h1,s1]=airThermalProperties(fs.temperature);
[~,he,se]=airThermalProperties(thermalEdge.temperature);
assert(abs(he+thermalEdge.speed^2/2-h1-fs.speed^2/2)<1e-5);
assert(abs(se-s1-fs.R*log(thermalEdge.pressure/fs.pressure)-thermalEdge.entropyIncrease)<1e-6);
assert(abs(thermalEdge.massCoordinate-thermalEdge.density*thermalEdge.speed*.2*stationT.massThickness)<1e-10);
fprintf('NASA7 manufactured root: Psi %.12g kg/s, Gamma %.8f deg, T %.8f K, U %.8f m/s.\n', ...
    thermalEdge.massCoordinate,thermalEdge.shockCrossingAngleDeg,thermalEdge.temperature,thermalEdge.speed);

%% 5. 多根必须拒绝；候选根包括表段内部，不能只检查表节点
multiTable=table;
multiTable.massCoordinate=[0 .07 .1 .14];
multiTable.shockCrossingAngleDeg=[90 90 asind(1/fs.Mach) asind(1/fs.Mach)];
multi=makeStation(fs.pressure,1,.0012);
em=recoverEntropyLayerEdgeState(fs,multi,multiTable,cfg);
assert(em.statusCode==7 && em.rootCountDetected==3);
assert(~em.numericalUniquenessChecked); assertMissing(em);
fprintf('Multiple-root refusal: %d candidates detected at %d internal scan samples.\n', ...
    em.rootCountDetected,em.scanSampleCount);

%% 6. 表外、非单调、物性/压力/来源/分离和迭代失败保持NaN真实状态
outside=makeStation(fs.pressure,.2,10);
eo=recoverEntropyLayerEdgeState(fs,outside,table,cfg);
assert(eo.statusCode==6); assertMissing(eo);
bad=makeStation(-1,.2,.01); e=recoverEntropyLayerEdgeState(fs,bad,table,cfg);
assert(e.statusCode==1); assertMissing(e);
bad=makeStation(fs.pressure,.2,-.01); e=recoverEntropyLayerEdgeState(fs,bad,table,cfg);
assert(e.statusCode==1); assertMissing(e);
bad=makeStation(fs.pressure,NaN,.01); e=recoverEntropyLayerEdgeState(fs,bad,table,cfg);
assert(e.statusCode==1); assertMissing(e);
bad=makeStation(fs.pressure,.2,.01); bad=rmfield(bad,'massThicknessSource');
e=recoverEntropyLayerEdgeState(fs,bad,table,cfg); assert(e.statusCode==2); assertMissing(e);
bad=makeStation(fs.pressure,.2,.01); bad.separated=true;
e=recoverEntropyLayerEdgeState(fs,bad,table,cfg); assert(e.statusCode==5); assertMissing(e);
badTables={table,table,table,table,table};
badTables{1}.massCoordinate=[0 .04 .02 .06];
badTables{2}.shockCrossingAngleDeg=[70 60 65 40];
badTables{3}.massCoordinateSource='';
badTables{4}.upstreamState.speed=fs.speed*1.01;
badTables{5}.thermalPropertyMode='thermallyPerfect';
for j=1:numel(badTables)
    e=recoverEntropyLayerEdgeState(fs,column,badTables{j},cfg);
    assert(all(e.statusCode==3)); assertMissing(e);
end
bad=makeStation(1e99,.2,.01); e=recoverEntropyLayerEdgeState(fs,bad,table,cfg);
assert(e.statusCode==4); assertMissing(e);
bad=makeStation(1e-20,.2,.01); e=recoverEntropyLayerEdgeState(fs,bad,tableT,thermal);
assert(e.statusCode==4); assertMissing(e);
tooHot=makeState(220,5000,10); hotTable=makeTable(tooHot,'thermallyPerfect');
e=recoverEntropyLayerEdgeState(tooHot,stationT,hotTable,thermal);
assert(e.statusCode==4); assertMissing(e);
few=cfg; few.entropyLayer=struct('maxIterations',1);
one=makeStation(fs.pressure,.2,station.massThickness(3));
e=recoverEntropyLayerEdgeState(fs,one,table,few);
assert(e.statusCode==8); assertMissing(e);
bad=row; bad.streamlineMetric=ones(4,1);
expectError(@()recoverEntropyLayerEdgeState(fs,bad,table,cfg),'RocketAeroHeat:EntropyLayer:Size');
fprintf('Entropy-layer mass matching passed: independent caloric/NASA7 roots, conservation, metric covariance, endpoints and zero limits, shape/source contracts, multiple-root and domain refusal.\n');
fprintf('Standalone numerical checks only; manufactured shock/mass inputs do not validate generated geometry, boundary-layer coupling, full-area heating or experimental accuracy.\n');
end

function fs=makeState(T,p,M)
gamma=1.4; R=287.05;
fs=struct('Mach',M,'gamma',gamma,'R',R,'temperature',T,'pressure',p, ...
    'density',p/(R*T),'speed',M*sqrt(gamma*R*T));
end

function table=makeTable(fs,mode)
table=struct('massCoordinate',[0 .02 .04 .06], ...
    'shockCrossingAngleDeg',[70 60 50 40], ...
    'massCoordinateSource','独立制造质量坐标，非自算激波流管', ...
    'shockShapeSource','独立制造穿越角表，非自算波形', ...
    'upstreamState',fs,'upstreamStateSource','明确指定的单工况上游解析输入', ...
    'thermalPropertyMode',mode);
end

function station=makeStation(p,h,delta)
station=struct('physicalPressure',p,'streamlineMetric',h,'massThickness',delta,'separated',false, ...
    'physicalPressureSource','制造物理静压，未求表面压力分布', ...
    'streamlineMetricSource','制造同一流管的横向度量h，非Dijkstra路径', ...
    'massThicknessSource','制造物理质量厚度，不用Thwaites系数换算', ...
    'separationSource','制造附着流状态，非真实分离判据');
end

function ref=caloricReference(fs,Gamma,p)
% 完全气体法向跃迁闭式公式；以完整速度的总焓恢复指定物理压力。
g=fs.gamma; cp=g*fs.R/(g-1); Mn=fs.Mach*sind(Gamma);
r=(g+1)*Mn.^2./((g-1)*Mn.^2+2);
p2=fs.pressure*(1+2*g/(g+1)*(Mn.^2-1));
T2=fs.temperature*(p2/fs.pressure)./r;
T0=fs.temperature+fs.speed^2/(2*cp);
p0=p2.*(T0./T2).^(g/(g-1));
T=T0*(p./p0).^((g-1)/g); U=sqrt(2*cp*(T0-T));
ref=struct('p0',p0,'T',T,'U',U,'rho',p./(fs.R*T));
end

function ref=thermalReference(fs,Gamma,p)
% 独立以温度为未知的法向守恒根，排除上游平凡根；不复用被测激波核。
Un=fs.speed*sind(Gamma); [~,h1]=airThermalProperties(fs.temperature);
normalH0=h1+Un^2/2; normalT0=airRecoveryTemperature(fs.temperature,Un,1);
T2=fzero(@(T)shockMomentum(T,fs,Un,normalH0), ...
    [fs.temperature*(1+1e-4),normalT0*(1-1e-6)]);
[~,h2,s2]=airThermalProperties(T2);
Un2=sqrt(2*(normalH0-h2)); rho2=fs.density*Un/Un2; p2=rho2*fs.R*T2;
T0=airRecoveryTemperature(fs.temperature,fs.speed,1); [~,h0,s0]=airThermalProperties(T0);
p0=p2*exp((s0-s2)/fs.R);
T=fzero(@(value)entropyEquation(value,s0,fs.R*log(p0/p)),[80 T0]);
[~,he]=airThermalProperties(T); U=sqrt(2*(h0-he));
ref=struct('p0',p0,'T',T,'U',U,'rho',p/(fs.R*T));
end

function value=shockMomentum(T,fs,Un,h0)
[~,h]=airThermalProperties(T); Un2=sqrt(2*(h0-h));
rho=fs.density*Un/Un2;
value=rho*fs.R*T+rho*Un2^2-fs.pressure-fs.density*Un^2;
end

function value=entropyEquation(T,s0,target)
[~,~,s]=airThermalProperties(T); value=s0-s-target;
end

function assertMissing(edge)
assert(~edge.available && ~any(edge.algebraicallyValid(:)) && ~any(edge.massMatchingSolved(:)));
fields=fieldnames(edge.fieldMetadata);
for j=1:numel(fields)
    name=fields{j};
    % 未收敛时保留已检查的数值括区作为诊断，不当作求解场量。
    diagnostics={'scanBracketLower','scanBracketUpper','finalBracketLower','finalBracketUpper'};
    if ~ismember(name,diagnostics), assert(all(isnan(edge.(name)(:)))); end
    assert(isequal(edge.fieldMetadata.(name).valid,edge.algebraicallyValid));
    assert(isequal(edge.fieldMetadata.(name).statusCode,edge.statusCode));
end
end

function expectError(action,id)
failed=false;
try, action(); catch exception, failed=strcmp(exception.identifier,id); end
assert(failed,'结构错误应返回明确标识');
end
