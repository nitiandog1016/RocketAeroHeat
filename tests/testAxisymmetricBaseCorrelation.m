function testAxisymmetricBaseCorrelation()
%TESTAXISYMMETRICBASECORRELATION 底压/代表底热的公式、物理域和源图校核。
% 模块用途：独立复算原报告公式，检查量纲/符号、连续流拒绝、输入
% 与求解量区分；调用内部直锥外流封装。源图比较使用原实验散点，
% 但散点属于关联式建模数据，仅验证公开关联复算，不称独立实验验证。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
out=fullfile(root,'outputs','v53_validation','base_region');
if ~exist(out,'dir'), mkdir(out); end
fs=makeState(180,5000,5);
geo=struct('family','cone','boundaryLayerRegime','laminar', ...
    'baseRadius',0.2,'halfAngleDeg',9,'noseBluntnessRatio',0, ...
    'hasRearSupport',false,'hasGasInjection',false,'flatCircularBase',true, ...
    'isolatedAxisymmetricBody',true,'source','解析尖头直圆锥校核几何');
Tw=300;

%% 1. 内部肩部外流与独立底压/底热公式；正负热流及给定壁温身份
r=solveSharpConeBaseFlow(fs,geo,Tw); assert(r.available && r.pressureAvailable);
e=r.shoulder; cp=fs.gamma*fs.R/(fs.gamma-1);
Re=e.density*e.speed*e.distance/e.dynamicViscosity;
C4=(1+2/(fs.gamma-1)/e.Mach^2)^(-2);
pbExpected=fs.pressure*3.05*(e.Mach^2/sqrt(Re))^.6/C4;
StExpected=1e-4*(pbExpected/fs.pressure*C4)^1.8;
qExpected=StExpected*fs.density*fs.speed*(cp*fs.temperature+fs.speed^2/2-cp*Tw);
assert(abs(r.basePressure/pbExpected-1)<1e-12);
assert(abs(r.baseStantonNumber/StExpected-1)<1e-12);
assert(abs(r.heatFluxRepresentative/qExpected-1)<1e-12 && r.heatFluxRepresentative>0);
assert(r.basePressure<e.pressure && r.internalExternalFlowComputed);
assert(~r.eligibleForSurfaceCoverage && ~r.radialHeatingDistributionSolved);
assert(isnan(r.baseEdgeTemperature) && isnan(r.baseVelocity) && isnan(r.solvedWallTemperature));
assert(strcmp(r.fieldMetadata.wallTemperatureInput.valueKind,'input'));
assert(strcmp(r.fieldMetadata.basePressure.valueKind,'computed'));
assert(strcmp(r.fieldMetadata.boundaryLayerRegime.valueKind,'input'));
hot=solveAxisymmetricBaseCorrelation(fs,e,geo,r.heatRecoveryReferenceTemperature+10);
equal=solveAxisymmetricBaseCorrelation(fs,e,geo,r.heatRecoveryReferenceTemperature);
assert(hot.heatFluxRepresentative<0 && equal.heatFluxRepresentative==0);
assert(equal.heatFluxAvailable && equal.fieldMetadata.heatFluxRepresentative.valid);
unknown=solveAxisymmetricBaseCorrelation(fs,e,geo,NaN);
assert(unknown.pressureAvailable && unknown.stantonAvailable && ~unknown.heatFluxAvailable);
assert(isnan(unknown.heatFluxRepresentative));

%% 2. 湍流 Eq.(11)/Eq.(18) 的整体开平方；圆柱只计算底压
gt=geo; gt.boundaryLayerRegime='turbulent';
t=solveAxisymmetricBaseCorrelation(fs,e,gt,Tw); assert(t.available);
J=1.7/log(21/gt.halfAngleDeg);
ratio=(e.pressure/fs.pressure)^2*(.025+.906/(1+.2*e.Mach^2))^J;
StT=8.3e-4*sqrt(ratio*(fs.pressure/e.pressure)^2)^1.85;
assert(abs(t.basePressure/(ratio*fs.pressure)-1)<1e-12);
assert(abs(t.baseStantonNumber/StT-1)<1e-12);
% 若把 Eq.(18) 的平方根写成只作用于最后的压力比，本检查会失败。
wrong=8.3e-4*(ratio*fs.pressure/e.pressure)^1.85;
assert(abs(t.baseStantonNumber/wrong-1)>.05);
gc=gt; gc.family='cylinder';
c=solveAxisymmetricBaseCorrelation(fs,e,gc,Tw);
assert(c.pressureAvailable && ~c.stantonAvailable && ~c.heatFluxAvailable);
assert(abs(c.basePressure/e.pressure-(.05+.967/(1+.2*e.Mach^2)))<1e-12);
assert(isnan(c.heatFluxRepresentative));

%% 3. 尺度与域门限：不能放宽流态、攻角、支撑或热物性补值
scaledFS=fs; scaledFS.pressure=fs.pressure*4; scaledFS.density=fs.density*4;
scaledE=e; scaledE.pressure=e.pressure*4; scaledE.density=e.density*4;
scaled=solveAxisymmetricBaseCorrelation(scaledFS,scaledE,gt,Tw);
assert(scaled.available && abs(scaled.basePressure/t.basePressure-4)<1e-12);
assert(abs(scaled.heatFluxRepresentative/t.heatFluxRepresentative-4)<1e-12);
cases={};
f=fs; f.alpha=1; cases{end+1}={f,e,geo};
g=geo; g.hasRearSupport=true; cases{end+1}={fs,e,g};
g=geo; g.hasGasInjection=true; cases{end+1}={fs,e,g};
g=geo; g.flatCircularBase=false; cases{end+1}={fs,e,g};
g=geo; g.halfAngleDeg=15; cases{end+1}={fs,e,g};
s=e; s.distance=1e-8; cases{end+1}={fs,s,geo};
s=e; s.speed=s.speed*1.02; cases{end+1}={fs,s,geo};
f=fs; f.gamma=1.3; cases{end+1}={f,e,geo};
g=geo; g.baseRadius=1e-9; cases{end+1}={fs,e,g};
for k=1:numel(cases)
    a=cases{k}; invalid=solveAxisymmetricBaseCorrelation(a{1},a{2},a{3},Tw);
    assert(~invalid.pressureAvailable && ~invalid.heatFluxAvailable);
    assert(isnan(invalid.basePressure) && isnan(invalid.heatFluxRepresentative));
end

%% 4. 源报告实验散点复算：保留所有数字化圆点，不用拟合线当实验
file=fullfile(root,'tests','data','base_correlation','sand93_fig22_23_points.csv');
data=readtable(file,'TextType','string');
lam=data.regime=="laminar"; tur=data.regime=="turbulent";
assert(nnz(lam)==14 && nnz(tur)==13);
assert(all(data.pixel_x>0 & data.pixel_y>0 & data.measured_stanton>0));
correlation=nan(height(data),1);
correlation(lam)=1e-4*data.parameter(lam).^1.8;
correlation(tur)=8.3e-4*data.parameter(tur).^1.85;
assert(max(abs(correlation./data.equation_stanton-1))<1e-12);
rel=correlation./data.measured_stanton-1;
summary=struct('laminarPoints',nnz(lam),'turbulentPoints',nnz(tur), ...
    'laminarRmsRelativeResidual',sqrt(mean(rel(lam).^2)), ...
    'turbulentRmsRelativeResidual',sqrt(mean(rel(tur).^2)), ...
    'laminarMaxAbsoluteRelativeResidual',max(abs(rel(lam))), ...
    'turbulentMaxAbsoluteRelativeResidual',max(abs(rel(tur))), ...
    'independentExperimentValidation',false, ...
    'comparisonKind','公开关联的建模实验散点复算；无原始工况/径向测点重建', ...
    'basePressurePa',r.basePressure,'baseHeatFluxRepresentative',r.heatFluxRepresentative, ...
    'baseStanton',r.baseStantonNumber,'heatSurfaceCoverageEligible',false);
% 宽于像素读取误差的差异属于源数据离散，记录实值，不调模型系数。
assert(summary.laminarRmsRelativeResidual<.3 && summary.turbulentRmsRelativeResidual<.3, ...
    '源图的散点和报告关联无法复算，请检查数字化或公式');
fid=fopen(fullfile(out,'matlab_summary.json'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(summary)); clear cleanup;
save(fullfile(out,'base_correlation_results.mat'),'r','t','c','data','summary');
fprintf('testAxisymmetricBaseCorrelation 通过：源图层流14/湍流13点，RMS %.4f/%.4f。\n', ...
    summary.laminarRmsRelativeResidual,summary.turbulentRmsRelativeResidual);
end

function fs=makeState(T,p,M)
% 模块用途：生成满足同一理想气体/总焓定义的独立物理校核工况。
g=1.4; R=287.05287;
fs=struct('temperature',T,'pressure',p,'Mach',M,'gamma',g,'R',R, ...
    'density',p/(R*T),'speed',M*sqrt(g*R*T),'alpha',0,'beta',0, ...
    'totalTemperature',T*(1+(g-1)/2*M^2), ...
    'totalPressure',p*(1+(g-1)/2*M^2)^(g/(g-1)));
end
