function summary=testWarmWallAttachedMassProfile(outputFolder)
%TESTWARMWALLATTACHEDMASSPROFILE 暖壁有利压梯附着支及原拒绝边界校核。
% 模块作用：用独立SciPy直接绝热BVP数据验证速度过冲、恢复温度、
% 剪应力及收敛，检查冷/暖壁真实热流符号和同剖面物理壁导数。
% 同方程独立数值校核不是实验；不改变球鼻几何/Kn/Re/薄层筛选。
if nargin<1, outputFolder='outputs/v54_validation/nose_branch'; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
folder=fullfile(root,'tests/data/warm_wall');
input=jsondecode(fileread(fullfile(folder,'independent_inputs.json')));
reference=jsondecode(fileread(fullfile(folder,'independent_adiabatic_summary.json')));
fs=struct('R',287.05287,'gamma',1.4);
cfg=struct('heating',struct('airThermalModel','thermallyPerfect','Prandtl',.72));
T0=882.1651199161178; results=cell(numel(input),1);
maxTawError=0; maxShearError=0; maxOvershootError=0; maxPhysicalHeatError=0;

%% 1. 独立直接绝热边界：有过冲仍是正剪应力、无逆流的附着支
for j=1:numel(input)
    p=input(j); b=reference(j).independent;
    [edge,station]=makeStation(p,fs);
    ad=solveAxisymmetricLaminarMassProfile(fs,edge,b.Taw,station,cfg);
    assert(ad.available,ad.status);
    assert(ad.velocityAudit.positiveWallShearSlope && ad.velocityAudit.noBackflow);
    assert(ad.velocityAudit.overshootRatio>0);
    assert(ad.velocityAudit.minimumEtaDerivative<0);
    assert(abs(ad.heatFlux)<2e-3,'独立绝热Taw输入须有很小的原始热流残差');
    maxShearError=max(maxShearError,abs(ad.wallShearStress/b.tau-1));
    maxOvershootError=max(maxOvershootError,abs(ad.velocityAudit.maximumRatio-b.maxF));
    assert(maxShearError<2e-7 && maxOvershootError<2e-8);
    assert(ad.edgeChecks.velocityDefect<=.005+1e-9);
    assert(ad.edgeChecks.relativeTemperatureDefect<=.005+1e-9);
    assert(ad.edgeChecks.relativeDensityDefect<=.005+1e-9);
    assert(ad.odeResidual<ad.numericalOptions.profileResidualTolerance);
    assert(ad.boundaryResidual<1e-7 && ad.thetaPhysical>0 && ad.massThickness>0);
    % 45度及以上真实支必定触发旧单调门限；回归必须覆盖这个原缺陷。
    if p.angle>=45
        assert(ad.velocityAudit.maximumRatio>1+1e-6 || ...
            ad.velocityAudit.minimumEtaDerivative<-1e-6);
    end

    %% 2. 冷壁加热/暖壁冷却；物理坐标焓壁导数独立复算热流
    cold=solveAxisymmetricLaminarMassProfile(fs,edge,300,station,cfg);
    hot=solveAxisymmetricLaminarMassProfile(fs,edge,T0,station,cfg);
    assert(cold.available && hot.available && cold.heatFlux>0 && hot.heatFlux<0);
    for profile={cold,hot}
        r=profile{1}; y=r.profile.physicalY;
        [~,h]=airThermalProperties(r.profile.temperature);
        scale=y(5); c=polyfit(y(1:5)/scale,h(1:5)-h(1),4);
        q=mu(r.wallTemperatureInput)/.72*c(end-1)/scale;
        errorQ=abs(q/r.heatFlux-1); maxPhysicalHeatError=max(maxPhysicalHeatError,errorQ);
        assert(errorQ<2e-4,'同剖面物理壁导数与热流不一致');
    end
    % 实际在两端均有效的物理温度区间求根；不能由输入Taw直接宣称复算。
    taw=fzero(@heatAtWall,[edge.temperature T0],optimset('TolX',1e-7));
    maxTawError=max(maxTawError,abs(taw-b.Taw)); assert(maxTawError<2e-5);
    results{j}=struct('angle',p.angle,'adiabatic',ad,'cold',cold,'hot',hot,'Taw',taw);
    fprintf('Warm-wall beta %.9f / angle %.6f: Taw %.12g K, overshoot %.9g ppm, tau %.12g Pa.\n', ...
        station.beta,p.angle,taw,1e6*ad.velocityAudit.overshootRatio,ad.wallShearStress);
end

%% 3. 数值外界与BVP容差收敛；Pr=1的总焓常量绝热控制
p=input(end); b=reference(end).independent; [edge,station]=makeStation(p,fs);
extended=cfg; extended.massProfile=struct('etaMax',16);
tight=cfg; tight.massProfile=struct('relativeTolerance',1e-8,'absoluteTolerance',1e-10);
r0=results{end}.adiabatic;
r16=solveAxisymmetricLaminarMassProfile(fs,edge,b.Taw,station,extended);
rt=solveAxisymmetricLaminarMassProfile(fs,edge,b.Taw,station,tight);
assert(r16.available && rt.available);
assert(abs(r16.wallShearStress/r0.wallShearStress-1)<1e-6);
assert(abs(rt.wallShearStress/r0.wallShearStress-1)<1e-6);
assert(abs(r16.thetaPhysical/r0.thetaPhysical-1)<1e-6);
assert(abs(r16.velocityAudit.maximumRatio-r0.velocityAudit.maximumRatio)<1e-8);
pr1=cfg; pr1.heating.Prandtl=1;
r1=solveAxisymmetricLaminarMassProfile(fs,edge,T0,station,pr1);
assert(r1.available && abs(r1.heatFlux)<1e-3);
assert(max(abs(r1.profile.states(4,:)-1))<1e-10);
assert(r1.velocityAudit.maximumRatio<=1+1e-9);
% 经典热壁Pr1支也可过冲；不能把Pr1普遍当成不允许过冲的理由。
heated1=solveAxisymmetricLaminarMassProfile(fs,edge,T0+200,station,pr1);
assert(heated1.available && heated1.heatFlux<0 && ...
    heated1.velocityAudit.maximumRatio>1+1e-6);

%% 4. 原物性、动量兼容、剖面残差与有限外缘拒绝仍有效
bad=station; bad.pressureGradient=0;
rejected=solveAxisymmetricLaminarMassProfile(fs,edge,300,bad,cfg);
assert(~rejected.available && rejected.statusCode==2);checkBlank(rejected);
rejected=solveAxisymmetricLaminarMassProfile(fs,edge,70,station,cfg);
assert(~rejected.available && rejected.statusCode==1);checkBlank(rejected);
short=cfg; short.massProfile=struct('etaMax',4,'edgeBuffer',1);
rejected=solveAxisymmetricLaminarMassProfile(fs,edge,b.Taw,station,short);
assert(~rejected.available && rejected.statusCode==5);checkBlank(rejected);
strict=cfg; strict.massProfile=struct('profileResidualTolerance',1e-15);
rejected=solveAxisymmetricLaminarMassProfile(fs,edge,300,station,strict);
assert(~rejected.available && rejected.statusCode==6);checkBlank(rejected);
summary=struct('tested',true,'experimentalValidation',false,'productionModelChanged', ...
    'attached-branch admissibility only; local similarity equations unchanged', ...
    'maxRecoveryTemperatureErrorK',maxTawError,'maxShearRelativeError',maxShearError, ...
    'maxOvershootRatioDifference',maxOvershootError,'maxPhysicalWallHeatRelativeError',maxPhysicalHeatError);
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end
save(fullfile(outputFolder,'warm_kernel_regression.mat'),'summary','results','r16','rt','r1','heated1');
fid=fopen(fullfile(outputFolder,'warm_kernel_regression_summary.json'),'w');
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));clear guard;
fprintf('Warm-wall attached kernel regression passed; independent same-equation numerical verification only.\n');

    function q=heatAtWall(t)
        r=solveAxisymmetricLaminarMassProfile(fs,edge,t,station,cfg);
        assert(r.available,r.status);q=r.heatFlux;
    end
end
function [e,s]=makeStation(p,fs)
e=struct('temperature',p.temperature,'pressure',p.pressure,'density',p.density, ...
    'speed',p.speed,'totalEnthalpy',p.Hphysical,'sources',struct);
for key={'temperature','pressure','density','speed','totalEnthalpy'}
    e.sources.(key{1})='M4/20km球面正激波驻点熵工程外缘；独立剖面数值校核输入';
end
dU=p.beta*p.density*mu(p.temperature)*p.speed^2*p.metric^2/(2*p.xi);
s=struct('streamlineMetric',p.metric,'xi',p.xi,'beta',p.beta, ...
    'pressureGradient',-p.density*p.speed*dU,'velocityGradient',dU,'sources',struct);
for key={'streamlineMetric','xi','beta','pressureGradient','velocityGradient','wallTemperature'}
    s.sources.(key{1})='已保存球面局部站数值基准；给定热边界，不是独立外流或实验验证';
end
assert(abs(e.density*fs.R*e.temperature/e.pressure-1)<1e-12);
end
function value=mu(T)
value=1.458e-6*T.^1.5/(T+110.4);
end
function checkBlank(r)
assert(isnan(r.heatFlux) && isnan(r.massThickness) && isnan(r.thetaPhysical));
assert(~r.fieldMetadata.heatFlux.valid && ~r.fieldMetadata.massThickness.valid);
end
