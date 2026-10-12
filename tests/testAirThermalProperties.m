function testAirThermalProperties()
% 干空气变比热：核对基准量、焓守恒反解、壁温反馈和适用域。
addpath(fileparts(fileparts(mfilename('fullpath'))));

%% 1. NASA7 比热与焓导数；行列向量都应保持尺寸
[cp,h]=airThermalProperties([300 1000 1500]);
assert(isequal(size(cp),[1 3]));
assert(abs(cp(1)-1005)<5 && cp(2)>1100 && cp(3)>cp(2));
assert(all(diff(h)>0));
step=0.01;
[~,above]=airThermalProperties(1200+step);
[~,below]=airThermalProperties(1200-step);
localCp=airThermalProperties(1200);
assert(abs((above-below)/(2*step)/localCp-1)<1e-6);
[~,~,entropyPlus]=airThermalProperties(1200+step);
[~,~,entropyMinus]=airThermalProperties(1200-step);
assert(abs((entropyPlus-entropyMinus)/(2*step)*1200/localCp-1)<1e-6);
[cpColumn,~]=airThermalProperties([300;1000]);
assert(isequal(size(cpColumn),[2 1]));

%% 2. 80~150 K 低温冻结段的焓、熵与 150 K 连续
[coldCp,coldH,coldS]=airThermalProperties([80 100 149.999 150]);
assert(all(isfinite([coldCp coldH coldS])));
assert(max(abs(coldCp(1:3)-coldCp(4)))<1e-10);
assert(abs(coldCp(1)/(3.5*287.05)-1)<0.01);
assert(abs(coldH(3)-coldH(4))<1.1);
assert(abs(coldS(3)-coldS(4))<0.02);
assert(abs(coldH(2)-coldH(4)-coldCp(4)*(100-150))<1e-8);
assert(abs(coldS(2)-coldS(4)-coldCp(4)*log(100/150))<1e-8);

%% 3. 恢复温度满足焓方程；超出物性范围保持无解
Te=[300;600]; U=[1500;900]; factor=0.85;
Taw=airRecoveryTemperature(Te,U,factor);
[~,edgeH]=airThermalProperties(Te);
[~,recoverH]=airThermalProperties(Taw);
assert(all(Taw>Te));
assert(all(abs(recoverH-edgeH-factor*U.^2/2)<1e-5));
assert(isnan(airRecoveryTemperature(300,5000,factor)));
[badCp,badH]=airThermalProperties([60 2600 NaN]);
assert(all(isnan(badCp)) && all(isnan(badH)));

%% 4. 等效比热及导数驱动壁温反馈，不改变给定的冷壁基准 h
cp0=airEnthalpySecantCp(1200,300);
[cpWarm,derivative]=airEnthalpySecantCp(1200,500);
assert(cpWarm>cp0 && derivative>0);
delta=0.01;
above=airEnthalpySecantCp(1200,500+delta);
below=airEnthalpySecantCp(1200,500-delta);
assert(abs((above-below)/(2*delta)-derivative)<1e-5);
law=struct('edgeTemperature',1200,'edgePressure',1, ...
    'exponent',0,'referenceViscosity0',1, ...
    'gasConstant',287.05,'thermalPropertyMode','thermallyPerfect', ...
    'referenceCp0',cp0);
[hWarm,dh]=evaluateWallConvection(100,1200,500,law);
assert(abs(hWarm/100-cpWarm/cp0)<1e-10);
assert(abs(dh/100-derivative/cp0)<1e-10);

%% 5. 非干空气不能误用本模型；旧版定比热可明确选择
fs=struct('R',287.05,'gamma',1.4);
[mode,scope]=selectAirThermalModel(fs,struct());
assert(strcmp(mode,'thermallyPerfect') && contains(scope,'轴向推进'));
cfg.heating.airThermalModel='caloricallyPerfect';
assert(strcmp(selectAirThermalModel(fs,cfg),'caloricallyPerfect'));
fs.R=320;
failed=false;
try, selectAirThermalModel(fs,struct()); catch, failed=true; end
assert(failed);
fprintf('Air thermal properties passed: cp(300)=%.1f, cp(1500)=%.1f J/(kg K).\n', ...
    cp(1),cp(3));
end
