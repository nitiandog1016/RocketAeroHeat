function result = testIdealAirThermodynamicState(folder)
%TESTIDEALAIRTHERMODYNAMICSTATE 新物性入口的原生对照、反解和兼容性定向检查。
% 只核验物性模块，输出小型CSV/JSON，流场与热流资格另行验证。
if nargin < 1, folder = fullfile('outputs','v153_validation','matlab01'); end
assert(~isfolder(folder),'已有记录，禁止覆盖测试。');
mkdir(folder);
fixture = readtable(fullfile('outputs','v153_validation','native03','core.csv'));
started = tic;

%% 1. 读取旧入口基线，核对新入口全部字段及四种密度和温度边界
legacyTemperature = [80 150 200 300 1000 2000 2500 NaN];
[oldCp,oldH,oldS] = airThermalProperties(legacyTemperature);
value = idealAirThermodynamicState('rhoT',fixture.rho,fixture.T);
fields = fixture.Properties.VariableNames(3:end);
export = table;
nativeError = 0;
for index = 1:numel(fields)
    name = fields{index};
    nativeError = max(nativeError,relativeError(value.(name),fixture.(name)));
    export.(name) = value.(name);
end
assert(nativeError<=2e-11,'原生物性对照失败。');
assert(strcmp(value.model,'Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'));
assert(~value.productionQualified && isequal(value.temperatureRangeK,[60 2000]));
writetable(export,fullfile(folder,'states.csv'));

%% 2. 从每个状态通过六种公开状态对恢复原状态
inputs = {{'rhoT',value.rho,value.T},{'PT',value.p,value.T}, ...
    {'Prho',value.p,value.rho},{'rhoe',value.rho,value.e}, ...
    {'hs',value.h,value.s},{'Ps',value.p,value.s}};
inverseError = 0;
for index = 1:numel(inputs)
    recovered = idealAirThermodynamicState(inputs{index}{:});
    for fieldIndex = 1:numel(fields)
        name = fields{fieldIndex};
        inverseError = max(inverseError,relativeError(recovered.(name),value.(name)));
    end
end
assert(inverseError<=2e-10,'状态恢复失败。');
matrixState = idealAirThermodynamicState('rhoT',.02,[80 300;1000 2000]);
assert(isequal(size(matrixState.cp),[2 2]));
columnState = idealAirThermodynamicState('PT',[100;1000],300);
assert(isequal(size(columnState.cp),[2 1]));

%% 3. 含给定湍黏度的输运导数、能量导数和SST能量剥离代数
temperature = [80 300 1500];
rho = [.001 .02 .2];
eddy = [0 1e-5 1e-4];
step = temperature*1e-5;
base = idealAirThermodynamicState('rhoT',rho,temperature,eddy);
minus = idealAirThermodynamicState('rhoT',rho,temperature-step,eddy);
plus = idealAirThermodynamicState('rhoT',rho,temperature+step,eddy);
derivativeError = max([relativeError((plus.e-minus.e)./(2*step),base.cv), ...
    relativeError((plus.h-minus.h)./(2*step),base.cp), ...
    relativeError((plus.cp-minus.cp)./(2*step),base.dcp_dt), ...
    relativeError((plus.lambda-minus.lambda)./(2*step),base.dlambda_dt)]);
assert(derivativeError<=3e-6,'导数差分失败。');
velocitySquared = [1e4 1e6 4e6];
turbulentK = [0 .1 200];
rhoEnergy = rho.*(base.e+.5*velocitySquared+turbulentK);
recovered = idealAirThermodynamicState('rhoe',rho,rhoEnergy./rho-.5*velocitySquared-turbulentK);
assert(relativeError(recovered.T,temperature)<=2e-10);

%% 4. 非法状态整次拒绝、维度拒绝及旧工程调用结果核对
invalid = {{'rhoT',0,300},{'rhoT',1,59.999},{'rhoT',1,2000.001}, ...
    {'rhoT',1,[300 NaN]},{'PT',-1,300},{'Prho',1,0}, ...
    {'rhoe',1,min(value.e)-1},{'rhoe',1,max(value.e)+1}, ...
    {'hs',min(value.h)-1,1},{'Ps',1,1e10},{'wrong',1,1}, ...
    {'rhoT',[1 1],[80;300]},{'rhoT',1,300,-1},{'rhoT',1,300,Inf}, ...
    {'rhoT',1,300,1e308},{'rhoT',1,300,[0 0]},{'rhoT',Inf,300}, ...
    {'rhoT',1e308,300}};
for index = 1:numel(invalid)
    caught = false;
    sentinel = 37;
    try
        sentinel = idealAirThermodynamicState(invalid{index}{:});
    catch
        caught = true;
    end
    assert(caught && isequal(sentinel,37),'非法状态未完整拒绝。');
end
[newCp,newH,newS] = airThermalProperties(legacyTemperature);
assert(isequaln(oldCp,newCp) && isequaln(oldH,newH) && isequaln(oldS,newS));
assert(all(isfinite(oldCp(1:7))) && isnan(oldCp(8)));

%% 5. 保存可独立复算的结果和实际耗时
result = struct('states',height(fixture),'statePairs',6,'nativeMaximumError',nativeError, ...
    'inverseMaximumError',inverseError,'derivativeMaximumError',derivativeError, ...
    'invalidStatesRejected',numel(invalid),'legacyCompatible',true, ...
    'elapsedSeconds',toc(started),'numericalPassed',true,'newQualifiedFlowArea',0);
fid = fopen(fullfile(folder,'summary.json'),'w','n','UTF-8');
assert(fid>=0); closeFile = onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
disp(result);
end

function value = relativeError(actual,expected)
% 所有字段使用事前登记的分母下限1。
error = abs(actual-expected)./max(abs(expected),1);
value = max(error(:));
end
