function state=recoverThermallyPerfectEdgeState(pressure,totalPressure,totalTemperature,R)
%RECOVERTHERMALLYPERFECTEDGESTATE 由静压/总压恢复干空气等熵外缘状态。
% 模块作用：给定来自同一温度相关激波模型的 p、p0、T0，以 NASA7
% 标准熵差反解静温，再以总焓差恢复速度、局部音速与 Mach。
% 只在固定组分干空气、80~2500 K 范围使用；越界或 p>p0 留作 NaN。

%% 1. 校验形状与热物性适用范围
validateattributes(pressure,{'numeric'},{'real'});
validateattributes(totalPressure,{'numeric'},{'real'});
validateattributes(totalTemperature,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(R,{'numeric'},{'real','finite','scalar','positive'});
if ~isequal(size(pressure),size(totalPressure))
    error('RocketAeroHeat:ThermalEdge:Size','静压与总压的数组大小必须相同。');
end
if abs(R-287.05)>5
    error('RocketAeroHeat:ThermalEdge:Gas','温度相关恢复只适用于干空气。');
end
p=double(pressure); p0=double(totalPressure);
blank=nan(size(p));
state=struct('valid',false(size(p)), ...
    'method','固定组分干空气 NASA7/低温振动冻结标准熵与总焓反算', ...
    'scope','80~2500 K；80~150 K 为振动冻结近似；仅由同模型总状态恢复', ...
    'temperature',blank,'density',blank,'speed',blank, ...
    'Mach',blank,'speedOfSound',blank,'specificHeat',blank, ...
    'entropyResidual',blank);
[~,h0,s0]=airThermalProperties(totalTemperature);
[~,~,sMin]=airThermalProperties(80);
if ~isfinite(h0) || ~isfinite(s0), return; end
target=s0-R*log(p0./p);
valid=isfinite(p) & isfinite(p0) & p>0 & p0>=p & ...
    target>=sMin-1e-9 & target<=s0+1e-9;
if ~any(valid(:)), return; end

%% 2. 熵随温度单调增加，对所有有效面元一起二分反解静温
lo=80*ones(nnz(valid),1);
hi=totalTemperature*ones(nnz(valid),1);
wanted=target(valid); wanted=wanted(:);
for iteration=1:47
    mid=0.5*(lo+hi);
    [~,~,sMid]=airThermalProperties(mid);
    higher=sMid>wanted;
    hi(higher)=mid(higher);
    lo(~higher)=mid(~higher);
end
T=0.5*(lo+hi);
[cp,h,s]=airThermalProperties(T);
gamma=cp./(cp-R);
U=sqrt(max(0,2*(h0-h)));
a=sqrt(gamma.*R.*T);
good=isfinite(T) & isfinite(U) & isfinite(a) & cp>R & a>0;
ids=find(valid);
ids=ids(good);
state.valid(ids)=true;
state.temperature(ids)=T(good);
state.density(ids)=p(ids)./(R*T(good));
state.speed(ids)=U(good);
state.Mach(ids)=U(good)./a(good);
state.speedOfSound(ids)=a(good);
state.specificHeat(ids)=cp(good);
state.entropyResidual(ids)=s0-s(good)-R*log(p0(ids)./p(ids));
end
