function [cpMean,derivativeWall] = airEnthalpySecantCp(Taw,Twall)
%AIRENTHALPYSECANTCP 把空气焓差换算成当前温区的等效比热。
% 模块作用：cpMean=[h(Taw)-h(Twall)]/(Taw-Twall)，供 q=h_dT*(Taw-Twall)
% 形式的现有壁温求解器使用；第二输出是固定 Taw 时对壁温的导数。
% 两温度接近时使用端点比热平均，避免 0/0 和数值消减。

%% 1. 广播标量壁温，并计算两端焓和比热
delta=double(Taw)-double(Twall);
Taw=double(Taw)+zeros(size(delta));
Twall=double(Twall)+zeros(size(delta));
[cpAw,hAw]=airThermalProperties(Taw);
[cpWall,hWall]=airThermalProperties(Twall);
cpMean=(hAw-hWall)./delta;
derivativeWall=(cpMean-cpWall)./delta;

%% 2. 近等温极限；无解值仍保留 NaN
near=isfinite(delta) & abs(delta)<1e-4.*max(1,abs(Taw));
cpMean(near)=(cpAw(near)+cpWall(near))/2;
if nargout>1
    step=0.1;
    [cpPlus,~]=airThermalProperties(Taw(near)+step);
    [cpMinus,~]=airThermalProperties(Taw(near)-step);
    derivativeWall(near)=0.5*(cpPlus-cpMinus)/(2*step);
end
end
