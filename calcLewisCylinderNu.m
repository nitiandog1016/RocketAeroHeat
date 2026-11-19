function nu = calcLewisCylinderNu(centerDeg, widthDeg)
%CALCLEWISCYLINDERNU 计算 NASA Lewis 光滑圆柱低湍流工况的带平均 Nu_D。
% 模块一：迎流驻点为零度，角度取有符号的圆周角；宽度为整个测带的弧角。
% 输入宽度为零时求中心点值。整带超出 Frössling 文献的 ±55° 范围即返回 NaN。
validateattributes(centerDeg,{'numeric'},{'real','finite','nonempty'});
validateattributes(widthDeg,{'numeric'},{'real','scalar','finite','nonnegative'});
eligible = abs(double(centerDeg)) + double(widthDeg)/2 <= 55;
nu = nan(size(centerDeg));
% 模块二：NACA TM 1432 印刷页 20，空气 Pr=0.7，x=s/D=θ/2。
% 长方形铜带在周向等弧长加权；用二次和四次项的解析均值避免小带宽相减。
x = deg2rad(double(centerDeg(eligible)))/2;
d = deg2rad(double(widthDeg))/4;
meanSquare = x.^2 + d^2/3;
meanFourth = x.^4 + 2*x.^2*d^2 + d^4/5;
nu(eligible) = sqrt(177000)*(0.9449 - 0.5100*meanSquare - 0.5956*meanFourth);
end
