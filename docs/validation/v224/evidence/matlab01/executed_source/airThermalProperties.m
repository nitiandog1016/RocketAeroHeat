function [cp,enthalpy,standardEntropy] = airThermalProperties(T)
%AIRTHERMALPROPERTIES 干空气温度相关的定压比热、质量比焓和标准熵。
% 模块作用：用 N2/O2 的 NASA 七系数多项式及 Ar 单原子气体近似，
% 按固定质量分数混合。焓/熵的零点任意，只用差值求换热和激波。
% 适用 80~2500 K；150~200 K 为 NASA7 低温多项式的少量外推。
% 80~150 K 采用在 150 K 连续衔接的振动冻结定比热近似：
% h=h150+cp150*(T-150)，s=s150+cp150*ln(T/150)。在此区间
% N2/O2 转动自由度仍充分激发，振动自由度基本冻结。此近似仅适用
% 稀薄且保持气相的固定组分空气；低于 80 K 或高于 2500 K 返回 NaN。
% 系数取 NASA/Glenn 常用的 200~1000 与 1000~6000 K 区间形式。

%% 1. 校验温度并预留输出；已有无解面元的 NaN 原样传递
validateattributes(T,{'numeric'},{'real'});
T=double(T);
cp=nan(size(T)); enthalpy=cp; standardEntropy=cp;
valid=isfinite(T) & T>=150 & T<=2500;
belowAnchor=isfinite(T) & T>=80 & T<150;
if ~any(valid(:)) && ~any(belowAnchor(:)), return; end
if any(valid(:))
tOriginal=T(valid);
t=tOriginal(:);

%% 2. 分别求 N2、O2、Ar，再按干空气质量分数混合
R=8.31446261815324;
mass=[0.755 0.232 0.013];
mass=mass/sum(mass);
molecularMass=[0.0280134 0.0319988 0.039948];
low=[3.53100528 -1.23660987e-4 -5.02999433e-7 2.43530612e-9 -1.40881235e-12 -1046.97628 2.96747468; ...
     3.78245636 -2.99673416e-3  9.84730201e-6 -9.68129509e-9  3.24372837e-12 -1063.94356 3.65767573; ...
     2.5         0               0              0              0              0           0];
high=[2.95257626 1.39690040e-3 -4.92631603e-7 7.86010367e-11 -4.60755321e-15 -923.948688 5.87188762; ...
      3.28253784 1.48308754e-3 -7.57966669e-7 2.09470555e-10 -2.16717794e-14 -1088.45772 5.45323129; ...
      2.5        0              0              0              0               0           0];
sumCp=zeros(size(t)); sumH=sumCp; sumS=sumCp;
for gas=1:3
    a=repmat(low(gas,:),numel(t),1);
    a(t>=1000,:)=repmat(high(gas,:),nnz(t>=1000),1);
    cpOverR=a(:,1)+a(:,2).*t+a(:,3).*t.^2+ ...
        a(:,4).*t.^3+a(:,5).*t.^4;
    hOverRT=a(:,1)+a(:,2).*t/2+a(:,3).*t.^2/3+ ...
        a(:,4).*t.^3/4+a(:,5).*t.^4/5+a(:,6)./t;
    sOverR=a(:,1).*log(t)+a(:,2).*t+a(:,3).*t.^2/2+ ...
        a(:,4).*t.^3/3+a(:,5).*t.^4/4+a(:,7);
    factor=mass(gas)*R/molecularMass(gas);
    sumCp=sumCp+factor*cpOverR;
    sumH=sumH+factor*t.*hOverRT;
    sumS=sumS+factor*sOverR;
end
cp(valid)=reshape(sumCp,size(tOriginal));
enthalpy(valid)=reshape(sumH,size(tOriginal));
standardEntropy(valid)=reshape(sumS,size(tOriginal));
end

%% 3. 低温振动冻结近似；从 150 K 锚点积分焓和熵，保证连续
if any(belowAnchor(:))
    [cp150,h150,s150]=airThermalProperties(150);
    cp(belowAnchor)=cp150;
    enthalpy(belowAnchor)=h150+cp150*(T(belowAnchor)-150);
    standardEntropy(belowAnchor)=s150+cp150*log(T(belowAnchor)/150);
end
end
