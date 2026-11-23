function t=thermistorResistanceModel(c,internalK)
% 模块一：根据 B 珠表面积和附录的两个共焦椭球坐标，计算壳层的径向导热热阻。
% Ti 表示电阻标定的等效内部平均温度，Tw 表示模型推算的珠表面温度。
assert(c.speedMps==3 && c.airK==293.15 && internalK==305.15 && c.sourceAreaM2==3.094e-6);
xiW=1.80;xiI=1.25;L=.3430e-3;ks=8.2*4184/3600;
shape=xiW^2-1+xiW^2*sqrt(xiW^2-1)*asin(1/xiW);
focus=sqrt(c.sourceAreaM2/(2*pi*shape));
integral=log((xiI+1)*(xiW-1)/((xiI-1)*(xiW+1)));
shellResistance=integral/(8*pi*focus*ks);
% 模块二：只用附录式（1′）确定测温修正的换热系数，三维珠部换热系数另由面元积分。
re=c.speedMps*L*c.rho/c.mu;
h=c.k/L*(2/3+.3*sqrt(re));
ratio=shellResistance*h*c.sourceAreaM2;factor=1/(1+ratio);
surfaceK=c.airK+factor*(internalK-c.airK);
% 模块三：用长度确定焦半距和使用膜温物性分别估计输入选择的影响。
focusFromLength=3*L*shape/(4*xiW*(xiW^2-1));
filmRatio=ratio;
for j=1:4
    filmK=c.airK+(internalK-c.airK)/(2*(1+filmRatio));
    air=phaseAwareAirState('PT',c.airPa,filmK,'lemmonAnalytic');
    filmRe=c.speedMps*L*air.densityKgm3/air.dynamicViscosityPas;
    filmRatio=shellResistance*c.sourceAreaM2*air.conductivityWmK/L*(2/3+.3*sqrt(filmRe));
end
t=struct('method','坪内、佐藤附录式（1′）、（2′），共焦等温壳导热与电阻等效测温', ...
    'doi','10.1299/kikai1938.25.810','xiWall',xiW,'xiInternal',xiI,'assumedMinorAxisRatio',2, ...
    'semifocalM',focus,'lengthM',L,'areaM2',c.sourceAreaM2,'solidConductivityWmK',ks,'shellResistanceKW',shellResistance, ...
    'reynoldsL',re,'appendixHWm2K',h,'ratio',ratio,'scale',factor,'internalK',internalK,'ambientK',c.airK,'surfaceK',surfaceK, ...
    'focusFromLengthM',focusFromLength,'ratioFromLength',ratio*focusFromLength/focus, ...
    'filmK',filmK,'filmReynoldsL',filmRe,'filmRatio',filmRatio,'reportedInternalUncertaintyK',NaN, ...
    'scope','B 珠名义内部尺寸，圆盘材料导热率。图 12 为模型曲线，实际内部几何和单点测温误差未知。');
end
