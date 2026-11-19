function report = summarizeComputedWall(area, pressureComputed, heatComputed, heatFluxWm2)
%SUMMARIZECOMPUTEDWALL 汇总同一工况、同序原三角面上的压力与热流计算结果。
% 模块一：调用者按各量的来源提供计算标志。给定值与超域值不计入计算面积。
validateattributes(area,{'numeric'},{'real','vector','nonempty','finite','positive'});
area=double(area(:));
n=numel(area);
validateattributes(pressureComputed,{'logical'},{'vector','numel',n});
validateattributes(heatComputed,{'logical'},{'vector','numel',n});
validateattributes(heatFluxWm2,{'numeric'},{'real','vector','numel',n});
p=pressureComputed(:); q=heatComputed(:); flux=double(heatFluxWm2(:));
assert(~any(isinf(flux)) && all(isfinite(flux(q))), ...
    'RocketAeroHeat:Coverage:HeatStatus','标为已计算的热流须有限；缺失值用 NaN 表示。');
total=sum(area);
assert(isfinite(total),'RocketAeroHeat:Coverage:Area','原面总面积超出数值范围。');

% 模块二：逐原面求交。即使某个缺项面很小，也按面标记整表面尚未算全。
joint=p & q;
report.faceCount=n;
report.totalArea=total;
report.jointComputedFaceCount=nnz(joint);
report.jointComputedAreaM2=sum(area(joint));
report.jointMissingAreaM2=sum(area(~joint));
report.fullSurfaceComputed=all(joint);
report.pressureComputedAreaFraction=sum(area(p))/total;
report.heatComputedAreaM2=sum(area(q));
report.heatComputedAreaFraction=report.heatComputedAreaM2/total;

% 模块三：功率只依赖热流。部分积分可为负值或零；没有可积分面时保持未知。
report.calculatedHeatPowerW=NaN;
report.wholeSurfaceHeatPowerW=NaN;
if any(q)
    power=sum(area(q).*flux(q));
    if isfinite(power)
        report.calculatedHeatPowerW=power;
        if all(q),report.wholeSurfaceHeatPowerW=power;end
    end
end
if isfinite(report.wholeSurfaceHeatPowerW)
    powerNote=sprintf('整表面热功率 %.6g W',report.wholeSurfaceHeatPowerW);
elseif isfinite(report.calculatedHeatPowerW)
    powerNote=sprintf('已算面热功率 %.6g W，整表面热功率未知',report.calculatedHeatPowerW);
else
    powerNote='整表面热功率未知';
end
report.coverageMessage=sprintf(['壁压与热流同面已算 %d/%d 面，%.6g / %.6g m²（%.2f%%），' ...
    '缺项面积 %.6g m²。%s。'],report.jointComputedFaceCount,n,report.jointComputedAreaM2,total, ...
    100*(report.jointComputedAreaM2/total),report.jointMissingAreaM2,powerNote);
end
