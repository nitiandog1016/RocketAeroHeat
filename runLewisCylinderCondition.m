function [mesh,stations,report]=runLewisCylinderCondition(staticPressurePa,airTemperatureK,wallTemperatureK)
%RUNLEWISCYLINDERCONDITION 计算新给定热边界下的圆柱工程场，沿用原面及六测带。
c=lewisCylinderCondition(staticPressurePa,airTemperatureK,wallTemperatureK);
[mesh,stations,report]=runLewisCylinderBenchmark();
% 从真实三角形的边界角度确定每个面带，不假定面列表的排列顺序。
f=mesh.F; v=mesh.V;
angle=rad2deg(unwrap(atan2(reshape(v(f,3),size(f)),reshape(v(f,2),size(f))),[],2));
center=mod((max(angle,[],2)+min(angle,[],2))/2+180,360)-180;
valid=isfinite(mesh.lewisCylinder.nusselt);
gauge=nan(size(valid));gauge(valid)=c.dynamicPressurePa*meanPressureCoefficient(center(valid),1);
h=mesh.lewisCylinder.nusselt*c.conductivityWmK/c.diameterM;
q=h*(c.adiabaticWallK-c.wallTemperatureK);
mesh.lewisCondition=c;
mesh.lewisWall=struct('pressurePa',c.referencePressurePa+gauge,'gaugePa',gauge,'heatFluxWm2',q);
mesh.evidenceIdentity.scope='给定边界的低速圆柱；中央附着区工程估计，压力外缘和有量纲热流尚无同工况独立验证';
% 试验表中的 Nu 留作独立参照。以下有量纲量只属于当前给定的研究工况。
stations.pressurePa=c.referencePressurePa+c.dynamicPressurePa*meanPressureCoefficient(stations.angleDeg,report.bandWidthDeg);
stations.heatFluxWm2=stations.predictedNu*c.conductivityWmK/c.diameterM*(c.adiabaticWallK-c.wallTemperatureK);
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
report.wallCoverage=summarizeComputedWall(area,isfinite(gauge),isfinite(q),q);
report.absolutePressureHeatAreaM2=report.wallCoverage.jointComputedAreaM2;
end

function cp=meanPressureCoefficient(centerDeg,widthDeg)
% Hiemenz 外缘速度引自 NACA TM 1432 印刷页20；用伯努利壁压近似并对测带积分。
% x=s/D；对速度平方积分，保留迎流驻点两侧的对称性。
x=deg2rad(centerDeg)/2;d=deg2rad(widthDeg)/4;
velocity=[-1.5144 0 -2.1709 0 3.6314 0];
primitive=polyint(conv(velocity,velocity));
cp=1-(polyval(primitive,x+d)-polyval(primitive,x-d))/(2*d);
end
