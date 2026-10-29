function [cp,h,s,matching]=surfaceFlowAirThermalProperties(T)
%SURFACEFLOWAIRTHERMALPROPERTIES 原项目cp(T)的分段连续焓/标准熵积分约定。
% 模块一：原airThermalProperties保持不变；只在本流线模型显式选择此约定。
% 模块二：保持80–1000K原焓/熵零点；高温段只加常数，使h和s在1000K连续。
% 模块三：cp完全保留原值，不裁温、不用拟合误差放宽沿程守恒门槛。
% 原有限位数多项式的h/s积分常数并非数值精确衔接；跨段时必须先统一
% 积分约定。旧外流总状态不能直接混用，须从物理p/T/rho/U显式转换。
persistent offsets
if isempty(offsets)
    knot=1000;delta=eps(knot);
    [leftCp,leftH,leftS]=airThermalProperties(knot-delta);
    [~,rightH,rightS]=airThermalProperties(knot);
    offsets=struct('knotK',knot,'highEnthalpyOffsetJPerKg',leftH+leftCp*delta-rightH, ...
        'highEntropyOffsetJPerKgK',leftS+leftCp/knot*delta-rightS, ...
        'convention','continuous-project-cp-v1','originalProductionLibraryModified',false);
end
[cp,h,s]=airThermalProperties(T);
high=isfinite(T)&T>=offsets.knotK;
h(high)=h(high)+offsets.highEnthalpyOffsetJPerKg;
s(high)=s(high)+offsets.highEntropyOffsetJPerKgK;
matching=offsets;
end
