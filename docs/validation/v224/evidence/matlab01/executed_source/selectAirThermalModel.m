function [mode,scope] = selectAirThermalModel(freestream,cfg)
%SELECTAIRTHERMALMODEL 明确热流使用的空气热物性与适用条件。
% 模块作用：默认温度相关干空气焓；可设为 caloricallyPerfect
% 重现旧版定比热结果。本选择作用于鼻端法向激波、非零攻角局部
% 二维切楔激波、非零攻角背风等熵膨胀及受限前体闭合、驻点焓差
% 和壁面换热。零攻角切锥与轴向波系优先使用同一 NASA7 模型；
% 80~150 K 使用振动冻结定比热低温近似；若轴向状态低于 80 K，
% 则整条侧壁外流回退定比热并显式标记。

%% 1. 读取配置并检查气体常数
mode='thermallyPerfect';
if nargin>=2 && isfield(cfg,'heating') && ...
        isfield(cfg.heating,'airThermalModel')
    mode=char(cfg.heating.airThermalModel);
end
if ~ismember(mode,{'thermallyPerfect','caloricallyPerfect'})
    error('RocketAeroHeat:AirThermal:Mode', ...
        'heating.airThermalModel 必须是 thermallyPerfect 或 caloricallyPerfect。');
end
if strcmp(mode,'thermallyPerfect')
    if ~isfield(freestream,'R') || ~isfield(freestream,'gamma') || ...
            abs(freestream.R-287.05)>5 || abs(freestream.gamma-1.4)>0.05
        error('RocketAeroHeat:AirThermal:Gas', ...
            '温度相关热物性只适用于近似干空气的 R、gamma；其他气体请选择 caloricallyPerfect。');
    end
    scope=['干空气 N2/O2/Ar 固定组分 NASA7 比热/焓，150~2500 K；' ...
        '80~150 K 采用与 150 K 焓/熵连续的振动冻结定比热近似；' ...
        '150~200 K 少量外推；80~150 K 的 Sutherland 黏度也是工程外推；' ...
        '不含解离、电离；鼻端单点法向激波用变比热，' ...
        '非零攻角局部二维切楔激波/背风膨胀/前体闭合也用变比热；' ...
        '低于 80 K 的非零攻角局部膨胀仅以定比热外缘占位维持流线、热流留空；' ...
        '零攻角切锥/轴向推进优先使用同一干空气模型；轴向域外则侧壁整条退回定比热'];
else
    scope='定比热 cp=gamma*R/(gamma-1)，与旧版热流计算一致；高温误差未修正';
end
end
