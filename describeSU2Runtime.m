function info=describeSU2Runtime(receipt)
%DESCRIBESU2RUNTIME 从实际进程回执描述求解器版本和受控改动，不由当前默认配置猜测历史。
% 模块一：旧回执缺失身份时明确未知；新回执展示准确包类型。
% 模块二：实验恢复与只读诊断分开，任何运行包均不自动赋予物理资格。
info=struct('version','未记录','flavor','未记录','mode','未记录','physicsQualified',false);
info.physicsCaption='laminar NS / ideal gas / Sutherland / Pr0.72';
if isfield(receipt,'executableVersion'),info.version=receipt.executableVersion;end
if isfield(receipt,'runtimeFlavor'),info.flavor=receipt.runtimeFlavor;end
switch info.flavor
    case 'official'
        info.mode='native';info.method='官方原生更新';
    case 'rah-axial-wall-v1'
        info.mode='axial-wall-temperature';info.method='受控逐节点空间壁温；原生流场更新';
    case 'rah-coolprop-v1'
        info.mode='compressible-heos-air';info.method='受控CoolProp 7.2.0 HEOS Air；可压缩变比热；Sutherland黏度与k=μcp/Pr，原生流场更新';
        info.physicsCaption='层流 NS / HEOS Air / Sutherland / k=μcp(T,ρ)/Pr0.72';
    case 'rah-frozen-energy-v1'
        assert(isfield(receipt,'frozenEnergyMode'),'冻结能量运行时缺少真实模式回执。');
        info.mode=receipt.frozenEnergyMode;
        switch info.mode
            case 'diagnostic',info.method='冻结能量诊断：更新保持原生';
            case 'apply',info.method='实验正内能恢复：仅原生零松弛且端点内能增加时沿原增量更新';
            otherwise,error('未知冻结能量模式，不能显示为官方原生解。');
        end
    otherwise
        info.method=['求解器更新身份未核实：' info.flavor];
        info.physicsCaption='求解物理配置未核实';
end
info.caption=['SU2 ' info.version];
end
