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
        if isfield(receipt,'physicsModel')
            p=receipt.physicsModel;assert(strcmp(p.model,'SST-2003m'),'未知的官方湍流合同。');
            info.mode='sst-2003m';info.physicsModel=p;
            info.method=sprintf('官方SST-2003m；全湍流研究假设，I=%.4g，μt/μ=%.4g，Prt=%.4g；%s', ...
                p.intensity,p.viscosityRatio,p.Prt,p.provenance);
            info.physicsCaption=sprintf('RANS / SST-2003m / 定比热 / Sutherland / Prt%.4g；转捩未知',p.Prt);
            if isfield(p,'outputStateContract')
                assert(strcmp(p.outputStateContract,'native-noncompact-flow-primitive-before-turbulence-v1'));
                info.method=[info.method '；p/T/热流为原生流动阶段，保存k/ω为随后湍流更新态，分别审核'];
            end
            if isfield(p,'diagnosticOnly') && p.diagnosticOnly
                info.method=[info.method '；仅一次阶段诊断，热状态使用更新流场与输入k；源后验失败，未追溯放行'];
                info.physicsCaption=[info.physicsCaption ' / 单步阶段诊断'];
            end
            if isfield(p,'numericalStudy')
                % 同一物理模型下的推进参数是数值方案，不冒充新增物理方法或已收敛结果。
                study=p.numericalStudy;
                if strcmp(study.studyID,'v146-common-SST149-first-order-startup')
                    % 实际流动阶数进入探针与导出，不能把一阶启动场标成二阶热环境。
                    switch study.variant
                        case 'control',label='二阶流动重构控制';assert(study.flowSpatialOrder==2);
                        case 'first_order',label='一阶流动启动候选；仅作启动研究';assert(study.flowSpatialOrder==1);
                        case 'restore_second_order',label='恢复二阶流动重构；仍待物理验证';assert(study.flowSpatialOrder==2);
                        otherwise,error('未知SST启动对照方案。');
                    end
                elseif strcmp(study.studyID,'v132-tail-rim-viscous-gradient')
                    % 黏性/源项梯度同时影响流动与湍流，不能标成SST单项修复。
                    switch study.variant
                        case 'wls',label='WLS黏性/源项梯度控制；对流重构WLS';
                        case 'gg',label='Green–Gauss黏性/源项梯度候选；流动与SST共同改变，对流重构WLS';
                        otherwise,error('未知梯度对照方案。');
                    end
                else
                    assert(strcmp(study.studyID,'v127-common-SST49-CFL5-vs-0.5'));
                    switch study.variant
                        case 'control',label='CFL上限5控制';
                        case 'cap05',label='CFL上限0.5候选';
                        otherwise,error('未知SST数值对照方案。');
                    end
                end
                info.method=[info.method '；同源数值方案：' label '，伪时间重置'];
                info.physicsCaption=[info.physicsCaption ' / ' label];
            end
            if isfield(p,'initializationStudy')
                % 启发式初猜只改变迭代起点，不称为新的湍流物理模型或合格解。
                study=p.initializationStudy;
                if strcmp(study.studyID,'v128-base-omega-floor-logharmonic-initial-guess')
                    switch study.variant
                        case 'control',label='原始初场控制';
                        case 'repaired',label='尾缘146点ω正值调和初猜';
                        otherwise,error('未知SST初场对照。');
                    end
                else
                    assert(strcmp(study.studyID,'v130-near-wall-diffusion-destruction-initial-guess'));
                    switch study.variant
                        case 'control',label='原始初场控制';
                        case 'column',label='近壁k/ω扩散耗散初猜，rhoE补偿保持热内能';
                        otherwise,error('未知SST近壁初场对照。');
                    end
                end
                info.method=[info.method '；初场对照：' label '，物理/数值参数不改，伪时间重置；不代表物理闭合'];
                info.physicsCaption=[info.physicsCaption ' / ' label];
            end
        end
    case 'rah-bounded-conservative-v1'
        % 原生残差和物理源保持；仅本受控包改变有界散度的守恒Jacobian链式换算。
        p=receipt.physicsModel;assert(strcmp(p.model,'SST-2003m'));
        study=p.numericalStudy;assert(strcmp(study.studyID,'v131-SST-bounded-conservative-chain-rule') && strcmp(study.variant,'bounded'));
        info.mode='sst-bounded-conservative';info.physicsModel=p;
        info.method=sprintf('受控有界标量守恒Jacobian修正；SST-2003m，I=%.4g，μt/μ=%.4g，Prt=%.4g；复用参考初态，伪时间重置；未获物理资格', ...
            p.intensity,p.viscosityRatio,p.Prt);
        info.method=[info.method '；p/T/热流为原生流动阶段，保存k/ω为随后湍流更新态，分别审核'];
        info.physicsCaption=sprintf('RANS / SST-2003m / 定比热 / Sutherland / Prt%.4g；转捩未知 / 受控有界输运',p.Prt);
    case 'rah-axial-wall-v1'
        info.mode='axial-wall-temperature';info.method='受控逐节点空间壁温；原生流场更新';
    case {'rah-coolprop-v1','rah-coolprop-recon-v1','rah-coolprop-recon-v2','rah-coolprop-recon-v3'}
        info.mode='compressible-heos-air';info.method='受控CoolProp 7.2.0 HEOS Air；可压缩变比热；Sutherland黏度与k=μcp/Pr，原生流场更新';
        info.physicsCaption='层流 NS / HEOS Air / Sutherland / k=μcp(T,ρ)/Pr0.72';
        if startsWith(info.flavor,'rah-coolprop-recon-')
            info.method='受控CoolProp 7.2.0 HEOS Air；局部重构EOS异常进入原生回退；Sutherland黏度与k=μcp/Pr';
            info.physicsCaption=[info.physicsCaption ' / 重构异常局部回退'];
        end
        if strcmp(info.flavor,'rah-coolprop-recon-v2'),info.method=[info.method '；上游消息修正未实际链接，等同v1二进制'];end
        if strcmp(info.flavor,'rah-coolprop-recon-v3'),info.method=[info.method '；已核验上游超域消息缺参修正'];end
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
