function result=estimateSonicBoattailPressure(cylinderPressurePa,boattailAngleDeg)
%ESTIMATESONICBOATTAILPRESSURE 复算NASA TN D-2761的声速流线膨胀压力候选。
% 模块用途：为短尾锥厚湍流边界层提供独立工程压力对照，不接入默认生产链。
% 输入：同工况、收缩前柱面的静压[Pa]；直尾锥收缩半角[deg]，范围0到21度。
% 输出：由M=1流线等熵膨胀得到的候选压力比和候选壁压。
% 这里的M=1是边界层内的假想声速流线，不是外缘、远场或整片流场的Mach。
% 本模块不求热流、外缘/恢复温度、壁温、横流、分离泡或底面流动。

%% 模块一：固定原报告定比热空气曲线及来源，限制为条件研究候选
validateattributes(cylinderPressurePa,{'numeric'},{'real','finite','scalar','positive'}, ...
    mfilename,'cylinderPressurePa');
validateattributes(boattailAngleDeg,{'numeric'},{'real','finite','nonempty','vector','>=',0,'<=',21}, ...
    mfilename,'boattailAngleDeg');
gamma=1.4;angles=double(boattailAngleDeg);mach=ones(size(angles));
result=struct('method','NASA TN D-2761：声速流线二维Prandtl-Meyer膨胀压力近似', ...
    'sourceURL','https://ntrs.nasa.gov/api/citations/19650013004/downloads/19650013004.pdf', ...
    'sourcePrintedPages',[7 8 9 33 34], ...
    'status','conditional-engineering-correlation-not-production-qualified', ...
    'scope','原报告短轴对称直尾锥、厚湍流边界层的压力对照；未认证任意箭体、攻角或分离区', ...
    'gamma',gamma,'cylinderPressurePa',double(cylinderPressurePa), ...
    'boattailAngleDeg',angles,'sourceNominalTailLengthToDiameter',.35, ...
    'upstreamPressureIsGivenToThisModule',true,'upstreamSonicMachIsApproximation',true, ...
    'externalEdgeStateSolved',false,'heatFluxSolved',false,'wallTemperatureSolved',false, ...
    'productionQualified',false,'physicalValidationComplete',false,'newProductionSolvedAreaM2',0);

%% 模块二：从声速端点直接求根，不把M=1篡改为1+eps来绕过旧超声速接口
for k=1:numel(angles)
    if angles(k)==0,continue;end
    target=deg2rad(angles(k));upper=2;
    while pm(upper)<target,upper=2*upper;end
    mach(k)=fzero(@(m)pm(m)-target,[1 upper],optimset('TolX',1e-14,'Display','off'));
end
factor=(1+.5*(gamma-1))./(1+.5*(gamma-1)*mach.^2);
ratio=factor.^(gamma/(gamma-1));

%% 模块三：只发布候选压力；声速流线代数量不冒充实际外缘热力状态
result.sonicStreamlineMachProxy=mach;
result.candidatePressureRatio=ratio;
result.candidateWallPressurePa=double(cylinderPressurePa)*ratio;
result.edgeTemperatureK=nan(size(angles));
result.recoveryTemperatureK=nan(size(angles));
result.wallTemperatureSolvedK=nan(size(angles));
result.heatFluxWm2=nan(size(angles));
result.limitations={'柱面静压须由同工况测量或本程序上游解提供，不能换成任意远场压力', ...
    '原报告认为18/21度尾锥有分离迹象，本候选不解分离或底区', ...
    '没有推导通用厚度阈值，也没有认证真实20度攻角箭体', ...
    '关联式复算与实验比较、完整NS/RANS验证和生产面积资格分别验收'};
    function value=pm(m)
        root=sqrt((m-1)*(m+1));
        value=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1))*root)-atan(root);
    end
end
