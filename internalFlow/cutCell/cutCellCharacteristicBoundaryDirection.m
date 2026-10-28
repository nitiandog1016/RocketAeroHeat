function [direction,report]=cutCellCharacteristicBoundaryDirection(P,dP,b,normal,gas)
%CUTCELLCHARACTERISTICBOUNDARYDIRECTION 对原局部特征边界的当前分支求方向导数。
% 模块一：原解析器先校验基态/储库/声速转换，给定储库及压力保持不变。
% 模块二：传播J+、声速、熵和切向速度，包含压力出口及回流重分类。
% 模块三：零面流速沿用内部熵/切向分支，不称跨分支全局光滑导数。
validateattributes(dP,{'numeric'},{'real','finite','numel',5});dP=reshape(dP,1,5);P=reshape(P,1,5);
[face,tag]=resolveNSCharacteristicBoundary3D(P,b,normal,gas);normal=reshape(normal,1,3);
g=gas.gamma;gm=g-1;a=sqrt(g*gas.R*P(5));un=dot(P(2:4),normal);
dun=dot(dP(2:4),normal);da=.5*a*dP(5)/P(5);dJ=dun+2*da/gm;
dlogK=dP(5)/P(5)-gm*dP(1)/P(1);dTangent=dP(2:4)-dun*normal;
switch tag.branch
    case 'supersonic-outflow'
        direction=dP;
    case 'supersonic-inflow'
        direction=zeros(1,5);
    otherwise
        ab=sqrt(g*gas.R*face(5));
        if tag.pressureInputApplied
            dlogRho=-dlogK/g;dT=-face(5)*dlogRho;dab=.5*ab*dT/face(5);
            dub=dJ-2*dab/gm;drho=face(1)*dlogRho;dv=dTangent+dub*normal;
        else
            dub=.5*dJ;dab=.25*gm*dJ;dT=2*face(5)*dab/ab;
            if strcmp(tag.branch,'reservoir-backflow'),selectedEntropy=0;selectedTangent=zeros(1,3);
            else,selectedEntropy=dlogK;selectedTangent=dTangent;end
            drho=face(1)/gm*(2*dab/ab-selectedEntropy);dv=selectedTangent+dub*normal;
        end
        direction=[drho dv dT];
end
assert(all(isfinite(direction)),'CutDirection:Boundary','特征边界方向导数非有限。');
report=struct('basePrimitive',face,'branch',tag.branch,'isInflow',tag.isInflow, ...
    'pressureInputApplied',tag.pressureInputApplied,'reclassifiedByFaceVelocity',tag.reclassifiedByFaceVelocity, ...
    'zeroTangentialPolicy',tag.zeroTangentialPolicy,'selectedBranchOnly',true, ...
    'prescribedReservoirOrPressureDifferentiated',false,'fullNSJacobianClaimed',false);
end
