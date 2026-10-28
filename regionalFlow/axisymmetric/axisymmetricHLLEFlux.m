function [flux,audit]=axisymmetricHLLEFlux(QL,QR,normal,gas,wavePolicy)
%AXISYMMETRICHLLEFLUX 二维面法向Euler四量HLLE通量，逐实际中间态检查正性。
% 模块一：davis用两侧声学界；einfeldt用左/右声学界及Roe平均特征界。
% 模块二：Roe声速采用正热项加速度方差，避免总焓减大动能的消去误差。
% 模块三：全单向波取原物理通量，混合波取同一HLL平均态；不裁状态或换方法。
% 模块四：保存实际波速、取样分支和中间态；NaN标识根本未使用的中间态。
% 单Riemann正性不能替代NS+SA/黏性/几何源/隐式迭代的整体正性或物理验收。
assert((ischar(wavePolicy)&&isrow(wavePolicy))||(isstring(wavePolicy)&&isscalar(wavePolicy)&&~ismissing(wavePolicy)), ...
    'AxisFV:InteriorFluxPolicy','HLLE波速策略须为单一文本。');wavePolicy=char(wavePolicy);
assert(ismember(wavePolicy,{'davis','einfeldt'}),'AxisFV:InteriorFluxPolicy','只支持davis/einfeldt波速。');
gas=axisymmetricFVGas(gas);L=axisymmetricFVState(QL,gas);R=axisymmetricFVState(QR,gas);nf=size(QL,1);
assert(isequal(size(QR),[nf 4])&&isequal(size(normal),[nf 2])&&isnumeric(normal)&&isreal(normal)&& ...
    all(isfinite(normal),'all')&&all(abs(sum(normal.^2,2)-1)<1e-12), ...
    'AxisFV:InteriorFluxGeometry','两侧完整状态必须匹配有限单位法向。');
unL=sum(L(:,2:3).*normal,2);unR=sum(R(:,2:3).*normal,2);
if strcmp(wavePolicy,'einfeldt')
    a=sqrt(L(:,1));b=sqrt(R(:,1));den=a+b;wl=a./den;wr=b./den;
    roeVelocity=wl.*L(:,2:3)+wr.*R(:,2:3);
    variance=wl.*wr.*sum((L(:,2:3)-R(:,2:3)).^2,2);
    soundSquared=wl.*L(:,6).^2+wr.*R(:,6).^2+.5*(gas.gamma-1)*variance;
    assert(all(isfinite(soundSquared)&soundSquared>0),'AxisFV:HLLE','Roe特征声速非有限或非正，未裁值。');
    roeSound=sqrt(soundSquared);roeNormal=sum(roeVelocity.*normal,2);
    low=min(unL-L(:,6),roeNormal-roeSound);high=max(unR+R(:,6),roeNormal+roeSound);
else
    low=min(unL-L(:,6),unR-R(:,6));high=max(unL+L(:,6),unR+R(:,6));
    roeVelocity=nan(nf,2);roeSound=nan(nf,1);
end
assert(all(isfinite([low high]),'all')&&all(low<high),'AxisFV:HLLE','HLLE波速非有限或波序无效。');
FL=[L(:,1).*unL,QL(:,2:3).*unL+L(:,5).*normal,(QL(:,4)+L(:,5)).*unL];
FR=[R(:,1).*unR,QR(:,2:3).*unR+R(:,5).*normal,(QR(:,4)+R(:,5)).*unR];
flux=FL;right=high<=0;flux(right,:)=FR(right,:);middle=low<0&high>0;
average=nan(nf,4);internal=nan(nf,1);branch=ones(nf,1);branch(right)=-1;branch(middle)=0;
if any(middle)
    lo=low(middle);hi=high(middle);d=hi-lo;
    value=(hi.*QR(middle,:)-lo.*QL(middle,:)-(FR(middle,:)-FL(middle,:)))./d;
    energy=value(:,4)-sum(value(:,2:3).^2,2)./(2*value(:,1));
    assert(all(isfinite(value),'all')&&all(value(:,1)>0&energy>0), ...
        'AxisFV:HLLE','HLLE实际平均态rho/内能非正或非有限，没有裁值或回退其他通量。');
    average(middle,:)=value;internal(middle)=energy;
    flux(middle,:)=(hi.*FL(middle,:)-lo.*FR(middle,:)+lo.*hi.*(QR(middle,:)-QL(middle,:)))./d;
end
assert(all(isfinite(flux),'all'),'AxisFV:HLLE','HLLE实际通量非有限。');
audit=struct('kind','explicit-interior-HLLE-Euler-flux-v1','wavePolicy',wavePolicy, ...
    'waveSpeedsMPerS',[low high],'sampleBranch',branch,'middleMask',middle, ...
    'intermediateConservativeState',average,'intermediateInternalEnergyPerVolume',internal, ...
    'roeVelocityXR',roeVelocity,'roeSoundSpeed',roeSound,'fluxPerArea',flux, ...
    'intermediateStateMeaning','Riemann-solver average, not a cell/transport/wall state', ...
    'unusedIntermediateMarkedNaN',true,'stateClipped',false,'adaptiveFluxFallback',false, ...
    'generalViscousImplicitPositivityProved',false,'physicalSolution',false);
end
