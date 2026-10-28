function [flux,audit]=axisymmetricHLLEFluxAccurate(QL,QR,normal,gas,wavePolicy)
%AXISYMMETRICHLLEFLUXACCURATE 显式双double形成的二维四量HLLE，最终通量仍为double。
% 模块一：保留原输入/温域检查，从原Q及同一gamma/R/cv重建高精度声学量。
% 模块二：同一Davis/Einfeldt波速与HLL公式使用TwoSum/Dekker余项，非新物理模型。
% 模块三：只检查实际使用的混合波平均态正性；不裁态、不换波速或自动回退。
gas=axisymmetricFVGas(gas);axisymmetricFVState(QL,gas);axisymmetricFVState(QR,gas);
assert(isequal(size(QL),size(QR))&&isnumeric(normal)&&isreal(normal)&&isequal(size(normal),[size(QL,1),2])&& ...
    all(isfinite(normal),'all')&&all(abs(sum(normal.^2,2)-1)<1e-12),'AxisFV:InteriorFluxGeometry','需要匹配状态和单位法向。');
assert(((ischar(wavePolicy)&&isrow(wavePolicy))||(isstring(wavePolicy)&&isscalar(wavePolicy)))&& ...
    ismember(char(wavePolicy),{'davis','einfeldt'}),'AxisFV:InteriorFluxPolicy','只支持Davis/Einfeldt。');wavePolicy=char(wavePolicy);
L=primitive(QL,gas);R=primitive(QR,gas);n1=D(normal(:,1));n2=D(normal(:,2));
unL=add(mul(L.u,n1),mul(L.v,n2));unR=add(mul(R.u,n1),mul(R.v,n2));
if strcmp(wavePolicy,'einfeldt')
    al=rt(L.rho);ar=rt(R.rho);den=add(al,ar);wl=div(al,den);wr=div(ar,den);
    ru=add(mul(wl,L.u),mul(wr,R.u));rv=add(mul(wl,L.v),mul(wr,R.v));
    du=sub(L.u,R.u);dv=sub(L.v,R.v);variance=mul(mul(wl,wr),add(mul(du,du),mul(dv,dv)));
    sound2=add(add(mul(wl,mul(L.a,L.a)),mul(wr,mul(R.a,R.a))),mul(mul(D(.5),sub(D(gas.gamma),D(1))),variance));
    sound=rt(sound2);rn=add(mul(ru,n1),mul(rv,n2));
    low=choose(sub(unL,L.a),sub(rn,sound),true);high=choose(add(unR,R.a),add(rn,sound),false);
    roeVelocity=[value(ru),value(rv)];roeSound=value(sound);
else
    low=choose(sub(unL,L.a),sub(unR,R.a),true);high=choose(add(unL,L.a),add(unR,R.a),false);
    roeVelocity=NaN(size(QL,1),2);roeSound=NaN(size(QL,1),1);
end
assert(all(less(low,high)),'AxisFV:HLLE','波序错误。');
FL=physical(QL,L,unL,n1,n2);FR=physical(QR,R,unR,n1,n2);
middle=less(low,D(0))&less(D(0),high);right=~less(D(0),high);
average=NaN(size(QL));internal=NaN(size(QL,1),1);branch=ones(size(QL,1),1);branch(right)=-1;branch(middle)=0;
result=FL;result.hi(right,:)=FR.hi(right,:);result.lo(right,:)=FR.lo(right,:);
if any(middle)
    lo=rows(low,middle);hi=rows(high,middle);d=sub(hi,lo);fl=rows(FL,middle);fr=rows(FR,middle);
    ql=D(QL(middle,:));qr=D(QR(middle,:));
    formed=div(add(sub(mul(hi,fl),mul(lo,fr)),mul(mul(lo,hi),sub(qr,ql))),d);
    state=div(sub(sub(mul(hi,qr),mul(lo,ql)),sub(fr,fl)),d);
    kinetic=div(add(mul(column(state,2),column(state,2)),mul(column(state,3),column(state,3))),mul(D(2),column(state,1)));
    energy=sub(column(state,4),kinetic);
    assert(all(less(D(0),column(state,1))&less(D(0),energy)),'AxisFV:HLLE','实际平均态rho/内能非正，未裁值。');
    average(middle,:)=value(state);internal(middle)=value(energy);
    result.hi(middle,:)=formed.hi;result.lo(middle,:)=formed.lo;
end
flux=value(result);
audit=struct('kind','explicit-interior-HLLE-Euler-flux-two-double-v1','wavePolicy',wavePolicy, ...
    'waveSpeedsMPerS',[value(low),value(high)],'waveLowParts',[low.lo,high.lo], ...
    'sampleBranch',branch,'middleMask',middle,'intermediateConservativeState',average, ...
    'intermediateInternalEnergyPerVolume',internal,'roeVelocityXR',roeVelocity,'roeSoundSpeed',roeSound, ...
    'fluxPerArea',flux,'fluxFormationLowPart',result.lo,'arithmeticPolicy','two-double-Euler-formation-v1', ...
    'intermediateStateMeaning','Riemann-solver average, not a cell/transport/wall state', ...
    'unusedIntermediateMarkedNaN',true,'stateClipped',false,'adaptiveFluxFallback',false, ...
    'generalViscousImplicitPositivityProved',false,'physicalSolution',false);
end
function p=primitive(Q,gas)
p.rho=D(Q(:,1));p.u=div(D(Q(:,2)),p.rho);p.v=div(D(Q(:,3)),p.rho);
e=sub(div(D(Q(:,4)),p.rho),mul(D(.5),add(mul(p.u,p.u),mul(p.v,p.v))));T=div(e,D(gas.cv));
p.p=mul(mul(p.rho,D(gas.R)),T);p.a=rt(mul(mul(D(gas.gamma),D(gas.R)),T));
assert(all(less(D(0),p.p)&less(D(0),p.a)),'AxisFV:State','高精度EOS压力/声速非正。');
end
function f=physical(Q,p,un,n1,n2)
mass=add(mul(D(Q(:,2)),n1),mul(D(Q(:,3)),n2));
a=add(mul(D(Q(:,2)),un),mul(p.p,n1));b=add(mul(D(Q(:,3)),un),mul(p.p,n2));e=mul(add(D(Q(:,4)),p.p),un);
f=struct('hi',[mass.hi,a.hi,b.hi,e.hi],'lo',[mass.lo,a.lo,b.lo,e.lo]);
end
function out=choose(a,b,minimum)
mask=less(a,b);if ~minimum,mask=~mask;end
out=b;out.hi(mask)=a.hi(mask);out.lo(mask)=a.lo(mask);
end
function yes=less(a,b)
yes=(a.hi<b.hi)|(a.hi==b.hi&a.lo<b.lo);
end
function out=rows(a,i),out=struct('hi',a.hi(i,:),'lo',a.lo(i,:));end
function out=column(a,j),out=struct('hi',a.hi(:,j),'lo',a.lo(:,j));end
function a=D(x),a=axisymmetricDoubleDouble('from',x);end
function a=value(x),a=axisymmetricDoubleDouble('value',x);end
function a=add(x,y),a=axisymmetricDoubleDouble('add',x,y);end
function a=sub(x,y),a=axisymmetricDoubleDouble('sub',x,y);end
function a=mul(x,y),a=axisymmetricDoubleDouble('mul',x,y);end
function a=div(x,y),a=axisymmetricDoubleDouble('div',x,y);end
function a=rt(x),a=axisymmetricDoubleDouble('sqrt',x);end
