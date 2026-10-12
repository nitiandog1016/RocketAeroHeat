function [faces,sourcePressure,diagnostics] = reconstructEulerEnthalpy(primitive,reconstruction,gamma)
%RECONSTRUCTEULERENTHALPY 重构 rho/u/v/H，并由同一组状态恢复压力。
% 用途：常总焓场在所有面仍保持常总焓，不把独立线性 p/rho 与速度
% 拼成不一致的总焓。所有输出均为 [rho,u,v,p]，求解器接口不变。
% 该方法只适用于定 gamma 理想气体；不重设单元 H，不强迫场内常总焓。
validateattributes(gamma,{'numeric'},{'real','finite','scalar','>',1});
validateattributes(primitive,{'numeric'},{'real','finite','size',[reconstruction.nx,reconstruction.nr,4]});
if any(primitive(:,:,[1 4])<=0,'all')
 error('RocketAeroHeat:EulerReconstruction:State','重构入口的密度和压力必须为正。');
end
rho=primitive(:,:,1);u=primitive(:,:,2);v=primitive(:,:,3);p=primitive(:,:,4);
H=gamma/(gamma-1)*p./rho+.5*(u.^2+v.^2);
base=cat(3,rho,u,v,H);
[faces,~,diagnostics,source]=reconstructEulerPrimitive(base,reconstruction,true);
points=[faces,{source}];theta=ones(size(rho));staticH=H-.5*(u.^2+v.^2);
if any(staticH<=0,'all')
 error('RocketAeroHeat:EulerReconstruction:Resolution','入口静焓低于当前浮点总焓表示的分辨能力。');
end

%% 1. 缩小重构增量，保持所有面与源点的静焓为正
% 沿 theta 的静焓为凹二次函数，theta=0 时等于正的单元静焓；
% 因而其正值区间从零起连续，可用二分找保守上界，无需裁剪压力。
for j=1:5
 d=points{j}-base;h=points{j}(:,:,4)-.5*(points{j}(:,:,2).^2+points{j}(:,:,3).^2);
 % 静焓是总焓减动能，必须留出与这两个大量相称的舍入裕量。
 % 目标仍小于入口静焓，不用压力下限修改任何单元状态。
 floorH=min(.5*staticH,max(1e-10*staticH,64*eps*abs(H)));bad=h<=floorH;
 if any(bad,'all')
  du=d(:,:,2);dv=d(:,:,3);dh=d(:,:,4);
  lo=zeros(nnz(bad),1);hi=ones(size(lo));a=staticH(bad);b=dh(bad)-u(bad).*du(bad)-v(bad).*dv(bad);c=.5*(du(bad).^2+dv(bad).^2);
  for k=1:48
   mid=.5*(lo+hi);ok=a+mid.*b-mid.^2.*c>floorH(bad);
   lo(ok)=mid(ok);hi(~ok)=mid(~ok);
  end
  theta(bad)=min(theta(bad),lo*(1-64*eps));
 end
end

%% 2. 由热力学关系恢复正压力；共同 theta 保持同一单元各面相容
for j=1:5
 state=base+theta.*(points{j}-base);
 pressure=state(:,:,1).*(state(:,:,4)-.5*(state(:,:,2).^2+state(:,:,3).^2))*(gamma-1)/gamma;
 if any(~isfinite(pressure) | pressure<=0,'all')
  error('RocketAeroHeat:EulerReconstruction:Enthalpy','静焓斜率限制后仍出现非正压力。');
 end
 state(:,:,4)=pressure;
 if j<=4,faces{j}=state;else,sourcePressure=pressure;end
end
diagnostics.enthalpyPositivityLimitedCellCount=nnz(theta<1);
diagnostics.method=[diagnostics.method '；rho/u/v/H 重构；静焓正值二分斜率限制；同状态恢复 p'];
end
