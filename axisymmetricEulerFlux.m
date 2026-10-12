function [flux,rate,fallbackCount] = axisymmetricEulerFlux(left,right,sx,sr,gamma,scheme)
%AXISYMMETRICEULERFLUX 一阶 HLL 及上风总焓能量输运的保守面通量。
% 输入左右守恒量 [rho,rho*u,rho*v,E]、带面积权重的二维法向。
% HLLenthalpy 是明确标注的混合试验方案：质量/动量采用标准 HLL，
% 能量使用同一质量通量乘以上风总焓。不会重设单元总焓或温度。
% 两种方案均使用唯一公共面通量；不承诺精确解析接触/剪切间断。
% 第三个返回量保留接口兼容；当前两种方案不自动切换，恒为零。
if nargin<6,scheme='HLLenthalpy';end
[rl,ul,vl,pl,al]=primitive(left,gamma);[rr,ur,vr,pr,ar]=primitive(right,gamma);
area=hypot(sx,sr);nx=sx./max(area,realmin);nr=sr./max(area,realmin);
ql=ul.*nx+vl.*nr;qr=ur.*nx+vr.*nr;
low=min(min(ql-al,qr-ar),0);high=max(max(ql+al,qr+ar),0);
%% 1. 标准 HLL 两波保守通量
fl=cat(3,rl.*ql,rl.*ql.*ul+pl.*nx,rl.*ql.*vl+pl.*nr,(left(:,:,4)+pl).*ql);
fr=cat(3,rr.*qr,rr.*qr.*ur+pr.*nx,rr.*qr.*vr+pr.*nr,(right(:,:,4)+pr).*qr);
flux=(high.*fl-low.*fr+low.*high.*(right-left))./(high-low);
%% 2. 试验方案只改变能量面通量，并保持逐面严格守恒
if strcmp(scheme,'HLLenthalpy')
    hl=(left(:,:,4)+pl)./rl;hr=(right(:,:,4)+pr)./rr;
    h=hr;positive=flux(:,:,1)>=0;h(positive)=hl(positive);
    flux(:,:,4)=flux(:,:,1).*h;
elseif ~strcmp(scheme,'HLL')
    error('RocketAeroHeat:Euler:FluxScheme','数值通量必须为 HLL 或 HLLenthalpy。');
end
flux=area.*flux;rate=area.*max(abs(low),abs(high));fallbackCount=0;
end
function [rho,u,v,p,a]=primitive(Q,g)
rho=Q(:,:,1);u=Q(:,:,2)./rho;v=Q(:,:,3)./rho;
p=(g-1)*(Q(:,:,4)-.5*rho.*(u.^2+v.^2));a=sqrt(g*p./rho);
end
