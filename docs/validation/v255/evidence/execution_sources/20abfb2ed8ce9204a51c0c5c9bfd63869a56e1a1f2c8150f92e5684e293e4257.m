function [dh,ds]=phaseAwareAirCaloricDrop(T0,drop)
%PHASEAWAREAIRCALORICDROP 原Cp表的稳定焓差和熵差，不用两个大绝对量相减。
% 模块一：绑定原表SHA；按T0向下缓存分段线性Cp的解析差量积分。
% 模块二：温降作为独立变量，log1p保证近驻点极小温降仍有精度。
% 返回h(T0)-h(T0-drop)、s0(T0)-s0(T0-drop)；不授p/T物理资格。
persistent origin d cp slopes h s
validateattributes(T0,{'numeric'},{'real','scalar','finite','>=',60,'<=',2000});
validateattributes(drop,{'numeric'},{'real','finite','nonempty','>=',0,'<=',T0-60});
if isempty(origin)||origin~=T0
    c=phaseAwareAirCaloric(T0); % 原入口验证同一固定SHA。
    file=fullfile(fileparts(mfilename('fullpath')),'data','thermo','phase_aware_dry_air_v1.json');
    table=jsondecode(fileread(file));below=table.temperatureK(table.temperatureK<T0);
    temps=[T0;flipud(below(:))];d=T0-temps;
    cp=interp1(table.temperatureK,table.cpJkgK,temps,'linear');cp(1)=c.cp;
    slopes=-diff(cp)./diff(d);dd=diff(d);tb=temps(1:end-1);
    h=[0;cumsum(dd.*(cp(1:end-1)+cp(2:end))/2)];
    s=[0;cumsum(-(cp(1:end-1)-slopes.*tb).*log1p(-dd./tb)+slopes.*dd)];
    origin=T0;
end
if T0==60,dh=zeros(size(drop));ds=dh;return;end
index=discretize(drop,d);base=reshape(d(index),size(drop));local=drop-base;
cb=reshape(cp(index),size(drop));grad=reshape(slopes(index),size(drop));tb=T0-base;
dh=reshape(h(index),size(drop))+local.*(cb-.5*grad.*local);
ds=reshape(s(index),size(drop))-(cb-grad.*tb).*log1p(-local./tb)+grad.*local;
end
