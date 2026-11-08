function state=phaseAwareAirState(mode,p,value)
%PHASEAWAREAIRSTATE 相态有界标准干空气，供显式选择的新混合热链使用。
% 模块用途：PT、Ph、Pe三个入口分别以压力/温度、压力/相对焓、压力/热内能
% 恢复同一Cp表的状态。60—2000K；气相/低密度联合域；SI单位。
% 旧NASA7、NIST、OTH及工程默认不变。整次输入非法即报错，不补值/裁温。

%% 1. 绑定生成数据及具名物性身份，缓存单调积分表
persistent table knots cpNodes slopes hNodes sNodes digest
if isempty(table)
    path=fullfile(fileparts(mfilename('fullpath')),'data','thermo','phase_aware_dry_air_v1.json');
    fid=fopen(path,'rb'); assert(fid>=0,'RAH:PhaseAir:File','物性表不可读。');
    closer=onCleanup(@()fclose(fid)); bytes=fread(fid,Inf,'*uint8');
    md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(bytes,'int8'));
    digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    assert(strcmp(digest,'7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'), ...
        'RAH:PhaseAir:SHA','物性表与已验证版本不一致。');
    candidate=jsondecode(native2unicode(bytes.','UTF-8'));
    assert(strcmp(candidate.model,'rah-phase-aware-dry-air-of13-v1') && ...
        strcmp(candidate.domain,'rah-lemmon-gas-low-density-60-2000-v1'),'RAH:PhaseAir:Identity','气体身份不一致。');
    knots=candidate.temperatureK(:); cpNodes=candidate.cpJkgK(:); R=candidate.gasConstantJkgK;
    assert(knots(1)==60 && knots(end)==2000 && all(diff(knots)>0) && ...
        all(isfinite(cpNodes)&cpNodes>R),'RAH:PhaseAir:Table','量热表无效。');
    slopes=diff(cpNodes)./diff(knots);
    hNodes=[0;cumsum((cpNodes(1:end-1)+cpNodes(2:end)).*diff(knots)/2)];
    sNodes=[0;cumsum((cpNodes(1:end-1)-slopes.*knots(1:end-1)).*log(knots(2:end)./knots(1:end-1))+slopes.*diff(knots))];
    anchor=find(knots==298.15); assert(isscalar(anchor),'RAH:PhaseAir:Anchor','参考节点缺失。');
    hNodes=hNodes-hNodes(anchor); sNodes=sNodes-sNodes(anchor); table=candidate;
end
mode=validatestring(mode,{'PT','Ph','Pe'});
validateattributes(p,{'numeric'},{'real','finite','nonempty','positive'});
validateattributes(value,{'numeric'},{'real','finite','nonempty'});
[p,value]=match(double(p),double(value)); R=table.gasConstantJkgK;

%% 2. 焓/内能有界反解，只读数学量热表，不用旧压力判定中间试探温度
if strcmp(mode,'PT')
    T=value;
else
    low=60*ones(size(value)); high=2000*ones(size(value));
    [~,flo]=thermal(low); [~,fhi]=thermal(high);
    if strcmp(mode,'Pe'), flo=flo-R*low; fhi=fhi-R*high; end
    assert(all(value(:)>=flo(:)&value(:)<=fhi(:)), ...
        'RAH:PhaseAir:Energy','真实目标焓/内能超出来源温区。');
    % 60次使区间达到双精度分辨率，保留零参考焓的严格绝对残差门槛。
    for k=1:60
        mid=low+(high-low)/2; [~,f]=thermal(mid);
        if strcmp(mode,'Pe'), f=f-R*mid; end
        below=f<value; low(below)=mid(below); high(~below)=mid(~below);
    end
    T=low+(high-low)/2; T(value==flo)=60; T(value==fhi)=2000;
end

%% 3. 最终p/T联合域及完整热状态；相对熵不是绝对化学熵
upper=pressureUpper(T);
assert(all(p(:)<=upper(:)),'RAH:PhaseAir:Domain','最终压力温度超出气相/低密度范围。');
[cp,h,s0]=thermal(T); rho=p./(R*T); e=h-R*T; s=s0-R*(log(p)-log(100000));
assert(all(isfinite(rho(:))&rho(:)>0&isfinite(s(:))),'RAH:PhaseAir:Numerics','状态超出浮点数值范围。');
if ~strcmp(mode,'PT')
    actual=h; if strcmp(mode,'Pe'), actual=e; end
    assert(all(abs(actual(:)-value(:))<=2e-10*max(1,abs(value(:)))), ...
        'RAH:PhaseAir:Inverse','有界反算残差未通过。');
end

%% 4. Lemmon2004低密度输运与身份；连续流资格仍需局部长度/梯度检查
x=log(T/103.3); omega=exp(.431-.4623*x+.08406*x.^2+.005341*x.^3-.00331*x.^4);
mu=2.66958e-8*sqrt(28.9586*T)./(.36^2*omega); tau=132.6312./T;
conductivity=.001308*(mu*1e6)+.001405*tau.^(-1.1)-.001036*tau.^(-.3);
assert(all(isfinite(mu(:))&mu(:)>0&isfinite(conductivity(:))&conductivity(:)>0), ...
    'RAH:PhaseAir:Transport','输运状态无效。');
state=struct('pressurePa',p,'temperatureK',T,'cpJkgK',cp,'cvJkgK',cp-R, ...
    'enthalpyJkg',h,'sensibleEnthalpyJkg',h,'internalEnergyJkg',e,'entropyJkgK',s, ...
    'densityKgm3',rho,'soundSpeedMps',sqrt(cp./(cp-R).*R.*T), ...
    'dynamicViscosityPas',mu,'conductivityWmK',conductivity,'Prandtl',cp.*mu./conductivity, ...
    'meanFreePathViscosityDefinitionM',mu./p.*sqrt(pi*R*T/2), ...
    'model',table.model,'domain',table.domain,'sourceSHA256',digest, ...
    'temperatureRangeK',[60 2000],'pressureCeilingPa',upper,'gasConstantJkgK',R, ...
    'energyReference','hf=0; hs(298.15K)=0; es=hs-RT', ...
    'entropyReference','relative s(298.15K,100000Pa)=0', ...
    'transportMethod','Lemmon2004 dilute Eq.2/5; not Sutherland or constant Pr', ...
    'continuumQualified',false,'physicalQualification',false,'engineeringHeatReady',false);

    function [cp,h,s]=thermal(t)
        % 表内解析积分，不借助外推或温度截断。
        assert(all(isfinite(t(:))&t(:)>=60&t(:)<=2000),'RAH:PhaseAir:Temperature','温度超出60—2000K。');
        index=discretize(t,knots); grad=slopes(index);
        % 索引列向量按输入形状还原，避免行列隐式广播。
        baseCp=reshape(cpNodes(index),size(t)); grad=reshape(grad,size(t));
        dt=t-reshape(knots(index),size(t)); baseT=reshape(knots(index),size(t));
        cp=baseCp+grad.*dt; h=reshape(hNodes(index),size(t))+baseCp.*dt+grad.*dt.^2/2;
        s=reshape(sNodes(index),size(t))+(baseCp-grad.*baseT).*log(t./baseT)+grad.*dt;
    end
end

function upper=pressureUpper(t)
% Lemmon2000露点式加事前0.1气相距离；高于maxcondentherm不调用露点式。
assert(all(isfinite(t(:))&t(:)>=60&t(:)<=2000),'RAH:PhaseAir:Temperature','温度超出60—2000K。');
r=t/100; upper=min(1e6,1e4*((r.*r).*r)); cold=t<=132.6312;
theta=1-t(cold)/132.6312;
dew=3785020*exp(132.6312./t(cold).*(-.1567266*sqrt(theta)-5.539635*theta+.7567212*theta.^2.5-3.514322*theta.^4));
upper(cold)=min(upper(cold),.1*dew);
end

function [a,b]=match(a,b)
% 仅同尺寸或标量扩展，不允许行列向量意外构成二维状态网格。
if isscalar(a), a=repmat(a,size(b));
elseif isscalar(b), b=repmat(b,size(a));
else, assert(isequal(size(a),size(b)),'RAH:PhaseAir:Shape','状态尺寸不一致。'); end
end
