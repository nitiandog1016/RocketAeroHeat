function [flux,audit]=cutCellExactEulerFlux(QL,QR,normal,gas,options)
%CUTCELLEXACTEULERFLUX 定比热理想气体精确Riemann波系在xi=0的批量三维Godunov通量。
% 模块一：三维速度沿真实单位法向分解，压力曲线解两侧激波/稀疏波匹配。
% 模块二：保留严格正根括区，用有保护Newton迭代，不将负星压裁成正数。
% 模块三：真实波扇/星区采样及选中侧切向速度重建全部五量通量。
% 模块四：真空、根未收敛或采样温度越界明确拒绝；不是完整NS/热流或物理资格。
% 本文件为新数值候选，只有对应测试与实际范围通过后才可用于研究；原HLLC不回退到此。
if nargin<5,options=struct;end
defaults=struct('pressureRelativeTolerance',2e-13,'residualRelativeTolerance',2e-13, ...
    'maximumRootIterations',100,'maximumBracketExpansions',60);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'CutExactEuler:Options','未知精确Riemann选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
for key={'pressureRelativeTolerance','residualRelativeTolerance'},validateattributes(options.(key{1}),{'numeric'},{'real','finite','positive','scalar','<',1});end
for key={'maximumRootIterations','maximumBracketExpansions'},validateattributes(options.(key{1}),{'numeric'},{'real','finite','positive','integer','scalar'});end
N=size(QL,1);assert(N>0&&isequal(size(QL),[N 5])&&isequal(size(QR),[N 5])&&isequal(size(normal),[N 3])&& ...
    all(isfinite(normal),'all')&&all(abs(sum(normal.^2,2)-1)<1e-12),'CutExactEuler:Geometry','五量尺寸或真实单位法向不符。');
assert(isstruct(gas)&&all(isfield(gas,{'gamma','R','temperatureRange'}))&&isfinite(gas.gamma)&&gas.gamma>1&& ...
    isfinite(gas.R)&&gas.R>0&&numel(gas.temperatureRange)==2&&all(isfinite(gas.temperatureRange))&& ...
    gas.temperatureRange(1)>0&&gas.temperatureRange(2)>gas.temperatureRange(1),'CutExactEuler:Gas','需要明确常gamma及温域。');
g=gas.gamma;gm=g-1;[rhoL,vL,pL,TL]=state(QL);[rhoR,vR,pR,TR]=state(QR);
uL=sum(vL.*normal,2);uR=sum(vR.*normal,2);aL=sqrt(g*pL./rhoL);aR=sqrt(g*pR./rhoR);
vacuumMargin=2*(aL+aR)/gm-(uR-uL);
assert(all(vacuumMargin>0),'CutExactEuler:Vacuum','输入形成真实真空Riemann问题；不补造正压力。');

%% 模块一：以每面最大输入压力归一化根变量，避免任意有量纲压力下限
pressureScale=max(pL,pR);lo=zeros(N,1);hi=ones(N,1);expansions=zeros(N,1);
for k=1:options.maximumBracketExpansions
    f=pressureFunction(hi.*pressureScale);grow=f<0;
    if ~any(grow),break;end
    hi(grow)=2*hi(grow);expansions(grow)=expansions(grow)+1;
    assert(all(isfinite(hi)),'CutExactEuler:Bracket','压力括区溢出，不裁压力。');
end
assert(all(pressureFunction(hi.*pressureScale)>=0),'CutExactEuler:Bracket','正压上界未建立。');
x=sqrt(pL./pressureScale).*sqrt(pR./pressureScale);done=false(N,1);iterations=zeros(N,1);
velocityScale=aL+aR+abs(uR-uL);finalResidual=nan(N,1);lastStep=nan(N,1);
for k=1:options.maximumRootIterations
    [f,derivative]=pressureFunction(x.*pressureScale);dx=f./(derivative.*pressureScale);
    assert(all(isfinite(dx))&&all(derivative>0&isfinite(derivative)),'CutExactEuler:Root','压力曲线导数退化。');
    newly=~done&abs(f)<=options.residualRelativeTolerance*velocityScale& ...
        abs(dx)<=options.pressureRelativeTolerance*max(abs(x),eps);
    iterations(~done)=k;done=done|newly;finalResidual=f;lastStep=dx;
    if all(done),break;end
    low=~done&f<0;high=~done&f>=0;lo(low)=x(low);hi(high)=x(high);
    proposal=x-dx;unsafe=~isfinite(proposal)|proposal<=lo|proposal>=hi;
    proposal(unsafe)=.5*(lo(unsafe)+hi(unsafe));x(~done)=proposal(~done);
    assert(all(x>0&isfinite(x)),'CutExactEuler:Root','正根不能表示，拒绝而不加epsilon。');
end
assert(all(done),'CutExactEuler:Root','精确Riemann压力根未达到原指定容差。');
pStar=x.*pressureScale;fL=curve(pStar,pL,rhoL,aL);fR=curve(pStar,pR,rhoR,aR);
uStar=.5*(uL+uR+fR-fL);
assert(all(pStar>0&isfinite(pStar))&&all(isfinite(uStar)),'CutExactEuler:State','星状态非有限或压力非正。');

%% 模块二：按真实三波结构采样xi=0，不将星状态一律当作界面状态
left=uStar>=0;right=~left;rho=zeros(N,1);un=rho;pressure=rho;region=zeros(N,1);
rhoStarL=starDensity(rhoL,pL);rhoStarR=starDensity(rhoR,pR);
leftShock=left&pStar>pL;leftRare=left&~leftShock;rightShock=right&pStar>pR;rightRare=right&~rightShock;
waveL=uL-aL.*sqrt((g+1)/(2*g)*pStar./pL+gm/(2*g));
waveR=uR+aR.*sqrt((g+1)/(2*g)*pStar./pR+gm/(2*g));
headL=uL-aL;tailL=uStar-aL.*(pStar./pL).^(gm/(2*g));
headR=uR+aR;tailR=uStar+aR.*(pStar./pR).^(gm/(2*g));
set(leftShock&waveL>=0,rhoL,uL,pL,1);set(leftShock&waveL<0,rhoStarL,uStar,pStar,2);
set(leftRare&headL>=0,rhoL,uL,pL,1);set(leftRare&headL<0&tailL<=0,rhoStarL,uStar,pStar,2);
fan=leftRare&headL<0&tailL>0;sound=2/(g+1)*(aL+.5*gm*uL);
assert(all(sound(fan)>0),'CutExactEuler:State','左波扇采样声速非正。');
rhoFan=zeros(N,1);pFan=rhoFan;rhoFan(fan)=rhoL(fan).*(sound(fan)./aL(fan)).^(2/gm);pFan(fan)=pL(fan).*(sound(fan)./aL(fan)).^(2*g/gm);
set(fan,rhoFan,sound,pFan,3);
set(rightShock&waveR<=0,rhoR,uR,pR,4);set(rightShock&waveR>0,rhoStarR,uStar,pStar,5);
set(rightRare&headR<=0,rhoR,uR,pR,4);set(rightRare&headR>0&tailR>=0,rhoStarR,uStar,pStar,5);
fan=rightRare&headR>0&tailR<0;sound=2/(g+1)*(aR-.5*gm*uR);
assert(all(sound(fan)>0),'CutExactEuler:State','右波扇采样声速非正。');
rhoFan=zeros(N,1);pFan=rhoFan;rhoFan(fan)=rhoR(fan).*(sound(fan)./aR(fan)).^(2/gm);pFan(fan)=pR(fan).*(sound(fan)./aR(fan)).^(2*g/gm);
set(fan,rhoFan,-sound,pFan,6);
assert(all(region>0&rho>0&pressure>0),'CutExactEuler:State','波系采样不完整或正性失效。');
tangent=vL-uL.*normal;if any(right),tangent(right,:)=vR(right,:)-uR(right).*normal(right,:);end
velocity=tangent+un.*normal;temperature=pressure./(rho*gas.R);
bad=~isfinite(temperature)|temperature<gas.temperatureRange(1)|temperature>gas.temperatureRange(2);
if any(bad)
    first=find(bad,1);error('CutExactEuler:Temperature', ...
        '精确界面采样越过原物性域：批次面%d，T=%.17g K（整批%.17g—%.17gK）；不裁温或回退。', ...
        first,temperature(first),min(temperature),max(temperature));
end
energy=pressure/gm+.5*rho.*sum(velocity.^2,2);
flux=[rho.*un,rho.*velocity.*un+pressure.*normal,(energy+pressure).*un];
assert(all(isfinite(flux),'all'),'CutExactEuler:State','精确Godunov通量非有限。');
audit=struct('kind','exact-constant-gamma-Riemann-xi0-Godunov-three-dimensional-flux', ...
    'pressureStarPa',pStar,'normalVelocityStarMs',uStar,'samplePrimitive',[rho velocity temperature], ...
    'samplePressurePa',pressure,'sampleRegion',region,'pressureRootResidualMs',finalResidual, ...
    'normalizedPressureRootStep',lastStep,'pressureRootIterations',iterations,'bracketExpansions',expansions, ...
    'vacuumVelocityMarginMs',vacuumMargin,'inputTemperatureLeftK',TL,'inputTemperatureRightK',TR, ...
    'pressureClipped',false,'temperatureClipped',false,'automaticFallbackUsed',false,'physicalValidation',false, ...
    'regionNames',{{'left-state','left-star','left-rarefaction-fan','right-state','right-star','right-rarefaction-fan'}});
    function [rho,u,p,T]=state(q)
        assert(isreal(q)&&all(isfinite(q),'all')&&all(q(:,1)>0),'CutExactEuler:State','输入保守态/密度非法。');
        rho=q(:,1);u=q(:,2:4)./rho;internal=q(:,5)-.5*rho.*sum(u.^2,2);
        assert(all(internal>0&isfinite(internal)),'CutExactEuler:InternalEnergy','输入内能非正。');
        p=gm*internal;T=internal./(rho*gas.R/gm);
        assert(all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)),'CutExactEuler:Temperature','输入温度越过原模型域。');
    end
    function [value,derivative]=pressureFunction(p)
        [fl,dl]=curve(p,pL,rhoL,aL);[fr,dr]=curve(p,pR,rhoR,aR);value=fl+fr+uR-uL;derivative=dl+dr;
    end
    function [value,derivative]=curve(p,pk,rhok,ak)
        shock=p>pk;rare=~shock;value=zeros(N,1);derivative=value;
        if any(shock)
            A=2./((g+1)*rhok(shock));B=gm/(g+1)*pk(shock);factor=sqrt(A./(p(shock)+B));
            value(shock)=(p(shock)-pk(shock)).*factor;
            derivative(shock)=factor.*(1-.5*(p(shock)-pk(shock))./(p(shock)+B));
        end
        if any(rare)
            logRatio=log(p(rare)./pk(rare));value(rare)=2*ak(rare)/gm.*expm1(gm/(2*g)*logRatio);
            derivative(rare)=exp(-(g+1)/(2*g)*logRatio)./(rhok(rare).*ak(rare));
        end
    end
    function value=starDensity(rhok,pk)
        shock=pStar>pk;ratio=pStar./pk;value=rhok.*ratio.^(1/g);
        value(shock)=rhok(shock).*(ratio(shock)+gm/(g+1))./(gm/(g+1)*ratio(shock)+1);
    end
    function set(mask,rhoValue,uValue,pValue,label)
        rho(mask)=rhoValue(mask);un(mask)=uValue(mask);pressure(mask)=pValue(mask);region(mask)=label;
    end
end
