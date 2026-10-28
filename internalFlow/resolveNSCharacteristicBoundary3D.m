function [face,info]=resolveNSCharacteristicBoundary3D(cellState,b,n,gas)
%RESOLVENSCHARACTERISTICBOUNDARY3D 定比热局部法向开放特征边界草稿。
% 模块用途：由内部[rho u v w T]和明确的外部储库，构造开放边界面状态。
% n为指向计算域外的单位法向。所有坐标/变量采用SI；不裁温或翻转速度。
% b.mode='riemannFarfield'：入射J-由reservoirPrimitive给定，出射J+来自内部；
% 熵/切向速度按最终面流向选择。法向速度恰为零时显式保留内部熵/切向速度。
% b.mode='pressureOutlet'：出流只给pressure，保持内部J+/熵/切向速度；
% 回流必须由reservoirPrimitive给入射J-/熵/切向速度，不同时强加出口压力。
% 外推超声速出流；超声速入流须外部储库全部特征向内。声速转换明确拒绝。
% 来源：NASA/TM-2011-217181 §2.2式8--17的法向Riemann不变量；本实现
% 源报告印刷(14)声速系数与(8)矛盾；这里由(8)反解为(gamma-1)/4，
% 并用均匀状态校核，而不照抄该印刷项。没有FUN3D隐式Jacobian、
% 完整非反射NSCBC或熵稳定保证。压力出口组合由
% 同一特征关系推导，不使用该报告§2.5强制正出流速度的做法。
% 本函数只构造数值边界；给定p/储库不是已求外部气动环境。黏性边界另行指定。

%% 1. 校验显式气体、内部状态、外部储库与边界法向
validateattributes(cellState,{'numeric'},{'real','finite','numel',5});
validateattributes(n,{'numeric'},{'real','finite','numel',3});
n=double(n(:)');assert(abs(norm(n)-1)<1e-12, ...
    'RocketAeroHeat:NSDraft:CharacteristicNormal','开放边界法向须单位长度。');
validateattributes(gas.gamma,{'numeric'},{'real','finite','scalar','>',1});
validateattributes(gas.R,{'numeric'},{'real','finite','scalar','positive'});
range=double(gas.temperatureRange(:)');
validateattributes(range,{'double'},{'finite','positive','numel',2});
assert(range(2)>range(1),'RocketAeroHeat:NSDraft:CharacteristicRange','温度范围须递增。');
assert(isfield(b,'mode') && isfield(b,'reservoirPrimitive'), ...
    'RocketAeroHeat:NSDraft:CharacteristicInput','须明确mode及回流/远场储库，不隐式补来流。');
mode=validatestring(b.mode,{'riemannFarfield','pressureOutlet'});
c=state(double(cellState(:)'));
validateattributes(b.reservoirPrimitive,{'numeric'},{'real','finite','numel',5});
r=state(double(b.reservoirPrimitive(:)'));
if strcmp(mode,'pressureOutlet')
    assert(isfield(b,'pressure'),'RocketAeroHeat:NSDraft:CharacteristicInput','压力出口须给定pressure。');
    validateattributes(b.pressure,{'numeric'},{'real','finite','scalar','positive'});
end
uc=dot(c(2:4),n);ur=dot(r(2:4),n);
ac=sqrt(gas.gamma*gas.R*c(5));ar=sqrt(gas.gamma*gas.R*r(5));
pc=c(1)*gas.R*c(5);pr=r(1)*gas.R*r(5);
Kc=pc/c(1)^gas.gamma;Kr=pr/r(1)^gas.gamma;
Jplus=uc+2*ac/(gas.gamma-1);Jminus=ur-2*ar/(gas.gamma-1);
vc=c(2:4)-uc*n;vr=r(2:4)-ur*n;
margin=128*eps(max([1,ac,abs(uc)]));
assert(abs(abs(uc)-ac)>margin,'RocketAeroHeat:NSDraft:CharacteristicSonic', ...
    '内部状态处在法向声速转换；当前局部边界未定义转换分支。');
pressureInputApplied=false;reclassified=false;zeroTangentialPolicy='not-used';

%% 2. 按入射特征数构造边界：五入、四入、一入或零入
if uc>ac
    face=c;branch='supersonic-outflow';
elseif uc<-ac
    assert(ur<-ar,'RocketAeroHeat:NSDraft:CharacteristicRegime', ...
        '内部为超声速入流，但显式储库不能提供全部入射特征。');
    face=r;branch='supersonic-inflow';
elseif strcmp(mode,'riemannFarfield')
    [ub,ab]=acoustic();
    if ub<0
        face=fromAcoustic(ub,ab,Kr,vr);branch='reservoir-backflow';
    else
        face=fromAcoustic(ub,ab,Kc,vc);branch='subsonic-outflow';
        if ub==0,zeroTangentialPolicy='interior-entropy-and-tangential-velocity';end
    end
else
    % 先用内部速度选择出口/回流分支，再按实际边界面速度检验。
    % 若变号，只允许使用已经明确给出的另一组特征输入重新构造一次。
    if uc>=0
        trial=pressureTrace();ut=dot(trial(2:4),n);
        if ut>=0
            face=trial;branch='subsonic-outflow';pressureInputApplied=true;
        else
            [ub,ab]=acoustic();face=fromAcoustic(ub,ab,Kr,vr);
            assert(ub<0,'RocketAeroHeat:NSDraft:CharacteristicAmbiguous', ...
                '压力出口面变成回流，但储库回流分支不自洽；不强制翻转速度。');
            branch='reservoir-backflow';reclassified=true;
        end
    else
        [ub,ab]=acoustic();trial=fromAcoustic(ub,ab,Kr,vr);
        if ub<0
            face=trial;branch='reservoir-backflow';
        else
            face=pressureTrace();ub=dot(face(2:4),n);
            assert(ub>=0,'RocketAeroHeat:NSDraft:CharacteristicAmbiguous', ...
                '储库回流面变为出流，但压力出口分支不自洽。');
            branch='subsonic-outflow';pressureInputApplied=true;reclassified=true;
        end
    end
    if dot(face(2:4),n)==0,zeroTangentialPolicy='interior-entropy-and-tangential-velocity';end
end

%% 3. 核对最终面速度/声速、来源与数值黏性边界身份
face=state(face);un=dot(face(2:4),n);a=sqrt(gas.gamma*gas.R*face(5));
if ~startsWith(branch,'supersonic')
    assert(abs(un)<a-128*eps(max([1,a,abs(un)])), ...
        'RocketAeroHeat:NSDraft:CharacteristicRegime', ...
        '局部亚声速特征构造得到声速/超声速面；当前边界拒绝跨分支外推。');
end
incoming=double(un-a<0)+3*double(un<0)+double(un+a<0);
info=struct('branch',branch,'mode',mode,'normal',n,'normalVelocity',un, ...
    'normalMach',un/a,'incomingCharacteristics',incoming, ...
    'pressureInputApplied',pressureInputApplied,'pressureInput',NaN, ...
    'reclassifiedByFaceVelocity',reclassified,'zeroTangentialPolicy',zeroTangentialPolicy, ...
    'entropyParameter',face(1)*gas.R*face(5)/face(1)^gas.gamma, ...
    'Jplus',un+2*a/(gas.gamma-1),'Jminus',un-2*a/(gas.gamma-1), ...
    'reservoirPrimitiveInput',r,'isInflow',un<0, ...
    'method','常gamma局部法向Riemann特征边界；明确储库/出口压力', ...
    'scope','数值开放边界；无完整非反射/熵稳定保证，需域长度与边界敏感性检查', ...
    'temperatureStatus','由边界特征与显式输入重构；不是区域壁温解');
if strcmp(mode,'pressureOutlet'),info.pressureInput=b.pressure;end
if un<0
    info.viscousNormalGradient='characteristic-face-to-cell-half-spacing';
else
    info.viscousNormalGradient='prescribed-zero-normal-gradient-u-T';
end

    function p=state(p)
        assert(numel(p)==5 && isreal(p) && all(isfinite(p)) && p(1)>0 && ...
            p(5)>=range(1) && p(5)<=range(2), ...
            'RocketAeroHeat:NSDraft:CharacteristicState', ...
            '内部/储库/特征面密度或温度无效；不裁状态通过。');
    end
    function [u,a]=acoustic()
        u=.5*(Jplus+Jminus);a=(gas.gamma-1)*.25*(Jplus-Jminus);
        assert(isfinite(a) && a>0,'RocketAeroHeat:NSDraft:CharacteristicSound', ...
            '入射/出射声学不变量产生非正声速。');
    end
    function p=fromAcoustic(u,a,K,vt)
        rho=exp(log(a*a/(gas.gamma*K))/(gas.gamma-1));
        T=a*a/(gas.gamma*gas.R);p=[rho,vt+u*n,T];
    end
    function p=pressureTrace()
        rho=exp((log(b.pressure)-log(Kc))/gas.gamma);
        T=b.pressure/(rho*gas.R);a=sqrt(gas.gamma*gas.R*T);
        u=Jplus-2*a/(gas.gamma-1);p=state([rho,vc+u*n,T]);
    end
end
