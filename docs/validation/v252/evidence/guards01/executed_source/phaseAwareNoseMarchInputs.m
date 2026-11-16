function [edges,stations]=phaseAwareNoseMarchInputs(e,angles,R)
%PHASEAWARENOSEMARCHINPUTS 复用已有同熵偏置球鼻外缘公式，生成非相似链的连续入口。
% 模块用途：原压力/总焓/驻点不改变；同32点Gauss积分xi，不重做激波或物性表。
% 只代表已声明的旋转球鼻工程近似，不能生成全后体或真实三维横流外缘。
assert(angles(1)==0&&all(diff(angles)>0)&&R>0);
states=edgeStates(angles);k=(1:31)';b=k./sqrt(4*k.^2-1);[Q,L]=eig(diag(b,1)+diag(b,-1));
[gx,order]=sort(diag(L));gw=Q(1,order).^2;gx=(gx+1)/2;
qa=angles(:)*gx(:)';quad=edgeStates(qa);
xi=R^3*angles(:).*sum(quad.densityKgm3.*quad.dynamicViscosityPas.*quad.U.*sin(qa).^2.*gw,2);
stations=repmat(struct('thetaRad',0,'arcM',0,'metricM',0,'xi',0,'beta',0, ...
    'pressureGradientPam',0,'velocityGradientPerSecond',0,'pressureCurvaturePam2',0),1,numel(angles));
edges=repmat(struct('pressurePa',0,'temperatureK',0,'speedMps',0,'totalEnthalpyJkg',0),1,numel(angles));
for j=1:numel(angles)
    th=angles(j);U=states.U(j);metric=R*sin(th);
    if j==1,du=e.stagnation.velocityGradientPerSecond;dp=0;beta=.5;
    else
        dp=-(e.stagnation.pressurePa-e.condition.pressurePa)*sin(2*th)/R;
        du=-dp/(states.densityKgm3(j)*U);
        beta=2*xi(j)*du/(states.densityKgm3(j)*states.dynamicViscosityPas(j)*U^2*metric^2);
    end
    stations(j)=struct('thetaRad',th,'arcM',R*th,'metricM',metric,'xi',xi(j),'beta',beta, ...
        'pressureGradientPam',dp,'velocityGradientPerSecond',du,'pressureCurvaturePam2',e.stagnation.pressureCurvaturePam2);
    edges(j)=struct('pressurePa',states.pressurePa(j),'temperatureK',states.temperatureK(j), ...
        'speedMps',U,'totalEnthalpyJkg',e.stagnation.enthalpyJkg);
end
    function s=edgeStates(th)
        p=e.condition.pressurePa+(e.stagnation.pressurePa-e.condition.pressurePa)*cos(th).^2;
        c0=phaseAwareAirCaloric(e.stagnation.temperatureK);target=c0.s0+c0.R*log(p/e.stagnation.pressurePa);
        lo=60*ones(size(th));hi=e.stagnation.temperatureK*ones(size(th));
        for it=1:50
            mid=(lo+hi)/2;c=phaseAwareAirCaloric(mid);below=c.s0<target;lo(below)=mid(below);hi(~below)=mid(~below);
        end
        T=(lo+hi)/2;T(th==0)=e.stagnation.temperatureK;s=phaseAwareAirState('PT',p,T,'linearTable');
        dh=e.stagnation.enthalpyJkg-s.enthalpyJkg;dh(th==0)=0;assert(all(dh(:)>=0));s.U=sqrt(2*dh);
    end
end
