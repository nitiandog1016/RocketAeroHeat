function result = marchHeadBoundaryLayerSegment(theta0,H0,U0,U1,nu0,nu1,ds)
%MARCHHEADBOUNDARYLAYERSEGMENT 推进一段准二维 Head 湍流边界层。
% 模块一：用动量厚度 theta 和形状因子 H=delta*/theta 表示状态。
% 动量式 theta'=Cf/2-(H+2)*theta*U'/U。
% Head 卷吸式 (U*theta*H1)'=U*E，E=0.0306/(H1-3)^0.6169。
% H1(H) 使用 Head 的两段经验闭合，Cf 用 Ludwieg-Tillmann。
% 原拟合在 H=1.6 有 H1≈0.02255 的接缝；高 H 分支系数由 1.5501
% 小幅归一化，使 H1 连续，避免凭空产生卷吸通量跳变；指数不变。
% 模块二：段内 U、参考运动黏度线性插值；RK4 步长加倍误差控制。
% 模块三：H>=2.4 仅作疑似准二维分离；域外/数值失败单独返回，
% 不输出附着流 Cf。该式原为不可压缩工程闭合，参考黏度修正并非
% 完整可压缩变换，不能求解横流、分离泡、再附或激波干扰。
% status: 0=成功，1=疑似分离，2=Re_theta 相关式域外，3=数值/形状域外。

values=[theta0 H0 U0 U1 nu0 nu1 ds];
validateattributes(values,{'numeric'},{'real','finite','positive'});
result=struct('theta',NaN,'shapeFactor',NaN,'Cf',NaN, ...
    'status',uint8(3),'fraction',0,'acceptedSteps',0, ...
    'reynoldsRange',[U0*theta0/nu0 U0*theta0/nu0]);
highCoefficient=0.8234/(1.6-1.1)^1.287*(1.6-0.6778)^3.064;
y=[theta0;H0]; x=0; step=ds;
[~,status]=rhs(0,y);
if status~=0, result.status=status; return; end
for attempt=1:2000
    step=min(step,ds-x);
    [full,s1]=rk4(x,y,step);
    [half,s2,r2]=rk4(x,y,step/2);
    s3=s2; fine=half;
    r3=r2;
    if s2==0, [fine,s3,r3]=rk4(x+step/2,half,step/2); end
    bad=max([s1 s2 s3]);
    if bad~=0
        if step<=max(ds*1e-7,eps(ds)*32)
            result.status=bad; result.fraction=x/ds; return;
        end
        step=step/2; continue;
    end
    % H 的绝对容差和 theta 的相对容差；不靠裁剪负厚度维持运行。
    scale=[max(abs(fine(1)),1e-12)*2e-6;2e-6];
    err=max(abs(fine-full)./(15*scale));
    if err<=1
        y=fine; x=x+step;
        result.reynoldsRange=[min([result.reynoldsRange(1) r2(1) r3(1)]) ...
            max([result.reynoldsRange(2) r2(2) r3(2)])];
        result.acceptedSteps=result.acceptedSteps+1;
        if ds-x<=max(ds*1e-12,eps(ds)*32)
            [~,status,cf]=rhs(ds,y);
            result.status=status; result.fraction=1;
            if status==0
                result.theta=y(1); result.shapeFactor=y(2); result.Cf=cf;
            end
            return;
        end
    end
    step=step*min(2,max(0.2,0.9*max(err,1e-8)^(-0.2)));
    if step<=eps(ds)*32, break; end
end
result.fraction=x/ds;

    function [dy,status,cf,Re]=rhs(s,state)
        dy=[NaN;NaN]; cf=NaN; Re=NaN; status=uint8(3);
        th=state(1); H=state(2);
        if any(~isfinite(state)) || th<=0 || H<=1.1, return; end
        if H>=2.4, status=uint8(1); return; end
        U=U0+(U1-U0)*(s/ds); nu=nu0+(nu1-nu0)*(s/ds);
        Re=U*th/nu;
        if Re<300 || Re>1e6, status=uint8(2); return; end
        if H<=1.6
            H1=3.3+0.8234/(H-1.1)^1.287;
            dH1=-1.287*0.8234/(H-1.1)^2.287;
        else
            H1=3.3+highCoefficient/(H-0.6778)^3.064;
            dH1=-3.064*highCoefficient/(H-0.6778)^4.064;
        end
        cf=0.246*10^(-0.678*H)*Re^(-0.268);
        E=0.0306/(H1-3)^0.6169;
        du=(U1-U0)/ds;
        dy(1)=cf/2-(H+2)*th*du/U;
        dy(2)=(E-H1*cf/2+H1*(H+1)*th*du/U)/(th*dH1);
        status=uint8(0);
    end

    function [value,status,range]=rk4(s,state,h)
        value=state; range=[NaN NaN];
        [a,status,~,ra]=rhs(s,state); if status~=0, return; end
        [b,status,~,rb]=rhs(s+h/2,state+h*a/2); if status~=0, return; end
        [c,status,~,rc]=rhs(s+h/2,state+h*b/2); if status~=0, return; end
        [d,status,~,rd]=rhs(s+h,state+h*c); if status~=0, return; end
        value=state+h*(a+2*b+2*c+d)/6;
        [~,status,~,re]=rhs(s+h,value);
        range=[min([ra rb rc rd re]) max([ra rb rc rd re])];
    end
end
