function [faces,sourcePressure,diagnostics,sourceState] = reconstructEulerPrimitive(primitive,reconstruction,limited)
%RECONSTRUCTEULERPRIMITIVE 受限线性重构 rho、u、v、p，返回四个面的原始量。
% 所有面采用 Venkatakrishnan 型平滑限制，包括物理边界，避免无约束外推。
% 边界/极值处可能退化为一阶，不声称全域二阶精度。所有面的 rho/p
% 通过缩小梯度保持正值，不裁剪单元守恒量。
% limited=false 仅用于独立解析场检查，不应在含激波计算中关闭限制器。
% 总焓重构器也调用此核，将第四正变量换成 H；第四输出保留完整源点
% 状态，随后由调用者恢复压力并另外检查静焓正值。
if nargin<3,limited=true;end
nx=reconstruction.nx;nr=reconstruction.nr;
validateattributes(primitive,{'numeric'},{'real','finite','size',[nx,nr,4]});
q=reshape(primitive,[],4);
if any(q(:,[1 4])<=0,'all'),error('RocketAeroHeat:EulerReconstruction:State','密度和压力必须为正。');end
neighbors=reconstruction.neighbors;n=size(q,1);gx=zeros(n,4);gr=gx;qmin=q;qmax=q;
for j=1:4
    neighbor=q(neighbors(:,j),:);
    if isfield(reconstruction,'ghost')
        selected=reconstruction.ghost(:,j);
        gnX=reconstruction.ghostNormalX(selected,j);gnR=reconstruction.ghostNormalR(selected,j);
        vn=neighbor(selected,2).*gnX+neighbor(selected,3).*gnR;
        neighbor(selected,2)=neighbor(selected,2)-2*vn.*gnX;
        neighbor(selected,3)=neighbor(selected,3)-2*vn.*gnR;
    end
    dq=neighbor-q;
    gx=gx+reconstruction.coefficientX(:,j).*dq;gr=gr+reconstruction.coefficientR(:,j).*dq;
    qmin=min(qmin,neighbor);qmax=max(qmax,neighbor);
end
theta=ones(n,4);
if limited
    epsilonSquared=reconstruction.smoothingLengthCubed.*max(abs(q),[],1).^2;
    for j=1:4
        delta=gx.*reconstruction.deltaX(:,j)+gr.*reconstruction.deltaR(:,j);
        negative=delta<0;
        upper=qmax-q;lower=qmin-q;
        bound=upper;bound(negative)=lower(negative);
        ratio=(bound.^2+2*bound.*delta+epsilonSquared)./ ...
            (bound.^2+bound.*delta+2*delta.^2+epsilonSquared);
        ratio(delta==0)=1;
        % 平滑参数按网格尺度三次方减小，抑制极小扰动的限制器反复切换。
        % 允许平滑尺度内的弱超调；rho/p 正值由下面独立的斜率限制保证。
        theta=min(theta,max(0,min(1,ratio)));
    end
end
% 边界外推和压力源点也要保持正值；限制梯度而非修改返回的 rho/p。
dx=[reconstruction.deltaX,reconstruction.sourceDeltaX];dr=[reconstruction.deltaR,reconstruction.sourceDeltaR];
for j=1:5
    delta=gx.*dx(:,j)+gr.*dr(:,j);
    for k=[1 4]
        negative=delta(:,k)<0;
        bound=(1-1e-12)*q(negative,k)./(-delta(negative,k));
        theta(negative,k)=min(theta(negative,k),bound);
    end
end
gx=theta.*gx;gr=theta.*gr;faces=cell(1,4);
for j=1:4
    faces{j}=reshape(q+gx.*reconstruction.deltaX(:,j)+gr.*reconstruction.deltaR(:,j),nx,nr,4);
end
sourcePressure=reshape(q(:,4)+gx(:,4).*reconstruction.sourceDeltaX+gr(:,4).*reconstruction.sourceDeltaR,nx,nr);
if nargout>3,sourceState=reshape(q+gx.*reconstruction.sourceDeltaX+gr.*reconstruction.sourceDeltaR,nx,nr,4);end
diagnostics.minimumLimiter=min(theta,[],'all');
diagnostics.limitedCellCount=nnz(any(theta<1-1e-10,2));
diagnostics.boundaryMethod=reconstruction.boundaryMethod;
diagnostics.method='Venkatakrishnan 型有理平滑限制；epsilon²=(h/L)³ max|q|²，h=sqrt(子午面积)；正值斜率限制';
end
