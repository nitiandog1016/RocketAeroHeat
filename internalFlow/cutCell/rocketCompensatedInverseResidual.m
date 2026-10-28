function [residual,errorBound]=rocketCompensatedInverseResidual(A,X)
%ROCKETCOMPENSATEDINVERSERESIDUAL 补偿计算原双精度输入的I-A*X及保守舍入界。
% 模块一：Dekker乘积误差分解与TwoSum求和，减少近相消逆乘积的测量误差。
% 模块二：只提高残差评价精度，不修改A、X或1e-8接受阈值。
% 模块三：超出安全分裂范围明确拒绝，误差界含二阶舍入与下溢项。
assert(isa(A,'double')&&isa(X,'double')&&isreal(A)&&isreal(X)&&isequal(size(A),size(X))&&size(A,1)==size(A,2)&&all(isfinite(A),'all')&&all(isfinite(X),'all'), ...
    'InverseResidual:Input','需要同阶有限双精度实方阵。');
N=size(A,1);splitter=134217729;maximumA=max(abs(A),[],'all');maximumX=max(abs(X),[],'all');
assert(maximumA<realmax/splitter&&maximumX<realmax/splitter&&maximumA<realmax/max(1,N*maximumX), ...
    'InverseResidual:Range','误差分解的乘积或分裂可能溢出，不能给出证书。');
total=zeros(N);compensation=total;
for k=1:N
    a=A(:,k);b=X(k,:);ta=splitter*a;tb=splitter*b;ahi=ta-(ta-a);alo=a-ahi;bhi=tb-(tb-b);blo=b-bhi;
    product=a.*b;productError=((ahi.*bhi-product)+ahi.*blo+alo.*bhi)+alo.*blo;
    next=total+product;part=next-total;sumError=(total-(next-part))+(product-part);
    compensation=compensation+(productError+sumError);total=next;
end
residual=(eye(N)-total)-compensation;gamma=N*eps/(1-N*eps);
magnitudeBound=norm(A,inf)*maximumX/(1-gamma);errorBound=8*gamma^2*magnitudeBound+8*eps*max(abs(residual),[],'all')+64*N*realmin;
assert(all(isfinite(residual),'all')&&isfinite(errorBound),'InverseResidual:Nonfinite','补偿残差或误差界非有限。');
end
