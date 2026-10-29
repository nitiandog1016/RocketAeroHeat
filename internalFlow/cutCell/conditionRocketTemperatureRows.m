function [working,report]=conditionRocketTemperatureRows(H)
%CONDITIONROCKETTEMPERATUREROWS 对全部不同齐次行作可逆、二进制精确的表示预处理。
% 模块一：仅在支撑完全相同且逐项满足Sterbenz条件时，用当前行减该支撑首行。
% 模块二：差值按2的整数次幂缩放；逐值核对可逆，保留全部条件，不观察秩或删行。
% 模块三：原H由调用方保留；本函数只改善投影表示，不修改物理方程或原温域门槛。
assert(issparse(H)&&isreal(H)&&all(isfinite(nonzeros(H))), ...
    'TemperatureConditioning:Input','需要原有限实数稀疏条件。');
q=size(H,1);Ht=H';workingTranspose=Ht;supports=cell(q,1);values=cell(q,1);
changed=zeros(0,1);references=changed;powers=changed;
for j=1:q
    [indices,~,x]=find(Ht(:,j));supports{j}=indices;values{j}=x;
    k=0;
    for earlier=1:j-1
        if isequal(indices,supports{earlier}),k=earlier;break;end
    end
    if k==0,continue;end
    y=values{k};
    % 同号且一半至两倍以内，二进制浮点减法由Sterbenz引理保证精确。
    if ~all(sign(x)==sign(y)&abs(x)>=.5*abs(y)&abs(x)<=2*abs(y)),continue;end
    difference=x-y;
    assert(any(difference~=0),'TemperatureConditioning:Duplicates', ...
        '精确相同行应先由原重复映射表示，不能在本函数删去。');
    [~,exponent]=log2(max(abs(difference)));power=-exponent;
    scaled=pow2(difference,power);
    assert(all(isfinite(scaled))&&isequal(pow2(scaled,-power),difference)&& ...
        isequal(y+pow2(scaled,-power),x), ...
        'TemperatureConditioning:Exact','二进制差值或逆向重构不精确，拒绝而不近似合并。');
    workingTranspose(:,j)=sparse(indices,ones(numel(indices),1),scaled,size(H,2),1);
    changed(end+1,1)=j;references(end+1,1)=k;powers(end+1,1)=power; %#ok<AGROW>
end
working=workingTranspose';
report=struct('kind','exact-same-support-Sterbenz-difference-and-dyadic-row-preconditioning', ...
    'inputDistinctRowCount',q,'outputDistinctConditionCount',q, ...
    'transformedRowIndices',changed,'referenceRowIndices',references,'powerTwoExponents',powers, ...
    'originalHChanged',false,'allOriginalConditionsRetained',true, ...
    'rankBasedRowSelection',false,'rankTruncationUsed',false, ...
    'exactDifferenceBySterbenzVerified',true,'binaryInverseReconstructionVerified',true, ...
    'additionalArrayBoundBytes',64*nnz(H)+64*q*q+65536);
end
