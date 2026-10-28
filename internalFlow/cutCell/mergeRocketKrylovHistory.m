function [history,work]=mergeRocketKrylovHistory(previous,previousTotal,local,localTotal)
%MERGEROCKETKRYLOVHISTORY 合并接受状态历史，保留最后接受后的全部实际求解工作量。
% 模块一：接受状态的调用数保持其发生时刻，不把后来拒绝计算改写到旧行。
% 模块二：新段调用数以旧段实际总数偏移，包含旧段末尾没有新接受状态的调用。
% 模块三：总工作量独立求和；它不是迭代次数、物理时间或已解区域面积。
validateHistory(previous);validateHistory(local);
assert(isequal(sort(fieldnames(previous)),sort(fieldnames(local))), ...
    'KrylovHistory:Fields','前后历史字段必须一致。');
validateattributes(previousTotal,{'numeric'},{'scalar','real','finite','integer','nonnegative'});
validateattributes(localTotal,{'numeric'},{'scalar','real','finite','integer','nonnegative'});
assert(previousTotal>=previous.residualEvaluations(end)&&localTotal>=local.residualEvaluations(end), ...
    'KrylovHistory:Work','实际总调用数不能小于已记录状态时刻的调用数。');
assert(local.iteration(1)==0,'KrylovHistory:Iteration','本段历史必须从零开始。');
history=previous;priorAccepted=previous.iteration(end);
for name=fieldnames(local)'
    key=name{1};values=local.(key)(2:end);
    if strcmp(key,'iteration'),values=values+priorAccepted;end
    if strcmp(key,'residualEvaluations'),values=values+previousTotal;end
    history.(key)=[previous.(key),values];
end
work=struct('priorAcceptedIterations',priorAccepted,'priorActualResidualEvaluations',previousTotal, ...
    'localActualResidualEvaluations',localTotal,'totalActualResidualEvaluations',previousTotal+localTotal, ...
    'totalAcceptedIterations',history.iteration(end), ...
    'residualEvaluationsAfterLastAcceptedHistory',previousTotal+localTotal-history.residualEvaluations(end));
end
function validateHistory(history)
assert(isstruct(history)&&isscalar(history)&&all(isfield(history,{'iteration','residualEvaluations'})), ...
    'KrylovHistory:Fields','历史必须包含接受序号和实际调用数。');
n=numel(history.iteration);assert(n>0,'KrylovHistory:Shape','历史不能没有基线。');
for name=fieldnames(history)'
    value=history.(name{1});assert(isnumeric(value)&&isreal(value)&&isrow(value)&&numel(value)==n, ...
        'KrylovHistory:Shape','历史字段须为同长度实数行向量。');
end
assert(all(isfinite(history.iteration))&&all(history.iteration==fix(history.iteration))&& ...
    history.iteration(1)==0&&isequal(history.iteration,0:history.iteration(end)), ...
    'KrylovHistory:Iteration','接受序号必须从零连续递增。');
calls=history.residualEvaluations;
assert(all(isfinite(calls))&&all(calls>=0&calls==fix(calls))&&all(diff(calls)>=0), ...
    'KrylovHistory:Work','调用数必须为不递减的非负整数。');
end
