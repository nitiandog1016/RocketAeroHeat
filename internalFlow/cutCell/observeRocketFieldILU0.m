function [value,observation]=observeRocketFieldILU0(factor,rhs)
%OBSERVEROCKETFIELDILU0 重放固定块回代，在首个非有限中间量处保存真实行证据。
% 模块一：有限前缀与原solveRocketFieldILU0算术顺序一致，不补主元、不裁值。
% 模块二：观察到溢出立即标记整个输出不可用并返回全NaN，未访问行不冒充零解。
% 模块三：保留失败行的输入、系数与逆，支持独立复算；不改变原生产求解入口。
assert(strcmp(factor.kind,'factored-conservative-field-block-ILU0'),'FieldILU:Kind','需要已固定子块因子。');
n=factor.cellCount;d=factor.dimension;assert(isequal(size(rhs),[n d])&&isreal(rhs)&&all(isfinite(rhs),'all'),'FieldILU:RHS','子块右端非法。');
b=rhs';forward=zeros(d,n);solution=forward;value=nan(n,d);
observation=struct('status','in-progress-fixed-factor-replay','valueIsUsable',false,'completedForwardRows',0,'completedBackwardRows',0, ...
    'maximumFiniteIntermediateMagnitude',0,'firstFailure',struct([]),'factorArithmeticChanged',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
for row=1:n
    ids=double(factor.rowPointer(row)):double(factor.diagonalPosition(row))-1;current=b(:,row);before=current;
    if ~isempty(ids),current=current-sum(pagemtimes(factor.values(:,:,ids),reshape(forward(:,factor.column(ids)),d,1,[])),3);end
    forward(:,row)=current;
    if any(~isfinite(current)),fail('forward-lower-sum',current,before,forward,ids);return;end
    observation.completedForwardRows=row;observation.maximumFiniteIntermediateMagnitude=max(observation.maximumFiniteIntermediateMagnitude,max(abs(current)));
end
for row=n:-1:1
    ids=double(factor.diagonalPosition(row))+1:double(factor.rowPointer(row+1))-1;current=forward(:,row);before=current;
    if ~isempty(ids),current=current-sum(pagemtimes(factor.values(:,:,ids),reshape(solution(:,factor.column(ids)),d,1,[])),3);end
    if any(~isfinite(current)),fail('backward-upper-sum',current,before,solution,ids);return;end
    solution(:,row)=factor.inverseDiagonal(:,:,row)*current;
    if any(~isfinite(solution(:,row))),fail('backward-diagonal-multiply',solution(:,row),current,solution,ids);return;end
    observation.completedBackwardRows=observation.completedBackwardRows+1;
    observation.maximumFiniteIntermediateMagnitude=max(observation.maximumFiniteIntermediateMagnitude,max([abs(current);abs(solution(:,row))]));
end
value=solution';observation.status='finite-complete-fixed-factor-replay';observation.valueIsUsable=true;
    function fail(stage,after,initial,partial,ids)
        observation.status='rejected-first-nonfinite-fixed-factor-intermediate';
        observation.firstFailure=struct('stage',stage,'algebraicRow',row,'rhsAtRow',b(:,row),'initialRowValue',initial,'nonfiniteRowValue',after, ...
            'factorEntryIDs',uint32(ids(:)),'dependencyAlgebraicRows',factor.column(ids),'dependencyValues',partial(:,factor.column(ids)), ...
            'factorBlocks',factor.values(:,:,ids),'inverseDiagonal',factor.inverseDiagonal(:,:,row));
    end
end
