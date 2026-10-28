function [value,report,evidence]=finishRocketTemperatureActiveDirection(packed,factor,rhs,metric,rows,options,free,freeDetail,freeCoordinates)
%FINISHROCKETTEMPERATUREACTIVEDIRECTION 对已求出的当前自由方向执行原加约束循环。
% 模块一：与正常入口共用同一后半程，不重算自由方向，不自动授予检查点复用资格。
% 模块二：恢复调用者须先核对当前状态、几何、全部产物SHA、坐标提升和原C/W目标。
% 模块三：当前状态只增加朝域外条件，原秩/可行性检查与最多三次受限求解保持。
assert(isequal(size(free),size(rhs),size(freeCoordinates))&&all(isfinite(free),'all')&& ...
    all(isfinite(freeCoordinates),'all')&&freeDetail.nativeIterations<=options.maximumIterations&& ...
    options.maximumActiveSolves>=1&&options.maximumActiveSolves<=3, ...
    'TemperatureActive:FreePrefix','自由前缀或当前有界求解选项不符。');
value=free;detail=freeDetail;coordinates=freeCoordinates;[tangent,bound,violating]=inspect(value);active=zeros(0,1);cycles={};
if ~isempty(options.checkpointFolder)
    file=fullfile(options.checkpointFolder,'free_direction.mat');assert(~isfile(file));
    save(file,'free','freeDetail','freeCoordinates','tangent','bound','violating','-v7');
end
stopReason='free-direction-points-inside-temperature-band';
for cycle=1:options.maximumActiveSolves
    if isempty(violating),break;end
    added=setdiff(violating,active);if isempty(added),stopReason='active-row-arithmetic-feasibility-not-met';break;end
    active=unique([active;added]);selected=rows;selected.H=rows.H(active,:);selected.faceIDs=rows.faceIDs(active);
    projector=buildRocketTemperatureProjector(factor,selected);
    [value,detail,coordinates]=solveRocketTemperatureProjectedDirection(packed,factor,rhs,metric,projector,selected,options.maximumIterations);
    [tangent,bound,violating]=inspect(value);
    item=struct('activeRowIndices',active,'activeFaceIDs',rows.faceIDs(active),'newRowIndices',added, ...
        'projector',projector,'value',value,'detail',detail,'coordinates',coordinates, ...
        'temperatureDerivativeK',tangent,'localDotRoundoffBoundK',bound,'remainingViolatingRowIndices',violating);
    cycles{end+1}=item; %#ok<AGROW>
    if ~isempty(options.checkpointFolder)
        file=fullfile(options.checkpointFolder,sprintf('active_cycle%02d.mat',cycle));assert(~isfile(file));save(file,'item','-v7');
    end
    stopReason='current-band-direction-feasible';
    if ~isempty(violating),stopReason='bounded-active-solve-budget-exhausted';end
end
evidence=struct('free',free,'freeDetail',freeDetail,'freeCoordinates',freeCoordinates,'cycles',{cycles}, ...
    'selectedCoordinates',coordinates,'selectedDetail',detail,'rows',rows);
report=struct('method','current-temperature-band-add-only-feasible-direction-QR-LSQR-prototype', ...
    'options',options,'freeLinearReport',freeDetail,'linearReport',detail,'activeSolveCount',numel(cycles), ...
    'totalNativeIterations',freeDetail.nativeIterations+sum(cellfun(@(x)x.detail.nativeIterations,cycles)), ...
    'activeRowIndices',active,'activeFaceIDs',rows.faceIDs(active),'bandFaceIDs',rows.faceIDs, ...
    'temperatureBoundSigns',rows.temperatureBoundSigns,'temperatureDerivativeK',tangent, ...
    'localDotRoundoffBoundK',bound,'remainingViolatingRowIndices',violating,'bandDirectionFeasible',isempty(violating), ...
    'stopReason',stopReason,'oldActiveSetInherited',false,'coolingRowsMayRemainFree',true, ...
    'constrainedOptimumClaimed',false,'nonlinearTemperatureBoundGuaranteed',false,'physicalBoundaryConditionAdded',false, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'physicalSolution',false);
    function [t,b,v]=inspect(direction)
        t=rows.H*direction(:);b=128*eps*full(abs(rows.H)*abs(direction(:)));
        % 只容纳十项点积的参考舍入，原物理温度和NS门槛不变。
        v=find(rows.temperatureBoundSigns.*t>b);
    end
end
