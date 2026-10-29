function [value,report,evidence]=finishRocketTemperatureRepreviewDirection(packed,factor,rhs,metric,rows,options,free,freeDetail,freeCoordinates,preview,firstSelection,originalRows)
%FINISHROCKETTEMPERATUREREPREVIEWDIRECTION 在同一当前状态内有预算地复查投影后的实际守恒候选。
% 模块一：复用本次已计算自由方向，最多三次受限求解、四次候选温域预览，不重算自由段。
% 模块二：所有旧行/活跃面只增不减；每次补全全部新越温CV原态面及可检查的越温界面。
% 模块三：有限步已越域的必要面显式激活，即使一阶导数朝内也不能冒充有限步合法。
% 模块四：仅检查指定alpha的单元/内部界面温域；开放边界及原完整NS、下降/稳态仍另验。
assert(options.maximumActiveSolves>=1&&options.maximumActiveSolves<=3&& ...
    freeDetail.nativeIterations<=options.maximumIterations&&isequal(size(free),size(rhs),size(freeCoordinates)), ...
    'TemperatureRepreview:Options','需要当前自由前缀与最多三次受限预算。');
value=free;detail=freeDetail;coordinates=freeCoordinates;cycles={};previews={};history={};
activeFaces=zeros(0,1);allRequired=zeros(0,1);allCells=zeros(0,1);allPure=zeros(0,1);allInterfaces=zeros(0,1);
stopReason='bounded-repreview-budget-exhausted';targetFeasible=false;allEvaluable=true;
if ~isempty(options.checkpointFolder)
    file=fullfile(options.checkpointFolder,'free_direction.mat');assert(~isfile(file));
    save(file,'free','freeDetail','freeCoordinates','-v7');
end
for round=0:options.maximumActiveSolves
    availableBefore=rows.faceIDs;
    if round==0
        selection=firstSelection;selected=rows;
    else
        [selected,selection]=selectRocketTemperatureCellInterfaceRows(preview.Q,preview.info,preview.cut, ...
            preview.gas,preview.scale,value,preview.targetAlpha);
    end
    rows=mergeRocketTemperatureRequiredRows(rows,selected);
    evaluable=selection.previewCellStateEvaluable;allEvaluable=allEvaluable&&evaluable;
    required=selection.requiredOriginalStateFaceIDs;interfaces=zeros(0,1);assessed=0;blocked=0;fullyAssessed=false;
    if evaluable
        ip=selection.interfacePreview;interfaces=ip.predictedViolatingInterfaceFaceIDs;
        assessed=numel(ip.assessedFaceIDs);blocked=numel(ip.blockedFaceIDs);fullyAssessed=ip.allInteriorFacesAssessed;
    end
    allCells=union(allCells,selection.predictedViolatingCellIDs);allPure=union(allPure,required);
    allInterfaces=union(allInterfaces,interfaces);allRequired=union(allRequired,union(required,interfaces));
    [tangent,bound,violating]=inspect(value);necessary=union(allRequired,rows.faceIDs(violating));
    needed=union(activeFaces,necessary);newActive=setdiff(needed,activeFaces);
    targetFeasible=evaluable&&isempty(selection.predictedViolatingCellIDs)&&fullyAssessed&&isempty(interfaces);
    summary=struct('round',round,'afterProjectedSolves',numel(cycles),'previewCellStateEvaluable',evaluable, ...
        'predictedViolatingCellIDs',selection.predictedViolatingCellIDs,'requiredOriginalStateFaceIDs',required, ...
        'predictedViolatingInterfaceFaceIDs',interfaces,'assessedInterfaceCount',assessed,'blockedInterfaceCount',blocked, ...
        'allInteriorFacesAssessed',fullyAssessed,'addedFaceIDs',setdiff(rows.faceIDs,availableBefore), ...
        'selectedFaceIDs',rows.faceIDs,'forcedActiveFaceIDs',allRequired,'outwardFaceIDs',rows.faceIDs(violating), ...
        'activeFacesBeforeSolve',activeFaces,'newActiveFaceIDs',newActive,'targetCandidateThermallyFeasible',targetFeasible);
    if round==0,summary.addedFaceIDs=setdiff(rows.faceIDs,originalRows.faceIDs);end
    history{end+1}=summary;previews{end+1}=selection; %#ok<AGROW>
    if ~isempty(options.checkpointFolder)
        file=fullfile(options.checkpointFolder,sprintf('repreview%02d.mat',round));assert(~isfile(file));
        save(file,'value','selection','summary','rows','-v7');
    end
    if targetFeasible&&isempty(violating)&&all(ismember(allRequired,activeFaces))
        stopReason='specified-target-candidate-temperature-checked-with-all-current-rows-feasible';break;
    end
    if isempty(newActive)
        stopReason='repreview-has-no-new-actionable-original-condition';break;
    end
    if round==options.maximumActiveSolves,break;end
    activeFaces=needed;[found,active]=ismember(activeFaces,rows.faceIDs);assert(all(found));
    selected=rows;selected.H=rows.H(active,:);selected.faceIDs=activeFaces;
    projector=buildRocketTemperatureProjector(factor,selected);
    [value,detail,coordinates]=solveRocketTemperatureProjectedDirection(packed,factor,rhs,metric, ...
        projector,selected,options.maximumIterations);
    [tangent,bound,violating]=inspect(value);
    item=struct('availableFaceIDs',rows.faceIDs,'activeFaceIDs',activeFaces,'newActiveFaceIDs',newActive, ...
        'forcedActiveFaceIDs',allRequired,'previewHistoryIndex',numel(history), ...
        'projector',projector,'value',value,'detail',detail,'coordinates',coordinates, ...
        'temperatureDerivativeK',tangent,'localDotRoundoffBoundK',bound,'remainingViolatingFaceIDs',rows.faceIDs(violating));
    cycles{end+1}=item; %#ok<AGROW>
    if ~isempty(options.checkpointFolder)
        file=fullfile(options.checkpointFolder,sprintf('active_cycle%02d.mat',numel(cycles)));assert(~isfile(file));
        save(file,'item','-v7');
    end
end
% 最终集合的索引只在全部补全结束后建立，避免插入小面号造成旧循环索引错位。
for j=1:numel(cycles)
    [found,cycles{j}.activeRowIndices]=ismember(cycles{j}.activeFaceIDs,rows.faceIDs);assert(all(found));
    [found,cycles{j}.newRowIndices]=ismember(cycles{j}.newActiveFaceIDs,rows.faceIDs);assert(all(found));
    [found,cycles{j}.availableRowIndices]=ismember(cycles{j}.availableFaceIDs,rows.faceIDs);assert(all(found));
    [found,cycles{j}.remainingViolatingRowIndices]=ismember(cycles{j}.remainingViolatingFaceIDs,rows.faceIDs);assert(all(found));
end
[~,active]=ismember(activeFaces,rows.faceIDs);[tangent,bound,violating]=inspect(value);
forcedSolved=all(ismember(allRequired,activeFaces));
selection=firstSelection;
selection.selectionPolicy='original-and-iteratively-predicted-cell-and-interface-faces';
selection.initialSelectedBandFaceIDs=firstSelection.selectedBandFaceIDs;
selection.predictedViolatingCellIDs=allCells;selection.requiredOriginalStateFaceIDs=allPure;
selection.predictedViolatingInterfaceFaceIDs=allInterfaces;selection.selectedBandFaceIDs=rows.faceIDs;
selection.selectedTemperatureBandK=rows.options.temperatureBandK;
audit=preview.info.exactEulerFluxAudit;T=audit.samplePrimitive(:,5);
margin=min(T-preview.gas.temperatureRange(1),preview.gas.temperatureRange(2)-T);
selection.candidateExtendedBandFaceIDs=audit.faceID(margin<=rows.options.temperatureBandK);
selection.allPredictedViolatingCellsRepresented=allEvaluable;selection.allPredictedViolatingInterfacesRepresented=allEvaluable;
selection.originalRowsExactlyRetained=true;
rows.selectionPolicy='all-original-and-all-iteratively-predicted-required-current-temperature-rows';
evidence=struct('free',free,'freeDetail',freeDetail,'freeCoordinates',freeCoordinates,'cycles',{cycles}, ...
    'selectedCoordinates',coordinates,'selectedDetail',detail,'rows',rows,'previewSelections',{previews});
report=struct('method','bounded-current-candidate-repreview-add-only-QR-LSQR', ...
    'options',options,'freeLinearReport',freeDetail,'linearReport',detail,'activeSolveCount',numel(cycles), ...
    'totalNativeIterations',freeDetail.nativeIterations+sum(cellfun(@(x)x.detail.nativeIterations,cycles)), ...
    'activeRowIndices',active,'activeFaceIDs',activeFaces,'bandFaceIDs',rows.faceIDs, ...
    'temperatureBoundSigns',rows.temperatureBoundSigns,'temperatureDerivativeK',tangent, ...
    'localDotRoundoffBoundK',bound,'remainingViolatingRowIndices',violating, ...
    'bandDirectionFeasible',isempty(violating)&&forcedSolved,'stopReason',stopReason, ...
    'oldActiveSetInherited',false,'coolingRowsMayRemainFree',true,'constrainedOptimumClaimed',false, ...
    'nonlinearTemperatureBoundGuaranteed',false,'physicalBoundaryConditionAdded',false, ...
    'lookaheadSelection',selection,'repreview',struct('targetAlpha',preview.targetAlpha,'maximumPreviews',options.maximumActiveSolves+1, ...
        'history',{history},'previewCount',numel(history),'targetCandidateThermallyFeasible',targetFeasible, ...
        'checkedScope','all-CV-and-interior-Riemann-temperature-only-not-open-boundary-or-full-NS', ...
        'forcedRequiredRowsActuallySolved',forcedSolved,'allPreviousRowsRetained',true, ...
        'allPreviousActiveFacesRetained',true,'initialFreeDirectionRecomputed',false), ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'physicalSolution',false);
    function [t,b,v]=inspect(direction)
        t=rows.H*direction(:);b=128*eps*full(abs(rows.H)*abs(direction(:)));
        v=find(rows.temperatureBoundSigns.*t>b);
    end
end
