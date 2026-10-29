function report=testRocketTemperatureRepreview(folder)
%TESTROCKETTEMPERATUREREPREVIEW 以独立制造的耦合线性系统验证投影后新暴露单元的通用补全。
% 模块一：32CV两种代数排列；线性C/W是数值制造系统，不冒称此Q的物理NS导数。
% 模块二：首次候选和首次投影后分别出现不同越温CV；最终方向与稠密零空间解比较。
% 模块三：一轮预算不足、齐次约束仍无法保证有限步温域、零方向及原行冲突均明确记录。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
p=load(rocketTemperatureResearchFixture('parent_least_squares_unit_evidence.mat'),'cases');
m=load(rocketTemperatureResearchFixture('peak_weighted_unit_evidence.mat'),'cases');
g=load(rocketTemperatureResearchFixture('temperature_lookahead_unit.mat'),'cases');
fixture=load(fullfile(root,'tests/data/temperature_repreview/late_exposure_closes_moderate.mat'));
stressFixture=load(fullfile(root,'tests/data/temperature_repreview/late_exposure_closes.mat'));
stallFixture=load(fullfile(root,'tests/data/temperature_repreview/curvature_stall.mat'));
cut=g.cases{3}.cut;gas=g.cases{3}.gas;scale=fixture.scale;Q=fixture.Q;n=size(Q,1);
ids=find(cut.openFaces.area>0&cut.openFaces.cellMinus>0&cut.openFaces.cellPlus>0);
[~,audit]=cutCellExactEulerFlux(Q(cut.openFaces.cellMinus(ids),:),Q(cut.openFaces.cellPlus(ids),:),cut.openFaces.normal(ids,:),gas);
audit.faceID=ids;info=struct('exactEulerFluxAudit',audit);base=selectRocketTemperatureBandRows(Q,info,cut,gas,scale);
preview=struct('Q',Q,'info',info,'cut',cut,'gas',gas,'scale',scale,'targetAlpha',fixture.targetAlpha, ...
    'selectionPolicy','original-and-iteratively-predicted-cell-and-interface-faces');cases=cell(2,1);
highQ=stressFixture.Q;
[~,highAudit]=cutCellExactEulerFlux(highQ(cut.openFaces.cellMinus(ids),:),highQ(cut.openFaces.cellPlus(ids),:),cut.openFaces.normal(ids,:),gas);
highAudit.faceID=ids;highInfo=struct('exactEulerFluxAudit',highAudit);
highBase=selectRocketTemperatureBandRows(highQ,highInfo,cut,gas,scale);
highPreview=preview;highPreview.Q=highQ;highPreview.info=highInfo;
for j=1:2
    a=p.cases{j};b=m.cases{j};factor=buildRocketAugmentedPatchQR(a.packed,a.factor,b.metric,b.detail.penaltyInOriginalResidualUnits);
    A=[b.base.C;factor.penaltyInOriginalResidualUnits*b.W];
    rhs=reshape((b.base.C')\((A'*A)*fixture.desiredFree(:)),n,5);
    [value,detail,evidence]=solveRocketTemperatureActiveDirection(a.packed,factor,rhs,b.metric,base, ...
        struct('maximumIterations',900,'lookaheadPreview',preview));
    save(fullfile(folder,sprintf('attempt_order%d.mat',j)),'value','detail','evidence','rhs','Q','scale','gas','-v7');
    assert(detail.activeSolveCount>=2&&detail.repreview.previewCount>=3&&detail.bandDirectionFeasible&& ...
        detail.repreview.targetCandidateThermallyFeasible&&detail.repreview.forcedRequiredRowsActuallySolved);
    assert(norm(evidence.free-fixture.desiredFree,'fro')/max(1,norm(fixture.desiredFree,'fro'))<1e-8);
    assert(~isempty(setdiff(detail.repreview.history{2}.predictedViolatingCellIDs,detail.repreview.history{1}.predictedViolatingCellIDs)), ...
        'TemperatureRepreview:Manufacture','制造场必须实际出现投影后新增越域单元。');
    denseChecks=zeros(numel(evidence.cycles),1);
    for k=1:numel(evidence.cycles)
        c=evidence.cycles{k};H=full(evidence.rows.H(c.activeRowIndices,:));Z=null(H);
        dense=Z*((A*Z)\[rhs(:);zeros(5*n,1)]);
        denseChecks(k)=norm(c.value(:)-dense)/max(1,norm(dense));assert(denseChecks(k)<1e-8);
        if k>1,assert(all(ismember(evidence.cycles{k-1}.activeFaceIDs,c.activeFaceIDs)));end
    end
    [found,old]=ismember(base.faceIDs,evidence.rows.faceIDs);assert(all(found)&&isequal(base.H,evidence.rows.H(old,:))&& ...
        isequal(base.temperatureBoundSigns,evidence.rows.temperatureBoundSigns(old)));
    [budgetValue,budget,budgetEvidence]=solveRocketTemperatureActiveDirection(a.packed,factor,rhs,b.metric,base, ...
        struct('maximumIterations',900,'maximumActiveSolves',1,'lookaheadPreview',preview));
    assert(budget.activeSolveCount==1&&budget.repreview.previewCount==2&&~budget.repreview.targetCandidateThermallyFeasible&& ...
        ~budget.repreview.forcedRequiredRowsActuallySolved&&~budget.bandDirectionFeasible);
    stallRhs=reshape((b.base.C')\((A'*A)*stallFixture.desiredFree(:)),n,5);
    [stallValue,stall,stallEvidence]=solveRocketTemperatureActiveDirection(a.packed,factor,stallRhs,b.metric,highBase, ...
        struct('maximumIterations',900,'lookaheadPreview',highPreview));
    assert(~stall.repreview.targetCandidateThermallyFeasible&& ...
        strcmp(stall.stopReason,'repreview-has-no-new-actionable-original-condition'));
    % 高动能制造场有限步温度可合法，但更严格的原128eps点积门槛可能仍拒绝；如实保留。
    stressRhs=reshape((b.base.C')\((A'*A)*stressFixture.desiredFree(:)),n,5);
    [stressValue,stress,stressEvidence]=solveRocketTemperatureActiveDirection(a.packed,factor,stressRhs,b.metric,highBase, ...
        struct('maximumIterations',900,'lookaheadPreview',highPreview));
    stressDerivative=stressEvidence.rows.H*stressValue(:);
    stressBound=128*eps*full(abs(stressEvidence.rows.H)*abs(stressValue(:)));
    expectedGate=all(stressEvidence.rows.temperatureBoundSigns.*stressDerivative<=stressBound)&& ...
        stress.repreview.forcedRequiredRowsActuallySolved;
    assert(stress.repreview.targetCandidateThermallyFeasible&&stress.bandDirectionFeasible==expectedGate);
    [zero,zeroReport]=solveRocketTemperatureActiveDirection(a.packed,factor,zeros(size(rhs)),b.metric,base, ...
        struct('maximumIterations',900,'lookaheadPreview',preview));
    assert(all(zero==0,'all')&&zeroReport.activeSolveCount==0&&zeroReport.repreview.targetCandidateThermallyFeasible);
    changed=evidence.rows;changed.temperatureBoundSigns(1)=-changed.temperatureBoundSigns(1);
    conflict=false;try,mergeRocketTemperatureRequiredRows(evidence.rows,changed);catch e,conflict=strcmp(e.identifier,'TemperatureRepreview:Retention');end
    assert(conflict,'TemperatureRepreview:Retention','共同面的上下侧冲突必须拒绝。');
    % 空集合也要来自同一气体与缩放；两个参数次序都拒绝错误身份。
    empty=base;
    for key={'faceIDs','ownerCV','neighbourCV','normal','QL','QR','H','leftRowsK','rightRowsK', ...
            'baseSamplePrimitive','baseTemperatureK','remainingUpperTemperatureK','sampleRegion','rowTwoNormK','temperatureBoundSigns'}
        empty.(key{1})=empty.(key{1})([],:);
    end
    emptyIdentityRejects=false(2,2);
    for changedField=1:2
        wrong=empty;
        if changedField==1,wrong.gas.R=wrong.gas.R*1.01;else,wrong.scale(1)=wrong.scale(1)*1.01;end
        for order=1:2
            try
                if order==1,mergeRocketTemperatureRequiredRows(wrong,base);else,mergeRocketTemperatureRequiredRows(base,wrong);end
            catch e,emptyIdentityRejects(changedField,order)=strcmp(e.identifier,'TemperatureRepreview:Identity');end
        end
    end
    assert(all(emptyIdentityRejects,'all'),'TemperatureRepreview:Identity','空集合不能跳过气体或缩放身份核对。');
    cases{j}=struct('packed',a.packed,'factor',factor,'metric',b.metric,'C',b.base.C,'W',b.W,'rhs',rhs, ...
        'Q',Q,'cut',cut,'gas',gas,'scale',scale,'info',info,'baseRows',base,'value',value,'detail',detail,'evidence',evidence, ...
        'denseChecks',denseChecks,'budgetValue',budgetValue,'budget',budget,'budgetEvidence',budgetEvidence, ...
        'stallValue',stallValue,'stall',stall,'stallEvidence',stallEvidence,'stallRhs',stallRhs, ...
        'highQ',highQ,'highInfo',highInfo,'highBaseRows',highBase, ...
        'stressValue',stressValue,'stress',stress,'stressEvidence',stressEvidence,'stressRhs',stressRhs, ...
        'zeroReport',zeroReport,'sideConflictRejected',conflict,'emptyIdentityRejects',emptyIdentityRejects);
end
report=struct('status','passed-two-orderings-late-exposure-dense-directions-budget-stall-zero-and-retention', ...
    'manufacturedLinearSystemNotPhysicalNSJacobian',true,'originalTemperatureRangeK',gas.temperatureRange, ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_repreview_unit.mat'),'cases','report','fixture','stallFixture','stressFixture','-v7');disp(report);
end
