function report=testRocketTemperatureTargetedSelection(folder)
%TESTROCKETTEMPERATURETARGETEDSELECTION 验证目标驱动选面及每项面元数据，不运行新NS。
% 模块一：上下温域制造场与真实Q35夹具的旧整带缺省输出逐项回放。
% 模块二：新集合严格等于原10K加全部预测所需面，H/符号/法向/两侧状态均对齐。
% 模块三：使用已保存Q37自由前缀，独立确认25单元均保留且不按秩选择条件。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
a=load(fullfile(root,'outputs/v90_validation/unit01/temperature_lookahead_unit.mat'),'cases');cases=cell(numel(a.cases),1);
for j=1:numel(a.cases)
    c=a.cases{j};
    [legacy,legacyReport]=selectRocketTemperatureLookaheadRows(c.Q,c.info,c.cut,c.gas,c.scale,c.value,1/8);
    assert(isequaln(legacy,c.rows)&&isequaln(legacyReport,c.detail),'TemperatureTargeted:Legacy','旧缺省选面发生变化。');
    [rows,detail]=selectRocketTemperatureLookaheadRows(c.Q,c.info,c.cut,c.gas,c.scale,c.value,1/8, ...
        struct('selectionPolicy','original-and-predicted-required-faces'));
    check(rows,detail,legacy);
    cases{j}=struct('source',c,'rows',rows,'detail',detail,'legacyDefaultExactlyReplayed',true);
end
step=fullfile(root,'outputs/v90_validation/continuous01/step01');
c=load(fullfile(step,'direction_context.mat'),'p','initialRows');
f=load(fullfile(step,'linear_checkpoints/lookahead_free_prefix.mat'),'free','originalRows');
assert(isequaln(c.initialRows,f.originalRows));
s=load(c.p.sourceBinding.sourceFile,'result','gas');d=load(c.p.geometryBinding.file,'domain');
[fullBand,fullReport]=selectRocketTemperatureLookaheadRows(s.result.Q,s.result.info,d.domain,s.gas,s.result.residualScale,f.free,1/8);
[realRows,realReport]=selectRocketTemperatureLookaheadRows(s.result.Q,s.result.info,d.domain,s.gas,s.result.residualScale,f.free,1/8, ...
    struct('selectionPolicy','original-and-predicted-required-faces'));
check(realRows,realReport,fullBand);
assert(numel(realReport.predictedViolatingCellIDs)==25&&numel(c.initialRows.faceIDs)==124&& ...
    numel(fullBand.faceIDs)==193&&numel(realRows.faceIDs)==131);
for key={'ownerCV','neighbourCV'}
    assert(all(realRows.(key{1})>=1&realRows.(key{1})<=size(s.result.Q,1)));
end
report=struct('status','passed-legacy-selector-replay-targeted-all-required-CVs-and-complete-row-metadata', ...
    'currentStateSHA256',rocketStreamingFileSHA256(c.p.sourceBinding.sourceFile), ...
    'currentFreePrefixSHA256',rocketStreamingFileSHA256(fullfile(step,'linear_checkpoints/lookahead_free_prefix.mat')), ...
    'currentOriginalRows',124,'currentPredictedViolatingCells',25,'currentTargetedRows',131,'candidateBandRows',193, ...
    'selectionUsesConstraintRank',false,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_targeted_unit.mat'),'cases','realRows','realReport','fullReport','report','-v7');disp(report);
end
function check(rows,detail,allRows)
expected=unique([detail.originalBandFaceIDs;detail.requiredOriginalStateFaceIDs]);
assert(isequal(rows.faceIDs,expected)&&detail.originalRowsExactlyRetained&& ...
    detail.allPredictedViolatingCellsRepresented&&~detail.constraintRankUsedToSelectFaces&& ...
    ~detail.allExtendedBandFacesSelected);
[found,slots]=ismember(rows.faceIDs,allRows.faceIDs);assert(all(found));
fields={'faceIDs','ownerCV','neighbourCV','normal','QL','QR','H','leftRowsK','rightRowsK', ...
    'baseSamplePrimitive','baseTemperatureK','remainingUpperTemperatureK','sampleRegion','rowTwoNormK','temperatureBoundSigns'};
for key=fields,assert(isequaln(rows.(key{1}),allRows.(key{1})(slots,:)),'TemperatureTargeted:Metadata','选面后元数据错位。');end
assert(isequal(detail.candidateExtendedBandFaceIDs,allRows.faceIDs));
end
