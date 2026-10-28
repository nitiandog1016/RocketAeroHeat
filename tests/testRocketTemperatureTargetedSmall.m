function report=testRocketTemperatureTargetedSmall(folder)
%TESTROCKETTEMPERATURETARGETEDSMALL 独立工程可运行的前瞻小域复现，不需要整箭矩阵缓存。
% 模块一：随包三组上下温域制造场只用来核对选面/元数据，不当作实际流场。
% 模块二：复算旧整带缺省与显式目标驱动集合，所有非计时对象须逐值回放。
% 模块三：不继承Q35/Q37方向资格，不运行完整NS或改变正式合格面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
input=rocketTemperatureResearchFixture('temperature_lookahead_unit.mat','outputs/v90_validation/unit01/temperature_lookahead_unit.mat');
a=load(input,'cases');cases=cell(numel(a.cases),1);
for j=1:numel(a.cases)
    c=a.cases{j};
    [legacy,legacyReport]=selectRocketTemperatureLookaheadRows(c.Q,c.info,c.cut,c.gas,c.scale,c.value,1/8);
    assert(isequaln(legacy,c.rows)&&isequaln(legacyReport,c.detail),'TemperatureTargeted:Legacy','旧整带缺省未逐项回放。');
    [rows,detail]=selectRocketTemperatureLookaheadRows(c.Q,c.info,c.cut,c.gas,c.scale,c.value,1/8, ...
        struct('selectionPolicy','original-and-predicted-required-faces'));
    assert(isequal(rows.faceIDs,unique([detail.originalBandFaceIDs;detail.requiredOriginalStateFaceIDs]))&& ...
        detail.originalRowsExactlyRetained&&detail.allPredictedViolatingCellsRepresented&&~detail.constraintRankUsedToSelectFaces);
    cases{j}=struct('source',c,'rows',rows,'detail',detail,'legacyDefaultExactlyReplayed',true);
end
report=struct('status','passed-isolatable-three-manufactured-targeted-selection-cases-with-legacy-replay', ...
    'caseCount',numel(cases),'fixtureSHA256',rocketStreamingFileSHA256(input), ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_targeted_small.mat'),'cases','report','-v7');disp(report);
end
