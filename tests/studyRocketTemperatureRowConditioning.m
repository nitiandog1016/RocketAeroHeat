function report=studyRocketTemperatureRowConditioning(folder,stepFolder)
%STUDYROCKETTEMPERATUREROWCONDITIONING 只重放已失败Q40的当前温度投影表示，不重算自由方向或NS。
% 模块一：绑定失败任务、上下文、预览及准备清单，保留全部487原行和187不同条件。
% 模块二：构造可逆差值表示并保存完整支撑基；原1e-10秩/正交与448MiB预算保持。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
contextFile=fullfile(stepFolder,'direction_context.mat');previewFile=fullfile(stepFolder,'linear_checkpoints/repreview00.mat');
c=load(contextFile,'factor','p','preparationFile','preparationSHA256');
assert(strcmp(rocketStreamingFileSHA256(c.preparationFile),c.preparationSHA256));
p=loadRocketPeakSystem(c.preparationFile,c.p.sourceBinding.sourceFile);assert(isequaln(p,c.p));
s=load(previewFile,'rows','summary');[found,selected]=ismember(s.summary.newActiveFaceIDs,s.rows.faceIDs);assert(all(found));
rows=s.rows;
for key={'faceIDs','ownerCV','neighbourCV','normal','QL','QR','H','leftRowsK','rightRowsK', ...
        'baseSamplePrimitive','baseTemperatureK','remainingUpperTemperatureK','sampleRegion','rowTwoNormK','temperatureBoundSigns'}
    rows.(key{1})=rows.(key{1})(selected,:);
end
projector=buildRocketTemperatureProjector(c.factor,rows);
[workingRows,conditioning]=conditionRocketTemperatureRows(rows.H(projector.exactRepresentativeRowIndices,:));
assert(projector.constraintCount==numel(s.summary.newActiveFaceIDs)&& ...
    projector.independentConstraintCount==size(workingRows,1)&&all(projector.normalizedQRDiagonal>1e-10));
report=struct('status','tested-real-Q40-invertible-row-conditioning-projector-only', ...
    'contextFile',contextFile,'contextSHA256',rocketStreamingFileSHA256(contextFile), ...
    'previewFile',previewFile,'previewSHA256',rocketStreamingFileSHA256(previewFile), ...
    'originalConstraintRows',projector.constraintCount,'allDistinctConditions',projector.independentConstraintCount, ...
    'transformedRows',numel(conditioning.transformedRowIndices),'minimumNormalizedQRDiagonal',min(projector.normalizedQRDiagonal), ...
    'newFreeIterations',0,'newProjectedIterations',0,'newActualNSEvaluations',0, ...
    'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_row_conditioning_real.mat'),'projector','workingRows','conditioning','rows','report','-v7');disp(report);
end
