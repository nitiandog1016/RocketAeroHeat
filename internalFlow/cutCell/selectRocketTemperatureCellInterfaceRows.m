function [rows,report]=selectRocketTemperatureCellInterfaceRows(Q,info,cut,gas,scale,directionU,targetAlpha)
%SELECTROCKETTEMPERATURECELLINTERFACEROWS 同时前瞻候选单元与真实内部Riemann采样温度。
% 模块一：完整复用原10K与全部预测越温CV的原态面，原H/符号不改。
% 模块二：原内核检查全部合法输入单元对，新增真正预测越温界面的当前解析H。
% 模块三：星区/扇区保持真实Riemann导数，不当作单元EOS；选面不读取秩或NS结果。
% 模块四：只做一次当前自由候选预览，不能保证投影后候选合法，原完整NS仍独立验收。
[rows,report]=selectRocketTemperatureLookaheadRows(Q,info,cut,gas,scale,directionU,targetAlpha, ...
    struct('selectionPolicy','original-and-predicted-required-faces'));
report.method='current-conservative-cell-and-original-Riemann-interface-temperature-preview-union';
report.selectionPolicy='original-and-predicted-cell-and-interface-faces';
report.interfacePreview=[];report.constraintRankUsedToSelectFaces=false;
if ~report.previewCellStateEvaluable
    report.stopReason='global-candidate-state-invalid-no-interface-preview';return;
end
candidate=rocketShellTrialState(Q,directionU,scale,targetAlpha,gas,'conservative');
preview=findRocketCandidateTemperatureFaces(candidate,cut,gas);report.interfacePreview=preview;
required=preview.predictedViolatingInterfaceFaceIDs;candidateT=preview.predictedInterfaceTemperatureK;
candidateSigns=ones(numel(required),1);candidateSigns(candidateT<gas.temperatureRange(1))=-1;
[included,oldSlots]=ismember(required,rows.faceIDs);
assert(all(rows.temperatureBoundSigns(oldSlots(included))==candidateSigns(included)), ...
    'TemperatureLookahead:OppositeSides','同一已有行同时遇到相反温域侧，受限原型拒绝而不改原符号。');
oldRows=rows;extraIDs=required(~included);
if ~isempty(extraIDs)
    extra=buildRocketInterfaceTemperatureRows(Q,cut,gas,scale,extraIDs);
    extra.temperatureBoundSigns=candidateSigns(~included);
    if isempty(rows.faceIDs)
        rows=extra;rows.options=oldRows.options;
        rows.currentCompleteInteriorFaceCount=oldRows.currentCompleteInteriorFaceCount;
        rows.oldActiveSetInherited=false;
    else
        keys={'faceIDs','ownerCV','neighbourCV','normal','QL','QR','H','leftRowsK','rightRowsK', ...
            'baseSamplePrimitive','baseTemperatureK','remainingUpperTemperatureK','sampleRegion','rowTwoNormK','temperatureBoundSigns'};
        [~,order]=sort([rows.faceIDs;extra.faceIDs]);
        for key=keys
            joined=[rows.(key{1});extra.(key{1})];rows.(key{1})=joined(order,:);
        end
    end
end
[retained,slots]=ismember(oldRows.faceIDs,rows.faceIDs);
assert(all(retained)&&isequal(oldRows.H,rows.H(slots,:))&& ...
    isequal(oldRows.temperatureBoundSigns,rows.temperatureBoundSigns(slots)), ...
    'TemperatureLookahead:Retention','原10K/CV所需行和符号未完整保留。');
assert(all(ismember(required,rows.faceIDs)),'TemperatureLookahead:Retention','有预测越温界面未保留。');
rows.selectionPolicy='all-original-10K-plus-all-predicted-CV-and-assessable-Riemann-interface-required-faces';
radius=max([10;min(rows.baseTemperatureK-gas.temperatureRange(1),gas.temperatureRange(2)-rows.baseTemperatureK)]);
rows.options.temperatureBandK=radius;
audit=info.exactEulerFluxAudit;T=audit.samplePrimitive(:,5);margin=min(T-gas.temperatureRange(1),gas.temperatureRange(2)-T);
report.selectedTemperatureBandK=radius;report.candidateExtendedBandFaceIDs=audit.faceID(margin<=radius);
report.selectedBandFaceIDs=rows.faceIDs;report.predictedViolatingInterfaceFaceIDs=required;
report.addedInterfaceFaceIDs=extraIDs;report.predictedInterfaceTemperatureK=candidateT;
report.allPredictedViolatingInterfacesRepresented=true;report.allPredictedViolatingCellsRepresented=true;
report.allExtendedBandFacesSelected=false;report.originalRowsExactlyRetained=true;
report.stopReason='retained-all-original-cell-required-and-all-assessable-predicted-interface-rows';
end
