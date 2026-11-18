function value=compactAerothermalEvidenceLedger(data)
% 将工作台各变量和区域的台账汇总写入 CSV 的 JSON 旁车。
ledger=data.evidenceLedger;
value=struct('schemaVersion',ledger.schemaVersion,'kind',ledger.kind,'identity',ledger.identity, ...
    'fields',struct([]),'regions',struct([]),'joint',withoutFaceIDs(ledger.joint), ...
    'acceptedEvidence',ledger.acceptedEvidence,'rejectedEvidence',ledger.rejectedEvidence);
for i=1:numel(ledger.fields)
    f=ledger.fields(i);entry=struct('key',f.key,'label',f.label,'unit',f.unit,'valueSHA256',f.valueSHA256,'areaSummary',f.areaSummary);
    if isempty(value.fields),value.fields=entry;else,value.fields(end+1)=entry;end
end
for i=1:numel(ledger.regions)
    r=ledger.regions(i);entry=struct('name',r.name,'area',r.area,'faceCount',numel(r.faceIDs),'fields',r.fields,'joint',withoutFaceIDs(r.joint));
    if isempty(value.regions),value.regions=entry;else,value.regions(end+1)=entry;end
end
% 导出方法说明与当前热流采用的积分面积。旋转曲面的逐面真面积只在此处附一次。
q=ledger.fields(strcmp({ledger.fields.key},'heatFlux'));
if ~isempty(q) && isfield(q.areaSummary,'heatPower')
    basis=q.areaSummary.heatPower.areaBasis;
    measure=struct('areaBasis',basis,'includedAvailabilityCodes',[1 2]);
    if strcmp(basis,'source-curved')
        measure.curvedAreaM2ByFaceID=data.wallBinding.curvedAreaPerTriangleM2(:);
    elseif strcmp(basis,'triangles')
        measure.csvAreaColumn='area_m2';
    end
    value.heatPowerMeasure=measure;
end
end
function value=withoutFaceIDs(value)
names=fieldnames(value);remove=names(endsWith(names,'FaceIDs'));if ~isempty(remove),value=rmfield(value,remove);end
end
