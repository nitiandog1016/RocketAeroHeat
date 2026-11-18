function replay(out,source)
addpath(out);d=jsondecode(fileread(source)); assert(isfield(d,'researchIdentity') && ~isfield(d,'evidenceIdentity'));
before=ledgerBeforeV263(d); after=buildAerothermalEvidenceLedger(d);
b=before.joint;a=after.joint;assert(strcmp(before.identity.origin,'engineering') && strcmp(after.identity.origin,'research'));
assert(b.pressureHeatFluxEngineeringArea>0 && a.pressureHeatFluxEngineeringArea==0);
assert(b.pressureHeatFluxAvailableArea==a.pressureHeatFluxAvailableArea && a.pressureHeatFluxQualifiedArea==0 && b.pressureHeatFluxQualifiedArea==0);
assert(isequaln(before.faceArea,after.faceArea) && isequaln({before.fields.values},{after.fields.values}));
r=struct('originalFaces',after.faceCount,'beforeM2',b.pressureHeatFluxEngineeringArea,'afterEngineeringM2',a.pressureHeatFluxEngineeringArea,'researchAvailableM2',a.pressureHeatFluxAvailableArea,'jointQualifiedM2',a.pressureHeatFluxQualifiedArea,'numericalFieldsUnchanged',true);
fid=fopen(fullfile(out,'replay.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(r));fclose(fid);
t=testAerothermalEvidenceLedger(fullfile(out,'fixture'));assert(isfile(fullfile(out,'fixture','fixture.json')));
end
