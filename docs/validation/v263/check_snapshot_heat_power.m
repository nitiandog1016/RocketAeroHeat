function check_snapshot_heat_power(folder)
% 核对部分覆盖的原面快照与完整覆盖的风洞壁面快照。不重新求解流场。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root,fullfile(root,'internalFlow','cutCell'));
mainPath=fullfile(root,'outputs','v176_validation','main02','engineering_data.json');
tunnelPath=fullfile(root,'outputs','v263_validation','tmr_resolved_view02','workbench','results.json');
assert(strcmp(rocketStreamingFileSHA256(mainPath),'9518d1d1f44abb83f60664a8c1d072d25797ce16d2182823900cfdf119110b1b'));
assert(strcmp(rocketStreamingFileSHA256(tunnelPath),'2672e38bcce5769277d5b696d25f5b4a9c9781626fa403b4f758772dfe5d0cfd'));

% 模块一：用存档的同一原面顺序复核热流功率、覆盖面积及各区加总。
main=jsondecode(fileread(mainPath));ledger=buildAerothermalEvidenceLedger(main);
q=ledger.fields(strcmp({ledger.fields.key},'heatFlux'));mainPower=q.areaSummary.heatPower;
mask=q.legacyStatus==1;
assert(strcmp(mainPower.status,'partial') && isnan(mainPower.wholeSurfacePowerW));
assert(nnz(mask)==13691 && abs(mainPower.calculatedAreaM2-sum(main.faceArea(mask)))<1e-12);
assert(abs(mainPower.calculatedPowerW-sum(q.values(mask).*main.faceArea(mask))*1000)<1e-7);
regions=ledger.regions;partial=0;covered=0;
for k=1:numel(regions)
    s=regions(k).fields.heatFlux.heatPower;covered=covered+s.calculatedAreaM2;
    if s.calculatedAreaM2>0,partial=partial+s.calculatedPowerW;end
end
assert(abs(partial-mainPower.calculatedPowerW)<1e-7 && abs(covered-mainPower.calculatedAreaM2)<1e-11);
main.evidenceLedger=ledger;
sidecar=compactAerothermalEvidenceLedger(main);
assert(strcmp(sidecar.heatPowerMeasure.csvAreaColumn,'area_m2') && isequal(sidecar.heatPowerMeasure.includedAvailabilityCodes,[1 2]));
writeJSON(fullfile(folder,'main_ledger.json'),sidecar);

% 模块二：以原生环带 CSV 的完整积分为独立参照，并检查缺失曲面积的旧快照。
wind=jsondecode(fileread(tunnelPath));ledger=buildAerothermalEvidenceLedger(wind);
q=ledger.fields(strcmp({ledger.fields.key},'heatFlux'));tunnelPower=q.areaSummary.heatPower;
csv=readtable(fullfile(root,'docs','validation','v263','evidence','tmr_resolved','wall.csv'));
reference=sum(csv.ringArea_m2.*csv.heatIntoWall_W_m2);
assert(strcmp(tunnelPower.status,'complete') && strcmp(tunnelPower.areaBasis,'source-curved'));
assert(abs(tunnelPower.wholeSurfacePowerW-reference)<1e-7 && ...
    abs(tunnelPower.integrationAreaM2-sum(csv.ringArea_m2))<1e-12);
assert(abs(reference-wind.surfaceIntegralAudit.displayFacetedHeatIntegralW)>1);
wind.evidenceLedger=ledger;sidecar=compactAerothermalEvidenceLedger(wind);
assert(isequal(sidecar.heatPowerMeasure.curvedAreaM2ByFaceID(:),wind.wallBinding.curvedAreaPerTriangleM2(:)));
writeJSON(fullfile(folder,'tunnel_ledger.json'),sidecar);
bad=wind;bad.wallBinding=rmfield(bad.wallBinding,'curvedAreaPerTriangleM2');
ledger=buildAerothermalEvidenceLedger(bad);
q=ledger.fields(strcmp({ledger.fields.key},'heatFlux'));assert(strcmp(q.areaSummary.heatPower.status,'area-unavailable'));
bad=wind;bad.wallBinding.curvedAreaPerTriangleM2(1)=bad.wallBinding.curvedAreaPerTriangleM2(1)/2;
ledger=buildAerothermalEvidenceLedger(bad);
q=ledger.fields(strcmp({ledger.fields.key},'heatFlux'));assert(strcmp(q.areaSummary.heatPower.status,'area-unavailable'));
report=struct('mainSnapshot',mainPower,'windTunnelSnapshot',tunnelPower,'windTunnelRingIntegralW',reference, ...
    'mainComputedFaces',nnz(mask),'mainOriginalFaces',numel(mask),'newJointAreaM2',0,'newFlowSolve',false);
writeJSON(fullfile(folder,'heat_power.json'),report);
end

function writeJSON(filename,data)
fid=fopen(filename,'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(data,'ConvertInfAndNaN',true));
end
