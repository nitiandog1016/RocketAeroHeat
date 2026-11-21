function mesh=attachThermistorOuterDomain(mesh)
% 模块一：把 30D 外域对照绑定到当前原面、20D 体网、冻结工况及真实参照帧。
root=fileparts(mfilename('fullpath'));key=jsondecode(fileread(fullfile(root,'cases','thermistor_b','outer30.json')));
[folder,guard]=openThermistorArchive('outer30'); %#ok<NASGU>
base=mesh.thermistorLead;if isfield(base,'referenceNative'),base=base.referenceNative;end
id=mesh.evidenceIdentity.sourceFrame;c=mesh.hybridCondition;
assert(strcmp(key.baseMeshSHA256,mesh.evidenceIdentity.meshSHA256) && strcmp(key.baseCaseSHA256,base.caseSHA256) && ...
    strcmp(key.bodySHA256,id.bodySHA256) && strcmp(base.sourceFrameSHA256,id.sha256) && strcmp(base.wallSHA256,id.wallSHA256));
bodyFile=fullfile(folder,'body.json');caseFile=fullfile(folder,'case.json');tableFile=fullfile(folder,'wall.csv');
s=jsondecode(fileread(fullfile(folder,'summary.json')));body=jsondecode(fileread(bodyFile));spec=jsondecode(fileread(caseFile));
assert(strcmp(cm011FileSHA256(bodyFile),s.bodySHA256) && strcmp(s.bodySHA256,key.bodySHA256) && ...
    strcmp(cm011FileSHA256(caseFile),s.caseSHA256) && strcmp(s.caseSHA256,key.caseSHA256) && ...
    strcmp(cm011FileSHA256(tableFile),s.wallSHA256) && strcmp(s.meshSHA256,key.meshSHA256) && strcmp(s.sourceFrameSHA256,key.sourceFrameSHA256));
assert(isequal(body.V,mesh.V) && isequal(body.F,mesh.F) && isequal(c.velocityMps,[-spec.speedMps 0 0]) && ...
    c.ambientPressurePa==spec.airPa && c.ambientTemperatureK==spec.airK && spec.wallK==305.15 && spec.outerRadiusOverD==30);
rows=readtable(tableFile);w=mesh.hybridWall;v=mesh.V;f=mesh.F;
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(isequal(rows.faceID,w.faceID) && isequal(rows.regionID,body.regionID) && all(abs(area-rows.areaM2)<1e-18) && ...
    abs(area'*rows.heatFluxWm2+s.totalOutwardW)<1e-12 && abs(s.energy.relativeDifference)<.005 && s.relativeChange20<.002);
pair=struct('available',true,'baseFrame',id.sha256,'expandedFrame',s.sourceFrameSHA256,'expandedReaction',s.sourceReactionSHA256, ...
    'expandedWallSHA256',s.wallSHA256,'expandedMeshSHA256',s.meshSHA256,'conductance20WK',base.totalOutwardW/(spec.wallK-spec.airK), ...
    'conductance30WK',s.totalOutwardW/(spec.wallK-spec.airK),'baseRadiusOverD',20,'expandedRadiusOverD',30);
% 模块二：同壁温直接比较，两次等功率反演各用本域热导，超出对照范围单独标记。
if isfield(c,'netHeatingPowerW')
    if c.netHeatingPowerW>s.totalOutwardW
        pair.available=false;pair.note=sprintf('30D 参照场支持的整件净功率上限为 %.9g mW。',1000*s.totalOutwardW);
        mesh.thermistorOuterDomain=pair;return;
    end
    [rows,~,~,wall30]=scaleThermistorWallField(rows,c,spec.airK,spec.wallK,s.totalOutwardW);
    wall20=mesh.hybridSolvedWallTemperatureK(1);pair.mode='givenPower';
    pair.scope=sprintf('同一整件净功率 %.9g mW，两域各用本域热导求壁温。',1000*c.netHeatingPowerW);
else
    [rows,~,~,wall30]=scaleThermistorWallField(rows,c,spec.airK,spec.wallK,s.totalOutwardW);
    wall20=c.wallTemperatureK;pair.mode='givenWall';
    pair.scope=sprintf('两域均给定全表面壁温 %.12g K，各用本域 %.6g K 参照场。',wall20,spec.wallK);
end
pair.wall20K=wall20;pair.wall30K=wall30;pair.outward20W=-area'*w.heatFluxWm2;pair.outward30W=-area'*rows.heatFluxWm2;
if strcmp(pair.mode,'givenPower')
    assert(max(abs([pair.outward20W pair.outward30W]-c.netHeatingPowerW))<=1e-11*c.netHeatingPowerW);
end
pair.heat30Wm2=rows.heatFluxWm2;pair.deltaHeatWm2=rows.heatFluxWm2-w.heatFluxWm2;
pair.deltaPressurePa=rows.pressurePa-w.pressurePa;
pair.deltaFirstK=rows.nearWallCellTemperatureK-w.nearWallCellTemperatureK;
% 原面 CSV 的配套来源回执也登记对照帧，便于离开工作台后追查两次求解。
mesh.evidenceIdentity.sourceFrame.outerDomainReference=struct('frameSHA256',s.sourceFrameSHA256, ...
    'reactionSHA256',s.sourceReactionSHA256,'meshSHA256',s.meshSHA256,'wallSHA256',s.wallSHA256,'mode',pair.mode);
for j=1:numel(s.regions)
    region=rows.regionID==j;
    pair.regions(j)=struct('name',s.regions(j).name,'faces',nnz(region), ...
        'outward20W',-area(region)'*w.heatFluxWm2(region),'outward30W',-area(region)'*rows.heatFluxWm2(region));
end
mesh.thermistorOuterDomain=pair;
end
