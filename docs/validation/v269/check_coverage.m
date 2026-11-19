function check_coverage(folder)
%CHECK_COVERAGE 在新的目录保存接口测试与三份既有计算结果的比对，所有大场只读。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root,fullfile(root,'tests'));
assert(isfolder(folder)&&~isfile(fullfile(folder,'result.json')),'结果文件已存在或目录无效。');
control=testComputedWallCoverage();
readJSON=@(name)jsondecode(fileread(fullfile(root,name)));
% 一、冻结自绘算例。原 MAT 经工作台导出，再同旧版独立原面分类及 v266 功率比对。
raw=load(fullfile(root,'outputs/v104_validation/baseline02/current_case.mat'),'current');
source=raw.current;clear raw;
p=source.pressure.fieldStatus(:)==1;q=source.heating.eligibleFace(:);
main=summarizeComputedWall(source.localGeometry.area,p,q,source.heating.heatFlux);
viewer=buildAerothermalViewerData(source);
assert(strcmp(main.coverageMessage,viewer.wallCoverage.coverageMessage) && contains(viewer.caseNote,main.coverageMessage));
joint=readJSON('docs/validation/v263/evidence/saved_surfaces.json');
power=readJSON('docs/validation/v263/evidence/snapshot_power/summary.json');
assert(main.faceCount==joint.faceCount && main.jointComputedFaceCount==joint.jointGroups.bothComputed.faces);
assert(abs(main.totalArea-joint.originalAreaM2)<1e-11 && abs(main.jointComputedAreaM2-joint.jointGroups.bothComputed.areaM2)<1e-11);
assert(abs(main.jointMissingAreaM2-(joint.originalAreaM2-joint.jointGroups.bothComputed.areaM2))<1e-11);
assert(abs(main.heatComputedAreaM2-power.mainSnapshot.calculatedAreaM2)<1e-11 && ~main.fullSurfaceComputed);
assert(abs(main.calculatedHeatPowerW/power.mainSnapshot.calculatedPowerW-1)<1e-11 && isnan(main.wholeSurfaceHeatPowerW));
clear source viewer;

% 二、以原工况、原网格和带逐变量方法的原面 CSV 还原微球，导出可直接打开的页面。
raw=strtrim(fileread(fullfile(root,'outputs/v268_validation/sphere03/fine/workbench/case-data.js')));
prefix='window.ROCKET_CASE=';assert(startsWith(raw,prefix)&&endsWith(raw,';'));
saved=jsondecode(raw(numel(prefix)+1:end-1));mesh=struct('V',saved.vertices,'F',saved.faces+1);
p=saved.fields(strcmp({saved.fields.key},'pressure'));q=saved.fields(strcmp({saved.fields.key},'heatFlux'));
assert(all(p.status==1 & q.status==1) && strcmp(p.unit,'kPa') && strcmp(q.unit,'kW/m²'));
mesh.hybridWall=readtable(fullfile(root,'outputs/v268_validation/sphere03/fine/wall.csv'),'TextType','string','Encoding','UTF-8');
assert(isequal(mesh.hybridWall.faceID,(0:numel(q.status)-1)'));
mesh.hybridCondition=saved.resultCondition.physics;
mesh.hybridWallTemperatureK=saved.fields(strcmp({saved.fields.key},'wallInput')).values;
mesh.localGeometry.region=saved.region;mesh.evidenceIdentity=saved.evidenceIdentity;
assert(~isfolder(fullfile(folder,'workbench')),'工作台已经存在，请使用新目录。');
page=exportAerothermalWeb(mesh,fullfile(folder,'workbench'));
exported=strtrim(fileread(fullfile(fileparts(page),'case-data.js')));
assert(startsWith(exported,prefix) && endsWith(exported,';'));
shown=jsondecode(exported(numel(prefix)+1:end-1));report=shown.wallCoverage;
assert(isequal(saved.faces,shown.faces) && isequal(saved.vertices,shown.vertices));
assert(isequal({saved.fields.key},{shown.fields.key}) && startsWith(shown.caseNote,report.coverageMessage));
for k=1:numel(saved.fields)
    old=saved.fields(k);new=shown.fields(k);known=old.status~=0;
    assert(isequal(old.status,new.status) && isequal(isnan(old.values),isnan(new.values)));
    assert(all(abs(old.values(known)-new.values(known))<1e-10));
    assert(isequal(old.methodNames,new.methodNames) && isequal(old.methodIndex,new.methodIndex));
end
sphere=readJSON('docs/validation/v268/evidence/audit_result.json');
assert(report.fullSurfaceComputed && report.jointComputedFaceCount==sphere.originalFaces && report.jointMissingAreaM2==0);
assert(abs(report.totalArea-sphere.originalAreaM2)<1e-15 && abs(report.wholeSurfaceHeatPowerW/sphere.originalTrianglePowerW-1)<1e-11);

% 三、风洞功率按源曲面计算。原工作台曲面积规则仍可复现原环带功率。
tunnel=readJSON('outputs/v263_validation/tmr_resolved_view02/workbench/results.json');
ledger=buildAerothermalEvidenceLedger(tunnel);
q=ledger.fields(strcmp({ledger.fields.key},'heatFlux'));curved=q.areaSummary.heatPower;
assert(strcmp(curved.areaBasis,'source-curved') && strcmp(curved.status,'complete'));
assert(abs(curved.wholeSurfacePowerW/power.windTunnelRingIntegralW-1)<1e-11);
result=struct('interfaceTests',control,'mainArchive',main,'sphereReplay',report, ...
    'tunnelHeatPower',curved,'newMainCoverageM2',0,'newValidationAreaM2',0,'browserInteractionTested',false);
fid=fopen(fullfile(folder,'result.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
end
