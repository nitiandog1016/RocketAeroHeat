function testSU2NormalImport(folder,out)
%TESTSU2NORMALIMPORT 回读已终态的新法向层实际启动，不再次运行求解器。
% 模块一：核对真实网格、原生方向、进程与来流；两步场不当收敛解。
% 模块二：将压力/热流/节点温度映射到全部原面，分开核对两种面积负荷。
assert(~isfolder(out));mkdir(out);
s=load(fullfile(folder,'normal_volume.mat'));native=fullfile(folder,'native_input');
receipt=jsondecode(fileread(fullfile(native,'engine_process.json')));
assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut);
assert(~receipt.request.dryRun && all(cellfun(@(v)isequal(v,true),struct2cell(receipt.inputsStable))));
log=fileread(fullfile(native,'engine_stdout.log'));
assert(contains(log,'58560 volume elements.') && ~contains(log,'Using Dummy Geometry'));
assert(contains(log,'All volume elements are correctly oriented.') && contains(log,'All surface elements are correctly oriented.'));
assert(~contains(lower(log),'re-orientation') && ~contains(lower(log),'negative volume'));
letters='XYZ';
for j=1:3
    value=regexp(log,['Velocity-' letters(j) '\|\s*([^|]+)\|'],'tokens','once');assert(~isempty(value));
    assert(abs(str2double(strtrim(value{1}))-s.info.expectedVelocity(j))<max(.01,abs(s.info.expectedVelocity(j))*1e-5));
end
data=readSU2SurfaceResult(fullfile(native,'surface_flow.vtk'));
result=mapSU2WallResult(data,s.projection);result.V=s.surface.V;result.F=s.surface.F;
assert(numel(result.heatFluxRaw)==2448 && all(isfinite([result.pressureRaw;result.heatFluxRaw;result.fluidNodeTemperatureMeanRaw])));
assert(~result.productionQualified && ~result.eligibleForSurfaceCoverage);
history=readtable(fullfile(native,'history.csv'),'VariableNamingRule','preserve');
names=strtrim(history.Properties.VariableNames);hf=history{end,strcmp(names,'HF')};energy=history{end,strcmp(names,'rms[RhoE]')};
assert(abs(result.conservativeHeatLoadPowerRaw/hf-1)<2e-6 && energy>-8);
report=struct('originalFaces',2448,'cells',58560,'actualImportNoReorientation',true,'lastEnergyResidual',energy, ...
    'engineHeatPowerRaw',hf,'physicalAreaHeatPowerRaw',result.targetHeatPowerRaw, ...
    'conservativeHeatLoadPowerRaw',result.conservativeHeatLoadPowerRaw,'physicsQualified',false,'newQualifiedArea',0);
save(fullfile(out,'mapped_startup.mat'),'result','report','-v7');
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Actual normal-grid import: native orientations accepted; all 2448 original faces mapped; startup only.\n');
end
