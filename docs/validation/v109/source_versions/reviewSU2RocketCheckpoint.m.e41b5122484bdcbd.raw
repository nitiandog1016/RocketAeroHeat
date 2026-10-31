function reviewSU2RocketCheckpoint(folder,out)
%REVIEWSU2ROCKETCHECKPOINT 复用预算终态保存的实际检查点，映射全部原面，不追加求解。
% 模块一：超时检查点与实际进程末次监测分开；它不等于最后一次迭代状态。
% 模块二：原面和两种负荷口径保持，未收敛身份不因文件可读而升级。
assert(~isfolder(out));mkdir(out);
s=load('outputs/v108_validation/normal01/normal_volume.mat','surface','projection');
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));assert(strcmp(receipt.status,'terminal'));
assert(all(cellfun(@(v)isequal(v,true),struct2cell(receipt.inputsStable))));
data=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));result=mapSU2WallResult(data,s.projection);
result.V=s.surface.V;result.F=s.surface.F;result.receipt=receipt;
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));result.flight=manifest.flight;
result.wallTemperatureInput=manifest.flight.wallTemperature;
result.resultStatus='saved-checkpoint-only; terminal-history-is-a-separate-time; physics-unqualified';
assert(all(isfinite(result.pressureRaw)) && all(isfinite(result.heatFluxRaw)) && ~result.productionQualified);
save(fullfile(out,'checkpoint_original_faces.mat'),'result','-v7');
report=struct('originalFaces',numel(result.pressureRaw),'originalArea',sum(result.faceArea), ...
    'physicalAreaHeatPower',result.targetHeatPowerRaw,'nativeHeatLoadPower',result.conservativeHeatLoadPowerRaw, ...
    'terminalTimedOut',receipt.timedOut,'checkpointIterationMustBeAuditedSeparately',true,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Saved complete-rocket checkpoint mapped to original faces; terminal history is separate; unqualified.\n');
end
