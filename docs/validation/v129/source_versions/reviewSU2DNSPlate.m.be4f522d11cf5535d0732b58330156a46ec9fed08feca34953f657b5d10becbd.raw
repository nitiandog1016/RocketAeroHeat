function report=reviewSU2DNSPlate(out,comparisonFile)
%REVIEWSU2DNSPLATE 只读既有终态，核对MATLAB字段与Python公开基准结果。
% 模块一：复用原生算例，绝不重新启动SU2或改写历史数值。
% 模块二：区分给定Tw、流体节点温度和未知求解Tw，保留SST算法身份。
folder=fullfile(out,'plate');target=fullfile(out,'matlab_receipt.mat');assert(~isfile(target));
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut);
info=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
data=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
assert(all(isfinite(data.pressureRaw)) && all(isfinite(data.heatFluxRaw)));
report=jsondecode(fileread(comparisonFile));
assert(report.interfaceAudit.checkpointLocalIteration==9999 && report.stationMatched);
assert(report.newQualifiedRocketAreaM2==0 && ~report.physicsQualified);
model=describeSU2Runtime(receipt);assert(strcmp(model.mode,'sst-2003m'));
prescribedWallTemperature=info.options.wallTemperature;assert(prescribedWallTemperature==97.5);
solvedWallTemperature=NaN(size(data.temperatureRaw));
assert(contains(data.temperatureIdentity,'fluid-node-state'));
assert(~report.exactDNSReplication && ~info.publicBenchmark.exactDNSReplication);
save(target,'info','receipt','data','report','model','prescribedWallTemperature','solvedWallTemperature','-v7');
fprintf('Readonly DNS comparison: matched=%d, residual gate=%d; original 9999 retained, rocket qualification=0.\n', ...
    report.stationMatched,report.residualConverged);
end
