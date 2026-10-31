function report=testSU2HybridVolume(out,meshFolder,inputFile)
%TESTSU2HYBRIDVOLUME 对原面混合网格独立审核并执行两步真实NS导入，不宣称收敛。
% 模块一：复用已有网格，不重复构造；核对界面、原面与完整外域体积。
% 模块二：导出并实际运行固定SU2，检查原壁面的p/q/流体温度及负荷口径。
assert(~isfolder(out));mkdir(out);d=load(inputFile,'V','F');s=struct('V',d.V,'F',d.F);
t=tic;volume=readSU2HybridVolume(meshFolder,s);projection=buildSU2WallProjection(volume,s);
assert(max(abs(full(sum(projection.weights,2))-1))<1e-14);
assert(nnz(projection.weights)==3*size(s.F,1));
flight=struct('Mach',5,'temperature',226.65,'pressure',1171.87,'wallTemperature',300, ...
    'alpha',20,'beta',0,'iterations',2,'outputEvery',2);
folder=fullfile(out,'engine');manifest=writeSU2RocketCase(volume,folder,flight);
receipt=runEmbeddedSU2Case(folder,struct('pythonExecutable','D:/anacondapython/python.exe','threads',2,'maximumSeconds',120));
data=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));r=mapSU2WallResult(data,projection);
r.V=s.V;r.F=s.F;r.flight=flight;r.wallTemperatureInput=300;r.receipt=receipt;r.resultDirectory=folder;
r.convergenceAudit=auditSU2ResultConvergence(folder,r);
assert(numel(r.pressureRaw)==size(s.F,1) && all(isfinite([r.pressureRaw;r.heatFluxRaw;r.fluidNodeTemperatureMeanRaw])));
assert(all(r.pressureRaw>0) && all(r.fluidNodeTemperatureMeanRaw>0) && ~r.productionQualified);
text=fileread(fullfile(folder,'engine_stdout.log'));
report=struct('originalFaces',size(s.F,1),'meshAudit',volume.topology,'meshQuality',volume.meshConstruction, ...
    'manifest',manifest,'receipt',receipt,'mappingAreaError',projection.maximumTargetAreaRelativeError, ...
    'physicalAreaHeatPower',r.targetHeatPowerRaw,'nativeNodeHeatPower',r.reconstructedEngineHeatPowerRaw, ...
    'solverOrientationMessages',{regexp(text,'[^\r\n]*[Oo]rient[^\r\n]*','match')}, ...
    'elapsedSeconds',toc(t),'newQualifiedArea',0,'physicsQualified',false);
save(fullfile(out,'mapped_result.mat'),'r','report','-v7');
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Hybrid original faces=%d, cells=%d, positive complete geometry, two actual NS steps; no physical qualification.\n',size(s.F,1),numel(volume.cells));
end
