function report=testSU2PyramidInterface(out)
%TESTSU2PYRAMIDINTERFACE 固定小型混合网格的原生方向、场量、原面投影检查。
% 模块一：读取完整中位体积和相容四边接口，直接拒绝翻转金字塔。
% 模块二：两步实际NS仅证明导入与字段链路，粗网格质量未过仍明确保留。
assert(~isfolder(out));mkdir(out);folder='outputs/v135_validation/fixture';d=load(fullfile(folder,'mesh_input.mat'));
s=struct('V',d.V,'F',d.F);volume=readSU2HybridVolume(folder,s);
assert(sum(volume.cellTypes==14)==8 && sum(volume.cellTypes==13)==16 && sum(volume.cellTypes==10)==68);
assert(volume.pyramidGeometry.minimumScaledJacobian>0 && ~volume.meshConstruction.qualityGatePassed);
bad=vertcat(volume.cells{volume.cellTypes==14});bad=bad(:,[1 4 3 2 5]);rejected=false;
try,su2PyramidGeometry(volume.vertices,bad);catch,rejected=true;end
assert(rejected);projection=buildSU2WallProjection(volume,s);
flight=struct('Mach',5,'temperature',226.65,'pressure',1171.87,'wallTemperature',300,'alpha',20,'beta',0,'iterations',2,'outputEvery',2);
engine=fullfile(out,'engine');manifest=writeSU2RocketCase(volume,engine,flight);assert(manifest.pyramids==8);
receipt=runEmbeddedSU2Case(engine,struct('pythonExecutable','D:/anacondapython/python.exe','threads',2,'maximumSeconds',120));
data=readSU2SurfaceResult(fullfile(engine,'surface_flow.vtk'));result=mapSU2WallResult(data,projection);
assert(numel(result.pressureRaw)==8 && all(isfinite([result.pressureRaw;result.heatFluxRaw;result.fluidNodeTemperatureMeanRaw])));
assert(all(result.pressureRaw>0) && all(result.fluidNodeTemperatureMeanRaw>0) && ~result.productionQualified);
log=fileread(fullfile(engine,'engine_stdout.log'));
assert(isempty(regexp(log,'re-orientation of [1-9][0-9]* PYRAMID volume elements','once')));
assert(~contains(log,'volume elements are distorted'));
save(fullfile(out,'interface.mat'),'volume','result','projection','manifest','receipt','-v7');
report=struct('originalFaces',8,'pyramids',8,'cellCount',92,'geometry',volume.topology, ...
    'pyramidGeometry',volume.pyramidGeometry,'invertedPyramidRejected',rejected,'nativeExitCode',receipt.exitCode, ...
    'nativeSeconds',receipt.elapsedSeconds,'pyramidReorientation',false,'finiteOriginalFields',true, ...
    'qualityGatePassed',false,'nativeFlowRuns',1,'realRocketFlowRuns',0,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
