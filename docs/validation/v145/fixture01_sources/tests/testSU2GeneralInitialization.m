function report=testSU2GeneralInitialization(out)
%TESTSU2GENERALINITIALIZATION 小型混合域两步实际SST接口，粗网格不授予质量或热资格。
% 模块一：复用已验证的小型棱柱/金字塔/四面体，不重复生成几何。
% 模块二：新壁距初态经同一后台准备、原生进程及同阶段热审核，保存真实终态。
assert(~isfolder(out));mkdir(out);d=load('outputs/v135_validation/fixture/mesh_input.mat');
s=struct('V',d.V,'F',d.F);v=readSU2HybridVolume('outputs/v135_validation/fixture',s);
folder=fullfile(out,'case');flight=struct('Mach',5,'temperature',226.5,'pressure',1171.8,'alpha',20,'beta',0, ...
    'wallTemperature',300,'iterations',2,'outputEvery',1,'namedCheckpoints',true);
writeSU2RocketCase(v,folder,flight);
spec=struct('model','SST-2003m','mode','wallDistance','thicknessM',.001,'intensity',.01,'viscosityRatio',10, ...
    'Prt',.9,'provenance','v145 explicit fully turbulent startup assumption; not measured transition');
job=embeddedSU2Job('start',folder,struct('threads',2,'maximumSeconds',120,'pythonExecutable','D:/anacondapython/python.exe','turbulenceInitialization',spec));
t=tic;while ~job.terminal && toc(t)<440,pause(.2);[job,r]=embeddedSU2Job('poll',job);end
assert(job.terminal && r.exitCode==0 && job.pipelineReceipt.exitCode==0,'小型SST流程未成功。');
a=readSU2SSTStateAudit(folder,r);assert(a.checkpointLocalIteration==1);
report=struct('passed',true,'nodes',size(v.vertices,1),'cells',numel(v.cells),'thermalAudit',a, ...
    'geometryQualityQualified',false,'physicsQualified',false,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
