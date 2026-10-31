function testSU2NativeCheckpoint(out)
%TESTSU2NATIVECHECKPOINT 复用实际具名结果，验证时刻身份、全原面映射和配置输出。
% 不启动新求解；缺对/篡改为显式合成夹具，不冒充真实进程。
assert(~isfolder(out));mkdir(out);s=load('outputs/v108_validation/normal01/normal_volume.mat','projection','volume');
names={'first_probe01','control01','cfl01'};reports=cell(3,1);
for j=1:3
    folder=fullfile('outputs/v111_validation',names{j});data=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
    result=mapSU2WallResult(data,s.projection);result.receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
    a=auditSU2ResultConvergence(folder,result);assert(a.checkpointIdentityConfirmed,a.reason);
    assert(strcmp(a.checkpointIdentityMethod,'native-indexed-pair-SHA-and-EOS') && ~a.checkpointResidualPassed && ~a.productionQualified);
    expected=599;if j==1,expected=0;end;assert(a.checkpointIteration==expected && a.nativeCheckpointIdentity.localIteration==expected);
    reports{j}=a;save(fullfile(out,[names{j} '_mapped.mat']),'result','a','-v7');
end
% 小型拒绝夹具只复制配置/回执与壁面，缺少体场时必须拒绝，不回退到HF推时刻。
folder='outputs/v111_validation/control01';fixture=fullfile(out,'missing_pair_fixture');mkdir(fixture);
copyfile(fullfile(folder,'surface_flow.vtk'),fullfile(fixture,'surface_flow.vtk'));
rr=jsondecode(fileread(fullfile(folder,'engine_process.json')));cfg=fileread(fullfile(folder,'case.cfg'));
try,identifySU2NativeCheckpoint(fixture,rr,cfg);error('unexpected-pass');catch ex,assert(~strcmp(ex.message,'unexpected-pass'));end
% 新输出选项显式启用；默认保持覆盖写，避免无意累积完整体场。
plain=fullfile(out,'writer_default');named=fullfile(out,'writer_named');
writeSU2RocketCase(s.volume,plain,struct('iterations',2));writeSU2RocketCase(s.volume,named,struct('iterations',2,'namedCheckpoints',true));
assert(contains(fileread(fullfile(plain,'case.cfg')),'WRT_RESTART_OVERWRITE= YES'));
assert(contains(fileread(fullfile(named,'case.cfg')),'WRT_RESTART_OVERWRITE= NO'));
report=struct('actualCases',{names},'audits',{reports},'missingPairRejected',true,'explicitWriterOptionPassed',true, ...
    'defaultDoesNotAccumulateNamedFiles',true,'solverInvoked',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Native checkpoint identity and original-face maps verified; no solver rerun.\n');
end
