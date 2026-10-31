function report=reviewSU2CoolPropCapability(out)
%REVIEWSU2COOLPROPCAPABILITY 复用真实引擎结果核对MATLAB读取，不重复流场求解。
% 模块一：真实守恒场与同库EOS输出绑定；确认保存迭代、压力及流体温度。
% 模块二：与既有NASA7作同温同速的焓增量对照，明确零点和模型差别。
% 模块三：旧定比热结果只读回归，保护工作台已有结果；不宣称新主流程UI验收。
assert(~isfile(out));folder='outputs/v121_validation/thermal20';
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
wall=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
w=jsondecode(fileread(fullfile(folder,'wall_audit.json')));
r=struct('receipt',receipt,'reconstructedEngineHeatPowerRaw',w.nativeIntegratedHeatPower);
a=auditSU2ResultConvergence(folder,r);assert(a.checkpointIdentityConfirmed,a.reason);
assert(a.checkpointIteration==19 && ~a.checkpointResidualPassed && ~a.productionQualified);
info=describeSU2Runtime(receipt);assert(strcmp(info.mode,'compressible-heos-air') && contains(info.method,'CoolProp'));
assert(contains(info.physicsCaption,'HEOS Air') && ~contains(info.physicsCaption,'ideal gas'));
assert(isequal(a.unavailableDiagnosticFields,{'Buffet'}));
assert(numel(wall.pressureRaw)==w.nodes && all(isfinite(wall.heatFluxRaw)));
p=readmatrix('outputs/v121_validation/properties/PT_output.csv','NumHeaderLines',1);
[nasaCp,nasaH]=airThermalProperties(p(:,3));
f=jsondecode(fileread('outputs/v121_validation/property_audit.json'));
freestream=readmatrix('outputs/v121_validation/properties/freestream_output.csv','NumHeaderLines',1);
[~,hinf]=airThermalProperties(freestream(3));
total=hinf+f.originalVelocityPreservedMps^2/2;
nasaT0=fzero(@(T)enthalpy(T)-total,[freestream(3) 2000]);
old=load('outputs/v119_validation/review01/result.mat','r');
% 旧结果路径由旧回执确定，不猜另一个运行目录。
oldAudit=auditSU2ResultConvergence(old.r.receipt.request.workDirectory,old.r);
assert(oldAudit.checkpointIdentityConfirmed,oldAudit.reason);
oldInfo=describeSU2Runtime(old.r.receipt);
assert(oldAudit.cumulativeCheckpointIteration==599 && strcmp(oldInfo.flavor,'official'));
% 明确标识的合成反例只改测试目录：必要残差缺失不能像未定义Buffet那样放行。
fixture=fullfile(fileparts(out),'nonfinite_residual_fixture');assert(~isfolder(fixture));mkdir(fixture);
for name={'case.cfg','surface_flow.vtk'},copyfile(fullfile(folder,name{1}),fullfile(fixture,name{1}));end
text=fileread(fullfile(folder,'history.csv'));lines=splitlines(string(text));row=split(lines(2),',');row(4)="nan";lines(2)=join(row,',');
writeText(fullfile(fixture,'history.csv'),char(join(lines,newline)));
raw=fileread(fullfile(folder,'engine_process.json'));raw=strrep(raw,receipt.outputSHA256.history_csv,sha(fullfile(fixture,'history.csv')));
raw=['{"syntheticTestFixture":true,"note":"no solver run; injected nonfinite density residual",' raw(2:end)];
writeText(fullfile(fixture,'engine_process.json'),raw);bad=r;bad.receipt=jsondecode(raw);
rejected=auditSU2ResultConvergence(fixture,bad);
assert(strcmp(rejected.status,'evidence-rejected') && ~rejected.checkpointIdentityConfirmed && contains(rejected.reason,'必要字段'));
report=struct('audit',a,'runtime',info,'oldConstantCpSnapshotPassed',true, ...
    'temperatureK',p(:,3),'pressurePa',p(:,2),'CoolPropCpJkgK',p(:,7),'NASA7CpJkgK',nasaCp, ...
    'CoolPropEnthalpyRiseJkg',p(:,6)-freestream(6),'NASA7EnthalpyRiseJkg',nasaH-hinf, ...
    'NASA7TotalTemperatureK',nasaT0,'CoolPropTotalTemperatureK',f.isentropicStagnationTemperatureK, ...
    'nonfiniteResidualRejected',rejected,'comparisonOnly',true,'newMainUIClaimed',false,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(out,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('CoolProp MATLAB readback and old snapshot passed; current startup field is not qualified.\n');
end
function h=enthalpy(T)
[~,h]=airThermalProperties(T);
end
function writeText(file,text)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));bytes=fread(fid,inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
