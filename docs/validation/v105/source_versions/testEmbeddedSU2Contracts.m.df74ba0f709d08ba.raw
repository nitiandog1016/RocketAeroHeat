function testEmbeddedSU2Contracts(source,out,pythonExecutable)
%TESTEMBEDDEDSU2CONTRACTS 使用真实壁面输出验证读取，并实际检查失败/超时终态。
% 模块一：原数据解析、缺失热流拒绝、负热流保留及截断文件拒绝。
% 模块二：SU2配置失败和1秒预算超时由真实进程执行，不计物理算例通过。
% 模块三：临时输入保留在新目录，禁止覆写已有证据；温度身份不与给定Tw混淆。
if nargin<3,pythonExecutable='python';end
assert(~exist(out,'dir'));mkdir(out);raw=fileread(source);
data=readSU2SurfaceResult(source);assert(numel(data.pressureRaw)==65);
assert(strcmp(data.temperatureIdentity,'raw-fluid-node-state; distinct-from-prescribed-wall-temperature'));
assert(~data.productionQualified && ~data.heatFluxSignConfirmed);
missing=strrep(raw,'SCALARS Heat_Flux double 1','SCALARS Missing_Flux double 1');
file=fullfile(out,'missing.vtk');write(file,missing);expectFailure(@()readSU2SurfaceResult(file),'RocketAeroHeat:SU2:Fields');
file=fullfile(out,'truncated.vtk');write(file,raw(1:100));expectFailure(@()readSU2SurfaceResult(file),'RocketAeroHeat:SU2:VTK');
negative=regexprep(raw,'(SCALARS Heat_Flux double 1\s+LOOKUP_TABLE default\s+)[^\s]+','$1-123.5','once');
file=fullfile(out,'negative.vtk');write(file,negative);other=readSU2SurfaceResult(file);assert(other.heatFluxRaw(1)==-123.5);
assert(isequal(other.heatFluxRaw(2:end),data.heatFluxRaw(2:end)));
bad=fullfile(out,'invalid_config');buildSU2FlatPlateCase(bad,struct('nx',24,'ny',24));
fid=fopen(fullfile(bad,'case.cfg'),'a');fprintf(fid,'RAH_INTENTIONAL_INVALID_OPTION= YES\n');fclose(fid);
expectFailure(@()runEmbeddedSU2Case(bad,struct('pythonExecutable',pythonExecutable,'maximumSeconds',30)), ...
    'RocketAeroHeat:SU2:Failed');
invalid=jsondecode(fileread(fullfile(bad,'engine_process.json')));
assert(strcmp(invalid.status,'terminal') && invalid.exitCode~=0 && ~invalid.timedOut);
limited=fullfile(out,'time_limit');buildSU2FlatPlateCase(limited,struct('iterations',1000000));
expectFailure(@()runEmbeddedSU2Case(limited,struct('pythonExecutable',pythonExecutable,'maximumSeconds',1)), ...
    'RocketAeroHeat:SU2:Failed');
timeout=jsondecode(fileread(fullfile(limited,'engine_process.json')));
assert(strcmp(timeout.status,'terminal') && timeout.timedOut);
report=struct('parserNodes',65,'missingHeatRejected',true,'truncatedRejected',true, ...
    'negativeHeatPreserved',true,'temperatureIdentityChecked',true, ...
    'realInvalidConfigTerminal',invalid.exitCode,'realTimeoutTerminal',timeout.exitCode, ...
    'physicsQualification',false);
write(fullfile(out,'report.json'),jsonencode(report));
fprintf('SU2 parser and actual failure/timeout lifecycle checks passed; no physical qualification.\n');
end

function write(file,text)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
function expectFailure(action,id)
try
    action();
catch exception
    assert(strcmp(exception.identifier,id),'意外错误类型：%s',exception.identifier);return;
end
error('RocketAeroHeat:SU2:Test','应拒绝的输入未报错：%s',id);
end
