function report=testSchaeferThermalContract(out)
%TESTSCHAEFERTHERMALCONTRACT 解析几何、总静态/单位Re换算及热边界拒绝测试，不启动数值引擎。
assert(~isfolder(out));mkdir(out);[s,m]=buildSchaeferFergusonBody();c=schaeferFergusonCondition();
assert(abs(m.facetedVolume/(m.exactSmoothVolume*m.facetedVolumeRatio)-1)<2e-13);
assert(abs(m.totalLength-17.57*.0254)<1e-15 && abs(m.actualNoseAngleDeg-15)>.1);
assert(isequal(sort(unique(m.faceRegion)),sort({'base';'cone';'cylinder';'flare'})));
assert(all(s.normal(strcmp(m.faceRegion,'base'),1)<-.999));
fs=c.atmosphere;U=c.Mach*fs.speedOfSound;
assert(abs(fs.density*U/fs.dynamicViscosity*.3048-1.6e6)<1e-8);
assert(abs(fs.temperature+.5*U^2/(1.4*fs.R/.4)-710*5/9)<1e-12);
iso=resolveSU2ThermalBoundary(s,struct('mode','isothermal','temperature',300));
ad=resolveSU2ThermalBoundary(s,struct('mode','adiabatic'));
assert(all(iso.givenTemperatureByFace==300) && all(isnan(ad.givenTemperatureByFace)) && all(ad.givenHeatFluxByFace==0));
profile=struct('mode','axialTemperature','x',[0 m.totalLength],'temperature',[80 140],'source','synthetic linear contract test; not experimental');
sp=resolveSU2ThermalBoundary(s,profile);
assert(max(abs(sp.givenTemperatureByFace-(80+60*s.center(:,1)/m.totalLength)))<1e-12 && ~sp.stockExecutableSupported);
bad={struct('mode','adiabatic','temperature',300),struct('mode','isothermal','temperature',-1), ...
    struct('mode','axialTemperature','x',[.1 .2],'temperature',[80 90]), ...
    struct('mode','axialTemperature','x',[0 0 m.totalLength],'temperature',[80 90 100])};
for j=1:numel(bad),rejected=false;try,resolveSU2ThermalBoundary(s,bad{j});catch,rejected=true;end;assert(rejected);end
% 必须在建目录前拒绝尚无实际调用器的空间Tw，不能变成等温成功运行。
dummy=struct('kind','complete-rocket-three-dimensional-su2-volume','physicalSolution',false, ...
    'sourceVertices',s.V,'sourceFaces',s.F);folder=fullfile(out,'unsupported');
rejected=false;try,writeSU2RocketCase(dummy,folder,struct('thermalBoundary',profile));catch ex,rejected=strcmp(ex.identifier,'RocketAeroHeat:SU2:SpatialThermalAdapter');end
assert(rejected && ~isfolder(folder));
% 语法检查只对本次改动文件执行；警告保留但不自动伪称为失败或物理通过。
files={'buildSchaeferFergusonBody.m','schaeferFergusonCondition.m','resolveSU2ThermalBoundary.m', ...
    'attachSU2ThermalBoundary.m','writeSU2RocketCase.m','createSU2ResearchView.m','rocketAeroHeatWorkbench.m','readSU2TimeoutCheckpoint.m'};
syntax=struct;for k=1:numel(files),issues=checkcode(files{k},'-id');assert(~any(strcmp({issues.id},'PARSE')));syntax.(matlab.lang.makeValidName(files{k}))=issues;end
report=struct('geometry',m,'condition',c,'syntax',syntax,'thermalContractPassed',true, ...
    'missingSpatialAdapterRejectedBeforeWrite',rejected,'solverInvoked',false,'newQualifiedArea',0);
save(fullfile(out,'geometry.mat'),'s','m','c','-v7');write(fullfile(out,'report.json'),report);
fprintf('Classic closed geometry, analytical volume, exact SI conversion and thermal contract passed.\n');
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
