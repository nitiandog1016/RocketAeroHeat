function testPhaseNonsimilarGuardsV252(folder)
%TESTPHASENONSIMILARGUARDSV252 一站修复回归与保存数据的映射完整性；不重跑三条链。
% 模块一：仅重算已绑定的第2站，核对晚期可用标记修复未改变数值。
d=load('outputs/v104_validation/baseline02/current_case.mat','current');mesh=d.current;
chain=jsondecode(fileread('outputs/v252_validation/matlab01/base.json'));
if ~iscell(chain.profiles),chain.profiles=num2cell(chain.profiles);end
a=chain.profiles{1};b=chain.profiles{2};r=solvePhaseAwareNonsimilarStation(b.edge,b.station,b.wallTemperature,a);
assert(r.available,'RAH:PhaseMarch:GuardRegression','修复后第2站拒绝。');
keys={'heatFlux','wallShearStress','delta','thetaPhysical','energyThickness','enthalpyDefectThickness'};difference=0;
for k=keys,difference=max(difference,abs(r.(k{1})-b.(k{1}))/max(1,abs(b.(k{1}))));end
assert(difference<1e-12);
% 模块二：从原准入输入恢复外缘；保存JSON的浮点回读可能不能恢复原MATLAB逐字节哈希。
e=jsondecode(fileread('outputs/v226_validation/matlab02/rocket_metadata.json'));
v=readmatrix('outputs/v226_validation/matlab02/rocket_faces.csv');
e.pressurePa=v(:,3);e.temperatureK=v(:,4);e.speedMps=v(:,5);e.densityKgm3=v(:,6);e.velocityMps=v(:,7:9);
e.Mach=v(:,10);e.recoveryTemperatureK=v(:,11);e.polarAngleRad=v(:,12);e.arcLengthM=v(:,13);
e.validFace=logical(v(:,14));e.reasonCode=uint8(v(:,15));e.faceAreaM2=v(:,2);
% 这是显式JSON重新序列化身份，仅供映射接口回归，不追改原链或原源摘要。
originalSHA=chain.sourceSHA256;
chain.sourceSHA256=phaseAwarePayloadSHA(chain.sourceText,jsonencode(chain.profiles));
actual=mapPhaseAwareMarchToSurface(mesh,e,chain);
old=jsondecode(fileread('outputs/v252_validation/matlab01/surface_base.json'));
assert(isequal(actual.validFace,old.validFace)&&max(abs(actual.heatFlux(actual.validFace)-old.heatFlux(actual.validFace)))<1e-9);
rejections=0;
bad=chain;bad.stations(2).thetaRad=bad.stations(2).thetaRad+.001;reject(bad);
bad=chain;bad.available(2)=false;reject(bad);
bad=chain;bad.wallTemperatureHistory(2)=301;reject(bad);
bad=chain;bad.edges(2).temperatureK=bad.edges(2).temperatureK+1;reject(bad);
bad=chain;bad.domain.radiusM=bad.domain.radiusM*2;reject(bad);
assert(rejections==5);
report=struct('passed',true,'singleStationRecomputed',2,'maximumScaledDifference',difference, ...
    'mappingRejections',rejections,'originalChainSHA256',originalSHA,'reserializedTestChainSHA256',chain.sourceSHA256, ...
    'originalResultsModified',false,'fullThreeChainsRerun',false,'physicalQualification',false);
fid=fopen(fullfile(folder,'report.json'),'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));disp(report);
    function reject(bad)
        rejected=false;try,mapPhaseAwareMarchToSurface(mesh,e,bad);catch err,rejected=strcmp(err.identifier,'RAH:PhaseMarch:Mapping');end
        assert(rejected);rejections=rejections+1;
    end
end
