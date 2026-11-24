function [mesh,stations,s]=loadShortSupportSurface(g,job,resolution,radiusD,target)
% 模块一：核对同一几何、温界、收敛端点和五区交叠文件，从存档或新算目录取数。
if nargin<2,job='';end
if nargin<3,resolution='coarse';end
if nargin<4,radiusD=14;end
if nargin<5,target=.25;end
root=fileparts(mfilename('fullpath'));cases=fullfile(root,'cases','short_sphere');
timed=ismember(resolution,{'timefull','timehalf'});
transferName='supportedTransfer65.npz';
powerMode=strcmp(resolution,'power');
if timed || powerMode || strcmp(resolution,'angular'),transferName='supportedAngularTransfer65.npz';end
if powerMode
    [cal,archive]=readShortSupportPowerCase(job);
elseif timed
    branch=resolution(5:end);
    if isempty(job)
        pair=readShortSupportTimePair(cases,radiusD,target);
        folder=cases;if radiusD==20,folder=fullfile(cases,'outer20');end
        if target==1,folder=fullfile(cases,'later1');end
        archive=fullfile(folder,['supportedTime65' upper(branch(1)) branch(2:end) '.zip']);
    else
        pair=jsondecode(fileread(fullfile(job,'export','pair.json')));archive=fullfile(job,'export',branch,'surface.zip');
        ack=jsondecode(fileread(fullfile(job,'completed.json')));
        assert(strcmp(cm011FileSHA256(fullfile(job,'export','pair.json')),ack.pairSHA256));
        validateShortTimeChain(pair,target);
    end
    assert(strcmp(cm011FileSHA256(archive),pair.([branch 'Archive'])));
elseif nargin<2 || isempty(job)
    if strcmp(resolution,'angular')
        pair=readShortSupportAngularPair(cases);archive=fullfile(cases,'supportedAngularSurface65.zip');
        assert(strcmp(cm011FileSHA256(archive),pair.fineArchiveSHA256));
    elseif strcmp(resolution,'fine')
        pair=readShortSupportGridPair(cases);archive=fullfile(cases,'supportedFine65.zip');
        assert(strcmp(cm011FileSHA256(archive),pair.fineArchiveSHA256));
    else
        archive=fullfile(cases,'supportedSurface65.zip');
        assert(strcmp(cm011FileSHA256(archive),'2e5b065e3b7bd7f15dc46a22bdc6f5a041bf57ebb9a158c31cb46426953c772b'));
    end
else
    ack=jsondecode(fileread(fullfile(job,'completed.json')));archive=fullfile(job,'surface','surface.zip');
    assert(strcmp(cm011FileSHA256(archive),ack.archiveSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(ack.source,'thermal.json')),ack.thermalSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(ack.source,'early','initial_flow.csv')),ack.frameSHA256));
end
temp=tempname;mkdir(temp);clean=onCleanup(@()discard(temp)); %#ok<NASGU>
unzip(archive,temp);s=jsondecode(fileread(fullfile(temp,'summary.json')));
assert(strcmp(cm011FileSHA256(fullfile(cases,'test65.json')),s.caseSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(cases,'supported65.zip')),s.bodySHA256) && ...
    strcmp(cm011FileSHA256(fullfile(cases,transferName)),s.transferSHA256));
if timed
    e=s.energy;assert(strcmp(s.schema,'short65-supported-bdf1-v1') && strcmp(s.branch,branch) && ...
        strcmp(s.frameSHA256,pair.([branch 'Frame'])) && strcmp(s.meshSHA256,pair.meshSHA256) && ...
        strcmp(s.seedFrame,pair.seed) && s.elapsedSeconds==pair.elapsedS && strcmp(e.frame,s.frameSHA256));
    assert(abs(e.relativeBalance)<.005 && abs(e.storageDifferenceW)<1e-9 && abs(e.assemblyGapW)<1e-9 && ...
        max(e.window)<-8 && e.rawMax<1e-8 && e.heatChange<.002 && e.releasePeakRelativeError<1e-5);
    power=e.wallW;meshSHA=s.meshSHA256;method='SU2 BDF1 · 双支管银球的物理时间端点';
    scope='固定温界的不可压层流瞬时解，已核对实际前一时间层、储热及独立壁反力。空间误差与实验精度待验证。';
else
    assert(s.thermal.passed && abs(s.thermal.early.physicalGapFraction)<.005 && ...
        abs(s.thermal.late.physicalGapFraction)<.005 && s.thermal.thermalChange<.002 && ...
        strcmp(s.thermal.early.sourceFrame,s.frameSHA256));
    power=s.thermal.early.wallW;meshSHA=s.thermal.meshSHA256;method='SU2 8.6.0 · 双支管银球有界稳态 · 统一远场质量通量';
    scope='名义一球径双支管、湿端盖及给定壁温。稳态层流研究解；空间误差和逐面实验精度待验证。';
end
for key={'wall','stations'}
    name=[key{1} '.csv'];assert(strcmp(cm011FileSHA256(fullfile(temp,name)),s.files.(matlab.lang.makeValidName(name))));
end
w=readtable(fullfile(temp,'wall.csv'));stations=readtable(fullfile(temp,'stations.csv'));
unzip(fullfile(cases,'supported65.zip'),temp);
geom=jsondecode(fileread(fullfile(temp,'geometry.json')));
assert(isequal(geom.F,g.F) && isequal(geom.regionID,g.regionID) && max(abs(geom.V-g.V),[],'all')<1e-14, ...
    'RAH:Short:FaceOrder','当前几何原面顺序与冻结解不一致。');
n=size(g.F,1);c=s.condition;
assert(isequal(w.faceID,(1:n)') && isequal(w.regionID,g.regionID) && max(abs(w.areaM2-g.areaM2))<1e-17);
assert(max(abs(w.absolutePa-w.gaugePa-c.airPa))<1e-8 && abs(w.areaM2'*w.fluxInWm2+power)<1e-10);
assert(isequal(stations.originalFaceID,g.stations.faceID) && s.jointExperimentalAreaM2==0);
later=[];if timed && target==1,later=pair;end
if radiusD==14 && ~powerMode,stations=attachShortStationGradients(stations,s,cases,job,later);end
% 模块二：用统一的原面适配器交付压力、热流、给定壁温及对应首层流体温度。
names=["外露银球";"+Y 支管外壁";"−Y 支管外壁";"+Y 湿端盖";"−Y 湿端盖"];
condition=struct('caseLabel','Short Test 65 · 双支管银球 · 稳态','velocityMps',[-c.speedMps 0 0], ...
    'ambientPressurePa',c.airPa,'ambientTemperatureK',c.airK,'reynolds',c.re,'alphaDeg',0,'betaDeg',0);
if timed
    condition.caseLabel=sprintf('Short Test 65 · 双支管银球 · 扰动后 %.6f ms · 步长 %.6f ms',1000*s.elapsedSeconds,1000*s.dtSeconds);
    if radiusD==20
        condition.caseLabel=[condition.caseLabel ' · 20D 远场'];condition.outerRadiusD=radiusD;
    end
    condition.elapsedPhysicalTimeS=s.elapsedSeconds;condition.physicalStepSizeS=s.dtSeconds;
end
setup=struct('condition',condition,'region',names(g.regionID),'wallTemperatureK',w.wallInputK, ...
    'engineeringSolver',@empty,'localSolver',@wall);
[mesh,~,coverage]=runHybridSurfaceCase(g,struct('hybridSurface',setup));
assert(coverage.fullSurfaceComputed && abs(coverage.wholeSurfaceHeatPowerW+power)<1e-10);
mesh.hybridWall.pressureGaugePa=w.gaugePa;mesh.hybridWall.pressureGaugeMethod=mesh.hybridWall.pressureMethod;
mesh.shortSupportSurface=true;
id=mesh.evidenceIdentity;id.meshSHA256=meshSHA;
id.sourceFrame=struct('kind',method,'sha256',s.frameSHA256, ...
    'caseSHA256',s.caseSHA256,'reactionSHA256',s.releaseSHA256,'wallSHA256',s.files.wall_csv);
id.sourceTraceable=true;id.productionQualified=false;
id.scope=scope;
if timed
    % 模块三：剪切、首层速度和原面流向绑定本分支物理终帧，与工作台共用时间标签。
    id.sourceFrame.previousSHA256=e.old;id.sourceFrame.shearVTKSHA256=e.vtk;
    path=fullfile(temp,'direction.csv');assert(strcmp(cm011FileSHA256(path),s.files.direction_csv));
    vel=readtable(path);assert(isequal(vel.faceID,w.faceID) && isequal(vel.regionID,w.regionID));
    vel.Properties.VariableNames(3:8)=strcat(resolution,{'ShearPa','ShearMagnitudePa','SpeedMps','SpeedMagnitudeMps','Status','Flags'});
    mesh.shortSupportDirection=shortSupportDirectionReasons(vel,s.flagNames,{resolution});
    mesh.shortSupportDirectionGrid=resolution;mesh.shortSupportTime=pair;
end
mesh.evidenceIdentity=id;s.coverage=coverage;s.outW=power;
if powerMode
    % 模块四：外层热平衡求出球温，支管和端盖温度由用户选定的同温假设给出。
    assert(strcmp(cal.balance.earlyFrame,s.frameSHA256) && abs(s.regions(1).outW-cal.targetSphereW)<.0002*cal.targetSphereW);
    assert(max(abs(w.wallInputK-cal.temperatureK))<1e-8 && max(abs([s.regions.outW]'-cal.balance.earlyRegionW(:)))<1e-10);
    sphere=g.regionID==1;mesh.hybridSolvedWallTemperatureK=nan(n,1);
    mesh.hybridSolvedWallTemperatureK(sphere)=cal.temperatureK;mesh.hybridWallTemperatureK(sphere)=NaN;
    experiment=jsondecode(fileread(fullfile(cases,'test65.json')));
    mesh.shortSupportPower=struct('targetSphereW',cal.targetSphereW,'sphereW',s.regions(1).outW, ...
        'supportW',sum([s.regions(2:end).outW]),'temperatureK',cal.temperatureK,'measuredMeanK',(experiment.meanSurfaceF-32)*5/9+273.15);
    mesh.evidenceIdentity.scope=sprintf('Short 65 名义等温银球，球面量热功率反求壁温。两管和端盖同温并独立供热。准稳态层流，沿用实验平均表温所定膜温导热系数 %.8f W/(m·K)。',c.kFilm);
    mesh.hybridCondition.caseLabel='Short Test 65 · 球体量热功率 · 准稳态';
end
if ~isempty(later),mesh.shortSupportTemporal=readShortSupportTimeDifference(job,cases,branch,pair,s);end
if radiusD==14 && timed && (target==1 || strcmp(resolution,'timehalf'))
    az=readShortSupportAzimuth(cases,job,s,stations,temp,later);
    if ~isempty(az),mesh.shortSupportAzimuth=az;end
end
if ismember('gradientNu',stations.Properties.VariableNames),mesh.shortSupportStationTable=stations;end
    function t=empty(~,~,~)
        t=table(zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1),'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
    end
    function out=wall(~,~,~,~,ids,~)
        j=ids+1;count=numel(j);
        t=table(ids,w.absolutePa(j),w.fluxInWm2(j),repmat("SU2 8.6.0 三维不可压层流",count,1), ...
            'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
        t.pressureMethod=repmat("SU2 壁节点表压按对偶面积和五区交叠加权，绝对压加环境压力",count,1);
        t.heatFluxMethod=repmat("同态释放壁反力与节点余量通量修正，按五区交叠守恒分配，正值入壁",count,1);
        t.nearWallCellTemperatureK=w.firstK(j);
        t.nearWallTemperatureMethod=repmat("对应壁节点的首层温差标量，由对偶面积和分区交叠换算到 K",count,1);
        reason='有界稳态、同态热量守恒及热稳定均通过';
        if timed
            reason='BDF1 实际时间步、全域储热与独立释放壁行均通过';
            t.heatFluxMethod=repmat("BDF1 同刻强温界前后壁行与逐边限幅修正，独立释放复核；正值入壁",count,1);
        end
        out=struct('accepted',true,'reason',reason,'wall',t);
    end
end
function discard(folder)
% 临时目录仅展开已核实的几何包及结果包，逐项释放。
for name={'summary.json','geometry.json','parents.csv','stations.csv','supported_m.stl','wall.csv','direction.csv'}
    p=fullfile(folder,name{1});if isfile(p),delete(p);end
end
rmdir(folder);
end
