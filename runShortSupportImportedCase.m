function [result,boundary,report]=runShortSupportImportedCase(input,cfg,progress)
% 模块一：核对用户确认的文件、面网、温界和所选时间层，然后调用试件原生入口。
radiusD=14;if isfield(cfg,'shortSupportRadiusD'),radiusD=cfg.shortSupportRadiusD;end
target=.25;if isfield(cfg.shortSupportBoundary,'targetDU'),target=cfg.shortSupportBoundary.targetDU;end
boundary=shortSupportImportBoundary(cfg.shortSupportBoundary.branch,radiusD,target);
powerMode=strcmp(boundary.branch,'power');
assert(isequal(cfg.shortSupportBoundary,boundary),'RAH:Short:Boundary','当前设置与 Test 65 的冻结工况不符。');
[b,g]=bindShortSupportSTL(input,input.shortSupportImport.source,cfg.geometry,radiusD);
assert(isequal(b,input.shortSupportImport),'RAH:Short:Mapping','建网后试件映射已改变。');
root=fileparts(mfilename('fullpath'));job=cfg.shortSupportJob;phaseRead=tic;lastPhase='';
switch cfg.shortSupportAction
    case 'solve'
        assert(~isempty(job) && ~isfolder(job),'RAH:Short:Job','新算需要独立的新任务目录。');
        % 全新检出可能还没有 outputs；启动 Python 前需先建立日志的父目录。
        parent=fileparts(job);if ~isfolder(parent),mkdir(parent);end
        script='shortSupportTime65.py';if radiusD==20,script='shortSupportOuter65.py';end
        if target==1,script='shortSupportLater65.py';end
        if powerMode,script='shortSupportPower65.py';end
        args={fullfile(root,'scripts',script),job};
        if target==1 && isfield(cfg,'shortSupportSourceJob') && ~isempty(cfg.shortSupportSourceJob)
            args=[args,{'--from-job',cfg.shortSupportSourceJob}];
        end
        runCM011Python(args,[job '.log'],@pulse,true);
    case 'archive'
        assert(isempty(job),'RAH:Short:Job','读取发布归档的目录设置无效。');
    case 'reopen'
        assert(radiusD==14 && isfolder(job) && strcmp(cm011FileSHA256(fullfile(job,'completed.json')),cfg.shortSupportSeal), ...
            'RAH:Short:Job','重开任务的路径或回执已改变。');
    otherwise
        error('RAH:Short:Action','请在计算页选择新算或读取归档。');
end
grid=['time' boundary.branch];if powerMode,grid='power';end
[result,stations,s]=loadShortSupportSurface(g,job,grid,radiusD,target);
if powerMode
    cal=readShortSupportPowerCase(job);seed=cal.balance.initialSHA256;
else
    seed=boundary.initialFrame;
end
if ~isempty(job)
    frame=shortSupportJobFrame(job,b,boundary);seed=frame.initial;
    assert(strcmp(frame.terminal,s.frameSHA256),'RAH:Short:Job','壁面源帧与所选任务不符。');
    if ~powerMode,assert(strcmp(frame.steady,s.baselineFrame),'RAH:Short:Job','稳态基线不符。');end
    result.evidenceIdentity.sourceFrame.completedReceiptSHA256=frame.receipt;
end
assert(strcmp(s.caseSHA256,boundary.caseSHA256));
if powerMode
    assert(strcmp(s.thermal.meshSHA256,b.meshSHA256) && strcmp(s.frameSHA256,cal.balance.earlyFrame));
else
    assert(strcmp(s.meshSHA256,b.meshSHA256) && strcmp(s.seedFrame,seed) && ...
        s.elapsedSeconds==boundary.elapsedS && s.dtSeconds==boundary.dtS);
end
% 模块二：面均压力、温度及方向按同一物理面取值；按导入面的实际面积折算热流密度。
order=b.referenceFaceID;area=b.areaM2;n=numel(order);inverse=zeros(n,1);inverse(order)=(1:n)';
if isfield(result,'shortSupportAzimuth')
    % 曲线坐标属于标准试件。每个命中位置的面号按本次导入映射，遮盖状态单独保留。
    ids=result.shortSupportAzimuth.points.referenceFaceID;found=ids>0;
    result.shortSupportAzimuth.points.fileFaceID(found)=inverse(ids(found));
    result.shortSupportAzimuth.stlSHA256=b.source.sha256;
end
w=result.hybridWall(order,:);facePower=w.heatFluxWm2.*b.referenceAreaM2;
w.faceID=(0:n-1)';w.heatFluxWm2=facePower./area;
w.heatFluxMethod=w.heatFluxMethod+"；按对应 STL 原面功率守恒换算";
assert(all(abs(w.heatFluxWm2.*area-facePower)<1e-12));
result.F=input.F;result.V=input.V;result.hybridWall=w;
result.hybridWallTemperatureK=result.hybridWallTemperatureK(order);
result.localGeometry.region=result.localGeometry.region(order);
result.shortSupportImport=b;
if powerMode
    result.hybridSolvedWallTemperatureK=result.hybridSolvedWallTemperatureK(order);
else
    direction=result.shortSupportDirection(order,:);direction.faceID=(1:n)';
    result.shortSupportDirection=direction;
    % 两步长的各项原始统计以参考网为测度；工作台只使用共同时刻、色标和本次文件编号的统计。
    pair=result.shortSupportTime;
    result.shortSupportTime=struct('elapsedS',pair.elapsedS,'flowRangeMax',pair.flowRangeMax, ...
        'heatRangeKWm2',pair.heatRangeKWm2,'gaugeRangePa',pair.gaugeRangePa);
end
report=summarizeComputedWall(area,isfinite(w.pressurePa),isfinite(w.heatFluxWm2),w.heatFluxWm2);
report.attempts=s.coverage.attempts;
for j=1:numel(report.attempts)
    report.attempts(j).targetFaceIDs=sort(inverse(report.attempts(j).targetFaceIDs+1)-1);
end
assert(report.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW+s.outW)<1e-10);
% 模块三：五区重新积分，计算文件网格上的峰值和反向区，挂接原实验测站。
regions=s.regions;
for j=1:numel(regions)
    faces=find(b.regionID==regions(j).id);q=w.heatFluxWm2(faces);
    outward=-sum(area(faces).*q);assert(abs(outward-regions(j).outW)<1e-10);
    [low,k]=min(q);regions(j).peakOutFace=faces(k);regions(j).referenceFaceID=order(faces(k));
    regions(j).minFluxIn=low;regions(j).maxFluxIn=max(q);regions(j).outW=outward;regions(j).areaM2=sum(area(faces));
end
[peak,face]=max(-w.heatFluxWm2);near=find(-w.heatFluxWm2>=.99*peak);
stations.referenceFaceID=stations.originalFaceID;stations.fileFaceID=inverse(stations.originalFaceID);
stations=removevars(stations,'originalFaceID');stations.faceMeanFluxInWm2=w.heatFluxWm2(stations.fileFaceID);
if ismember('gradientNu',stations.Properties.VariableNames),result.shortSupportStationTable=stations;end
if ~powerMode
    reverse=find(direction.([grid 'Status'])==-1);
    report.reverseFaces=reverse;report.reverseAreaM2=sum(area(reverse));
end
report.regions=regions;report.stations=stations;report.peakOutFace=face;report.peakOutWm2=peak;
report.peakReferenceFaceID=order(face);report.nearPeakFaces=near;
report.nativePeakReferenceFaceID=s.peakOutFace;report.nativePeakFileFaceID=inverse(s.peakOutFace);
report.nativePeakWm2=s.peakOutWm2;
report.nearPeakAreaM2=sum(area(near));
report.stationRMSE=s.stationRMSE;report.equivalentRelativeError=s.equivalentRelativeError;
report.jointExperimentalAreaM2=s.jointExperimentalAreaM2;
if powerMode
    report.energy=s.thermal;
    report.stations.powerFrameSHA256=repmat(string(s.frameSHA256),height(stations),1);
    report.stations.method=repmat("同态释放壁反力和逐边修正；用反求球温计算当地导热系数与 Nu",height(stations),1);
    result.shortSupportStationTable=report.stations;
else
    report.energy=s.energy;
end
if radiusD==20
    % 两个半径的面号、温界、物理时刻及步长逐一相配。差值沿用实际导入面积。
    [wall,measurements,ref]=referenceForDifference(root,boundary.branch);
    assert(ref.elapsedSeconds==s.elapsedSeconds && ref.dtSeconds==s.dtSeconds && ...
        isequal(ref.condition,s.condition) && strcmp(ref.bodySHA256,s.bodySHA256) && ...
        strcmp(ref.transferSHA256,s.transferSHA256) && isequal(wall.faceID,(1:n)') && isequal(wall.regionID(order),b.regionID) && ...
        max(abs(wall.areaM2(order)-b.referenceAreaM2))<1e-17 && isequal(wall.wallInputK(order),result.hybridWallTemperatureK));
    delta=w.heatFluxWm2-wall.fluxInWm2(order).*(b.referenceAreaM2./area);
    result.shortSupportOuter=struct('fluxChangeWm2',delta,'gaugeChangePa',w.pressureGaugePa-wall.gaugePa(order), ...
        'baseFrame',ref.frameSHA256,'baseMesh',ref.meshSHA256,'powerChangeW',s.outW-ref.energy.wallW);
    assert(abs(area'*delta+result.shortSupportOuter.powerChangeW)<1e-10);
    assert(isequal(stations.referenceFaceID,measurements.originalFaceID));
    report.stations.baseNu=measurements.computedNu;
    report.baseOutW=ref.energy.wallW;report.baseStationRMSE=ref.stationRMSE;report.baseEquivalentRelativeError=ref.equivalentRelativeError;
end
if target==1
    % 同面压力与温度按源面号取值，热流时间差按导入面的实际面积保持各面瓦数。
    d=result.shortSupportTemporal;d.gaugeChangePa=d.gaugeChangePa(order);d.firstChangeK=d.firstChangeK(order);
    areaRatio=b.referenceAreaM2./area;
    d.fluxChangeWm2=d.fluxChangeWm2(order).*areaRatio;
    d.sharedMax.timeFluxChange=d.sharedMax.timeFluxChange*max(areaRatio);
    assert(abs(area'*d.fluxChangeWm2+d.powerChangeW)<1e-10);result.shortSupportTemporal=d;
    old=pair.timeComparison.(boundary.branch);
    report.stations.earlyNu=old.earlyNu(:);report.earlySourceFrame=old.earlyFrame;
    report.temporalPowerChangeW=old.wallChangeW;report.earlyStationRMSE=old.earlyStationRMSE;
    if isfield(result,'shortSupportStationTable'),result.shortSupportStationTable=report.stations;end
end
result.shortSupportAction=cfg.shortSupportAction;result.shortSupportJob=job;
id=result.evidenceIdentity;id.reusedArchive=strcmp(cfg.shortSupportAction,'archive');
id.sourceFrame.stlSHA256=b.source.sha256;id.sourceFrame.stlScaleToM=b.source.scale;
id.sourceFrame.initialSHA256=seed;id.sourceFrame.transferSHA256=b.transferSHA256;
result.evidenceIdentity=id;
if ~powerMode,boundary.initialFrame=seed;end
assert(strcmp(cm011FileSHA256(b.source.file),b.source.sha256),'RAH:Short:Source','计算期间源 STL 已改变。');
    function pulse()
        if ~isempty(progress)
            if toc(phaseRead)>1
                phaseRead=tic;path=fullfile(job,'progress.json');
                if isfile(path)
                    entry=jsondecode(fileread(path));
                    if ~strcmp(entry.message,lastPhase),lastPhase=entry.message;progress(struct('status','started','message',entry.message));end
                end
            end
            progress(struct('status','running','message','双支管银球原生计算中'));
        end
    end
end

function [wall,stations,s]=referenceForDifference(root,branch)
% 差分只读已绑定的基线两表，避免再次构建 32,480 面工作台快照及其方法列。
cases=fullfile(root,'cases','short_sphere');p=readShortSupportTimePair(cases);
file=fullfile(cases,['supportedTime65' upper(branch(1)) branch(2:end) '.zip']);
assert(strcmp(cm011FileSHA256(file),p.([branch 'Archive'])));
folder=tempname;mkdir(folder);cleanup=onCleanup(@()discardReference(folder)); %#ok<NASGU>
unzip(file,folder);s=jsondecode(fileread(fullfile(folder,'summary.json')));
assert(strcmp(s.frameSHA256,p.([branch 'Frame'])) && strcmp(s.meshSHA256,p.meshSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(folder,'wall.csv')),s.files.wall_csv) && ...
    strcmp(cm011FileSHA256(fullfile(folder,'stations.csv')),s.files.stations_csv));
wall=readtable(fullfile(folder,'wall.csv'));stations=readtable(fullfile(folder,'stations.csv'));
end
function discardReference(folder)
for name={'summary.json','wall.csv','stations.csv','direction.csv'}
    file=fullfile(folder,name{1});if isfile(file),delete(file);end
end
rmdir(folder);
end
