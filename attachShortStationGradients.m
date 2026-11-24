function stations=attachShortStationGradients(stations,s,cases,job)
% 模块一：本地新任务优先使用随任务封存的七站剖面；旧任务和发布归档按实际终帧匹配。
local=false;
if ~isempty(job)
    receipt=jsondecode(fileread(fullfile(job,'completed.json')));
    local=isfield(receipt,'stationGradientsSHA256');
end
if local
    file=fullfile(job,'export','stationGradients.json');expected=receipt.stationGradientsSHA256;
else
    file=fullfile(cases,'stationGradients65.json');
    expected='8245b78b562d770419cbc1d5efaeacde94e9b612883a0ec9ab04e54eb953bd9d';
end
assert(strcmp(cm011FileSHA256(file),expected),'RAH:Short:Gradient','七站梯度文件与来源记录不一致。');
data=jsondecode(fileread(file));
assert(strcmp(data.schema,'short65-station-gradient-v1'),'RAH:Short:Gradient','七站梯度文件格式有误。');
if isfield(s,'meshSHA256'),mesh=s.meshSHA256;else,mesh=s.thermal.meshSHA256;end
matches=strcmp({data.sets.frameSHA256},s.frameSHA256) & strcmp({data.sets.meshSHA256},mesh);
if ~local && ~any(matches),return;end
assert(nnz(matches)==1,'RAH:Short:Gradient','本任务没有唯一匹配的七站体场。');g=data.sets(matches);
assert(strcmp(g.caseSHA256,s.caseSHA256) && strcmp(g.stationsSHA256,s.files.stations_csv), ...
    'RAH:Short:Gradient','七站位置或工况与主面表不符。');
if isfield(s,'energy'),assert(g.elapsedS==s.elapsedSeconds);else,assert(isempty(g.elapsedS) || isnan(g.elapsedS));end
% 模块二：保留测站面号和壁反力主值，只在具备法向剖面的站点填入独立对照值。
count=height(stations);assert(count==7 && numel(g.rows)==count);
fields={'gradientNu','gradientInWm2','linearNu','wideNu','linearInWm2','wideInWm2','firstDistanceM','layerCount','sectionNormalAngleDeg','joinPlaneClearanceM'};
for k=1:numel(fields),stations.(fields{k})=nan(count,1);end
stations.gradientStatus=strings(count,1);stations.wireGradientMagnitudeKPerM=nan(count,1);
stations.fluidGradientKPerM=nan(count,1);stations.secondDistanceM=nan(count,1);
stations.firstSampleK=nan(count,1);stations.secondSampleK=nan(count,1);
for j=1:count
    if iscell(g.rows),row=g.rows{j};else,row=g.rows(j);end
    assert(row.faceID==stations.originalFaceID(j) && row.angleDeg==stations.angleDeg(j) && row.experimentNu==stations.reportedNu(j));
    assert(abs(row.reactionInWm2-stations.fluxInWm2(j))<1e-9 && abs(row.reactionNu-stations.computedNu(j))<1e-10);
    stations.wireGradientMagnitudeKPerM(j)=row.wireGradientMagnitudeKPerM;
    if strcmp(row.status,'available')
        for k=1:numel(fields),stations.(fields{k})(j)=row.(fields{k});end
        stations.fluidGradientKPerM(j)=row.gradientInWm2/g.kFilm;
        stations.secondDistanceM(j)=row.samples(3).distanceM;
        stations.firstSampleK(j)=row.samples(2).temperatureK;stations.secondSampleK(j)=row.samples(3).temperatureK;
        stations.gradientStatus(j)="已重建";
    else
        assert(strcmp(row.status,'unavailable'));stations.gradientStatus(j)=string(row.reason);
    end
end
stations.gradientKFilmWmK=repmat(g.kFilm,count,1);
stations.gradientMethod=repmat("名义球面外法线，实际弦插值，壁点与前两层非等距二次导数；膜温导热系数，正热流入壁",count,1);
stations.gradientFrameSHA256=repmat(string(s.frameSHA256),count,1);
end
