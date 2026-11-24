function [result,boundary,report]=runShortSupportImportedCase(input,cfg,progress)
% 模块一：核对用户确认的文件、面网、温界和所选时间层，然后调用试件原生入口。
boundary=shortSupportImportBoundary(cfg.shortSupportBoundary.branch);
assert(isequal(cfg.shortSupportBoundary,boundary),'RAH:Short:Boundary','当前设置与 Test 65 的冻结工况不符。');
[b,g]=bindShortSupportSTL(input,input.shortSupportImport.source,cfg.geometry);
assert(isequal(b,input.shortSupportImport),'RAH:Short:Mapping','建网后试件映射已改变。');
root=fileparts(mfilename('fullpath'));job=cfg.shortSupportJob;
switch cfg.shortSupportAction
    case 'solve'
        assert(~isempty(job) && ~isfolder(job),'RAH:Short:Job','新算需要独立的新任务目录。');
        runCM011Python({fullfile(root,'scripts','shortSupportTime65.py'),job},[job '.log'],@pulse,true);
    case 'archive'
        assert(isempty(job),'RAH:Short:Job','读取发布归档的目录设置无效。');
    otherwise
        error('RAH:Short:Action','请在计算页选择新算或读取归档。');
end
grid=['time' boundary.branch];[result,stations,s]=loadShortSupportSurface(g,job,grid);
assert(strcmp(s.meshSHA256,b.meshSHA256) && strcmp(s.seedFrame,boundary.initialFrame) && ...
    strcmp(s.caseSHA256,boundary.caseSHA256) && s.elapsedSeconds==boundary.elapsedS && s.dtSeconds==boundary.dtS);
% 模块二：面均压力、温度及方向按同一物理面取值；按导入面的实际面积折算热流密度。
order=b.referenceFaceID;area=b.areaM2;n=numel(order);inverse=zeros(n,1);inverse(order)=(1:n)';
w=result.hybridWall(order,:);facePower=w.heatFluxWm2.*b.referenceAreaM2;
w.faceID=(0:n-1)';w.heatFluxWm2=facePower./area;
w.heatFluxMethod=w.heatFluxMethod+"；按对应 STL 原面功率守恒换算";
assert(all(abs(w.heatFluxWm2.*area-facePower)<1e-12));
result.F=input.F;result.V=input.V;result.hybridWall=w;
result.hybridWallTemperatureK=result.hybridWallTemperatureK(order);
result.localGeometry.region=result.localGeometry.region(order);
direction=result.shortSupportDirection(order,:);direction.faceID=(1:n)';
result.shortSupportDirection=direction;result.shortSupportImport=b;
% 两步长的各项原始统计以参考网为测度；工作台只使用共同时刻、色标和本次文件编号的统计。
pair=result.shortSupportTime;
result.shortSupportTime=struct('elapsedS',pair.elapsedS,'flowRangeMax',pair.flowRangeMax, ...
    'heatRangeKWm2',pair.heatRangeKWm2,'gaugeRangePa',pair.gaugeRangePa);
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
reverse=find(direction.([grid 'Status'])==-1);
report.regions=regions;report.stations=stations;report.peakOutFace=face;report.peakOutWm2=peak;
report.peakReferenceFaceID=order(face);report.nearPeakFaces=near;report.reverseFaces=reverse;
report.nativePeakReferenceFaceID=s.peakOutFace;report.nativePeakFileFaceID=inverse(s.peakOutFace);
report.nativePeakWm2=s.peakOutWm2;
report.reverseAreaM2=sum(area(reverse));report.nearPeakAreaM2=sum(area(near));
report.stationRMSE=s.stationRMSE;report.equivalentRelativeError=s.equivalentRelativeError;
report.jointExperimentalAreaM2=s.jointExperimentalAreaM2;report.energy=s.energy;
result.shortSupportAction=cfg.shortSupportAction;result.shortSupportJob=job;
id=result.evidenceIdentity;id.reusedArchive=strcmp(cfg.shortSupportAction,'archive');
id.sourceFrame.stlSHA256=b.source.sha256;id.sourceFrame.stlScaleToM=b.source.scale;
id.sourceFrame.initialSHA256=s.seedFrame;id.sourceFrame.transferSHA256=b.transferSHA256;
result.evidenceIdentity=id;
assert(strcmp(cm011FileSHA256(b.source.file),b.source.sha256),'RAH:Short:Source','计算期间源 STL 已改变。');
    function pulse()
        if ~isempty(progress)
            progress(struct('status','running','message','双支管银球原生计算中'));
        end
    end
end
