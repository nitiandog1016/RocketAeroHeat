function [result,report]=runThermistorImportedCase(input,cfg,progress)
% 模块一：先核对冻结文件、刚体姿态、实际面表及体网身份，再执行用户选择的计算。
p=input.thermistorPrepared;binding=input.thermistorImport;
[standard,source]=thermistorWorkbenchModel('leads');
assert(strcmp(source.bodySHA256,p.bodySHA256) && cfg.thermistorAngularCells==96 && strcmp(cfg.thermistorVariant,'sidewind'));
[ok,file]=fileattrib(cfg.geometry.stlFile);
assert(ok && strcmpi(file.Name,binding.source.file),'RAH:Bead:Source','冻结工况和准备时的 STL 路径不同。');
check=bindThermistorSTL(input,standard,binding.source,cfg.geometry);
assert(isequal(check,binding),'RAH:Bead:Mapping','导入面与建网时的对应已改变。');
standard.thermistorPrepared=p;
if isfield(cfg.thermistorBoundary,'leadThermalModel')
    [result,report]=runThermistorWireWorkbench(standard,cfg,progress);
else
    [result,report]=runThermistorMeasuredPoint(standard,cfg,progress);
end
assert(strcmp(cm011FileSHA256(binding.source.file),binding.source.sha256),'RAH:Bead:Source','计算期间 STL 源文件已改变。');
% 模块二：压力和温度取同一物理面；换算热流密度，使每个实际 STL 面的功率守恒。
order=binding.referenceFaceID;n=numel(order);inverse=zeros(n,1);inverse(order)=1:n;
area=binding.areaM2;ratio=binding.referenceAreaM2./area;
w=result.hybridWall(order,:);w.faceID=(0:n-1)';
power=w.heatFluxWm2.*binding.referenceAreaM2;
w.heatFluxWm2=w.heatFluxWm2.*ratio;
w.heatFluxMethod=w.heatFluxMethod+"；按对应三角面的面积比守恒换算至 STL";
assert(max(abs(area.*w.heatFluxWm2-power))<1e-12);
result.F=input.F;result.V=input.V;result.hybridWall=w;
result.hybridWallTemperatureK=result.hybridWallTemperatureK(order);
result.localGeometry.region=result.localGeometry.region(order);
for field={'thermistorWireFluxChange','thermometryFluxChange'}
    k=field{1};
    if isfield(result,k),q=result.(k);result.(k)=q(order).*ratio;end
end
% 模块三：用导入文件的面积、编号、区域及符号重建报告。历史对照的数值保留原测度。
coverage=summarizeComputedWall(area,isfinite(w.pressurePa),isfinite(w.heatFluxWm2),w.heatFluxWm2);
assert(coverage.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW-coverage.wholeSurfaceHeatPowerW)<1e-12);
for field=fieldnames(coverage)',k=field{1};report.(k)=coverage.(k);end
assert(report.engineeringFirstAreaM2==0);
report.regionalAddedAreaM2=coverage.jointComputedAreaM2;
for j=1:numel(report.attempts)
    report.attempts(j).targetFaceIDs=sort(inverse(report.attempts(j).targetFaceIDs+1)-1);
end
regions=report.regions;
for j=1:numel(regions)
    take=strcmp(result.localGeometry.region,regions(j).name);
    assert(nnz(take)==regions(j).faces && abs(sum(power(take))+regions(j).outwardW)<1e-12);
    regions(j).areaM2=sum(area(take));
end
[~,peak]=min(w.heatFluxWm2);
report.regions=regions;result.thermistorSidewind.regions=regions;result.thermistorSidewind.peakFaceID=peak;
if isfield(result,'thermistorWire')
    m=result.thermistorWire;m.peakFaceID=peak;m.regions=regions;m.areaM2=sum(area);
    result.thermistorWire=m;report.wire=m;
else
    m=result.thermistorMeasuredPoint;m.peakFaceID=peak;m.regions=regions;m.areaM2=sum(area);
    if isfield(m,'coldComparison')
        % 比较误差和峰值在标准网度量上取得，另存其对应的导入面号。
        c=m.coldComparison;c.faceNumbering='标准网（1 起点）；importedFaceID 为同一位置的 STL 面号';
        c.importedFaceID=inverse(c.peakFaceID);
        for field=fieldnames(c.references)'
            k=field{1};ref=c.references.(k);ref.importedPeakFaceID=inverse(ref.referencePeakFaceID);
            for measure={'pressurePa','fluxWm2','commonFirstK','relativeFlux'}
                key=measure{1};ref.(key).importedFaceID=inverse(ref.(key).faceID);
            end
            c.references.(k)=ref;
        end
        m.coldComparison=c;
    end
    result.thermistorMeasuredPoint=m;report.measuredPoint=m;
end
result.thermistorImport=binding;
result.evidenceIdentity.sourceFrame.stlSHA256=binding.source.sha256;
result.evidenceIdentity.sourceFrame.stlScaleToM=binding.source.scale;
end
