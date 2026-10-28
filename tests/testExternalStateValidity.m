function report=testExternalStateValidity(outputFolder,sourcePath)
%TESTEXTERNALSTATEVALIDITY 实际20度箭体的外缘状态真实性及热流整链路校核。
% 模块一：读取v55真实19752面快照；不重算STL、网格、外流或准流线。
% 模块二：以当前连续介质筛选/热流入口重新计算，核对域外占位不能
% 充当真实外缘或稀薄流。上游热力未知仍传播，不能为保覆盖率放宽门限。
% 模块三：独立图可达性与动理论公式校核；检查有效结果/输入不变并
% 保存当前实际工况MAT、JSON与离线Web例子。不是新的区域物理/实验验证。
% sourcePath显式指定时只接受该文件；省略时依次找本地v55输出和完整包
% 的继承v55快照。不会以fresh v56结果代替v55对照，不回退物理场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1 || isempty(outputFolder),outputFolder=fullfile(root,'outputs/v56_validation/state_validity');end
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end
if nargin<2,sourcePath=[];end
[source,sourceResolution]=resolveBaseline(root,sourcePath);
a=load(source);old=a.result;cfg=a.cfg;fs=old.boundaryLayerEdge.freestream;
nf=size(old.F,1);area=old.area(:);totalArea=sum(area);
assert(nf==19752 && a.meshFaces==nf && cfg.flight.alpha==20 && fs.alpha==20);
assert(cfg.flight.Mach==5 && fs.Mach==5 && fs.altitude==30000 && fs.beta==0);
assert(abs(old.heatingCoverage.summary.computedAreaFraction-.7746161231243949)<1e-12, ...
    'RocketAeroHeat:StateValidity:Baseline','必须使用v55的77.4616%%对照，不接受fresh v56结果作为旧基线。');
e=old.boundaryLayerEdge;path=old.surfaceStreamlines;
assert(isequaln(old.externalFlow.edge,e));
placeholder=logical(e.thermalExpansionOutOfRangeFace(:));
oldDirect=logical(old.continuumRegime.directUnresolvedFace(:));
oldDown=logical(old.continuumRegime.downstreamUnresolvedFace(:));
assert(nnz(oldDirect)==164 && nnz(oldDown)==3053);
assert(all(placeholder(oldDirect)) && nnz(placeholder)==336);
assert(all(e.temperature(oldDirect)<80));
assert(~old.axisymmetricMarch.applied,'本工况不应混入零攻角轴向推进');

%% 1. 当前入口重算：外流和路径全部保留，只有热流/有效性/覆盖按新合同评估
timer=tic;result=assessContinuumRegime(old,fs,cfg);
screen=result.continuumRegime;
result=calcSurfaceHeating(result,fs,cfg);elapsed=toc(timer);
c=result.continuumRegime;h=result.heating;coverage=result.heatingCoverage;
assert(isequaln(screen,c),'独立筛选与热流内部筛选必须一致');
assert(isequal(result.V,old.V) && isequal(result.F,old.F));
assert(isequaln(result.boundaryLayerEdge,old.boundaryLayerEdge));
assert(isequaln(result.surfaceStreamlines,old.surfaceStreamlines));
assert(isequaln(result.externalFlow,old.externalFlow));
assert(isequaln(result.pressure,old.pressure));
assert(isequal(h.wallTemperatureByFace,old.heating.wallTemperatureByFace));
assert(cfg.flight.Twall==old.heating.wallTemperature && h.wallTemperature==cfg.flight.Twall);
assert(isequal(c.directStateUnresolvedFace,placeholder));
assert(all(isnan(c.knudsen(placeholder))) && all(isnan(c.edgeKnudsen(placeholder))));
assert(all(isnan(c.edgeMeanFreePath(placeholder))) && all(c.regimeID(placeholder)==0));
assert(~any(c.knudsenUnresolvedFace & placeholder));
assert(all(c.stateUnresolvedFace(oldDirect|oldDown)));
assert(~any(c.knudsenUnresolvedFace(oldDirect|oldDown)), ...
    '旧164+3053占位路径不能继续冒充真实Kn超限');
assert(all(c.applicabilityID(placeholder)==6));
assert(all(c.applicabilityID(c.downstreamStateUnresolvedFace)==7));
assert(all(isnan(h.heatFlux(c.stateUnresolvedFace))));
assert(all(~h.eligibleFace(c.stateUnresolvedFace)));
assert(~any(coverage.flags.rarefaction & c.stateUnresolvedFace));
assert(all(coverage.flags.externalFlow(placeholder)));
assert(all(coverage.flags.upstreamExternalState(c.downstreamStateUnresolvedFace)));

%% 2. 独立图可达性核对两种风险传播；局部真实Kn仍按原物理门限筛选
physical=logical(e.algebraicallyValid(:)) & ~placeholder;
kinetic=e.dynamicViscosity(:)./e.pressure(:).*sqrt(pi*fs.R*e.temperature(:)/2);
independentKn=max(kinetic/c.characteristicLength,c.freestreamKnudsen);
independentKn(~physical)=NaN;
assert(max(abs(independentKn(physical)-c.knudsen(physical)))<1e-12);
realDirect=physical & isfinite(independentKn) & independentKn>=c.maxContinuumKnudsen;
assert(isequal(realDirect,c.directKnudsenUnresolvedFace));
stateReach=reachable(path,placeholder,nf);realKnReach=reachable(path,realDirect,nf);
assert(isequal(stateReach,c.stateUnresolvedFace));
assert(isequal(realKnReach,c.knudsenUnresolvedFace));
assert(isequal(stateReach|realKnReach,c.unresolvedFace));
assert(all(isnan(h.heatFlux(realKnReach))),'物理Kn超限仍须阻断热流');
assert(all(c.knudsen(physical & c.downstreamStateUnresolvedFace)<.01), ...
    '本次实际下游当地状态可在连续域，但上游未知不应放行热流');

%% 3. 未受影响的有效结果完整保留；新增未知传播损失必须如实记录
oldValid=logical(old.heating.eligibleFace(:));newValid=logical(h.eligibleFace(:));
retained=oldValid & newValid;removed=oldValid & ~newValid;added=~oldValid & newValid;
assert(~any(added),'状态真实性修复不能制造新的区域物理解');
assert(nnz(retained)>10000);
assert(isequal(h.faceMethodID(retained),old.heating.faceMethodID(retained)));
assert(isequaln(h.heatFlux(retained),old.heating.heatFlux(retained)));
assert(isequaln(h.recoveryTemperature(retained),old.heating.recoveryTemperature(retained)));
assert(isequaln(h.heatTransferCoefficient(retained),old.heating.heatTransferCoefficient(retained)));
assert(all(c.stateUnresolvedFace(removed)),'新增未解只能来自本轮真实上游状态阻断');
assert(isequal(isfinite(h.heatFlux),newValid));
assert(coverage.summary.computedAreaFraction<=old.heatingCoverage.summary.computedAreaFraction);

%% 4. 查看器逐场量合同：占位的外缘压力/静温不显示为已计算；给定量不变
viewer=buildAerothermalViewerData(result);oldViewer=buildAerothermalViewerData(old);
edgePressure=field(viewer,'edgePressure');edgeTemperature=field(viewer,'edgeTemperature');
assert(all(edgePressure.status(placeholder)==0) && all(isnan(edgePressure.values(placeholder))));
assert(all(edgeTemperature.status(placeholder)==0) && all(isnan(edgeTemperature.values(placeholder))));
assert(all(edgeTemperature.values(edgeTemperature.status==1)>=80));
assert(~any(edgeTemperature.status(oldDirect)>0));
pressure=field(viewer,'pressure');oldPressure=field(oldViewer,'pressure');
assert(isequaln(pressure.values,oldPressure.values) && isequal(pressure.status,oldPressure.status));
wall=field(viewer,'wallInput');assert(all(wall.status==2));
assert(isequaln(wall.values,old.heating.wallTemperatureByFace));
solved=field(viewer,'wallSolved');assert(all(solved.status==0) && all(isnan(solved.values)));
heat=field(viewer,'heatFlux');assert(isequal(heat.status,double(newValid)));
assert(all(isnan(heat.values(c.stateUnresolvedFace))));
assert(abs(heat.computedAreaFraction-sum(area(newValid))/totalArea)<1e-13);
assert(all(isfinite(edgePressure.values(physical))) && all(edgePressure.status(physical)==1));
assert(all(isfinite(edgeTemperature.values(physical))) && all(edgeTemperature.status(physical)==1));
assert(isequaln(edgeTemperature.values(physical),e.temperature(physical)));
assert(isequaln(edgePressure.values(physical),e.pressure(physical)/1000));

%% 5. 保存当前真实工况结果和精确面积，不把旧诊断例子冒充当前重算结果
report=struct('date','2026-10-14','sourceSnapshot',source,'sourceResolution',sourceResolution, ...
    'baselineVersion','v55','tested',true, ...
    'experimentalValidation',false,'fullAreaSolved',false, ...
    'meshRecomputed',false,'externalFlowRecomputed',false,'streamlinesRecomputed',false, ...
    'heatingAndValidityRecomputed',true,'seconds',elapsed,'flight',cfg.flight, ...
    'faces',nf,'totalAreaM2',totalArea,'oldComputed',stats(oldValid,area), ...
    'newComputed',stats(newValid,area),'retainedExactly',stats(retained,area), ...
    'removedByUnknownUpstream',stats(removed,area),'newSolvedHeatArea',stats(added,area), ...
    'oldDirectKn',stats(oldDirect,area),'oldDownstreamKn',stats(oldDown,area), ...
    'directUnknownThermal',stats(c.directStateUnresolvedFace,area), ...
    'downstreamUnknownThermal',stats(c.downstreamStateUnresolvedFace,area), ...
    'totalUnknownThermal',stats(c.stateUnresolvedFace,area), ...
    'directPhysicalKn',stats(c.directKnudsenUnresolvedFace,area), ...
    'downstreamPhysicalKn',stats(c.downstreamKnudsenUnresolvedFace,area), ...
    'validitySummary',c.summary,'reasonNames',{coverage.reasonNames}, ...
    'countByReason',coverage.countByReason,'areaByReasonM2',coverage.areaByReason, ...
    'scope',['只从同一20度实际快照重算热流/有效性，源外流/准流线没有重新求解；' ...
        '占位不算真实外缘或稀薄流，输入壁温不算求得壁温；' ...
        '31度收缩、横流、再压缩与分离区域物理仍未闭合，无新增已求解面积']);
webFolder=fullfile(outputFolder,'actual_rocket');exportAerothermalWeb(result,webFolder);
writeText(fullfile(outputFolder,'state_validity_actual_rocket.json'),jsonencode(report));
writeText(fullfile(webFolder,'results.json'),jsonencode(viewer,'ConvertInfAndNaN',true));
save(fullfile(outputFolder,'state_validity_actual_rocket.mat'),'result','report','cfg','fs', ...
    'placeholder','oldDirect','oldDown','retained','removed','added');
faceID=(1:nf)';oldHeat=old.heating.heatFlux;heatFlux=h.heatFlux;faceArea=area;
oldValid=logical(oldValid);newValid=logical(newValid);reasonCode=coverage.reasonCode;
physicalKn=c.knudsen;diagnosticKn=c.diagnosticKnudsen;unknownDirect=c.directStateUnresolvedFace;
unknownDownstream=c.downstreamStateUnresolvedFace;physicalKnBlocked=c.knudsenUnresolvedFace;
edgeTemperatureDisplayed=edgeTemperature.values;edgeTemperatureStatus=edgeTemperature.status;
edgePressureDisplayedKpa=edgePressure.values;edgePressureStatus=edgePressure.status;
writetable(table(faceID,faceArea,oldHeat,heatFlux,oldValid,newValid,removed,reasonCode,physicalKn, ...
    diagnosticKn,unknownDirect,unknownDownstream,physicalKnBlocked,edgeTemperatureDisplayed, ...
    edgeTemperatureStatus,edgePressureDisplayedKpa,edgePressureStatus), ...
    fullfile(outputFolder,'state_validity_faces.csv'));
fprintf('Actual rocket state validity passed: %d faces, alpha20; physical unknown %d direct + %d downstream; old heat %.9f%% -> new %.9f%%; removed %d faces / %.12g m2; retained %d exact; no new region solve.\n', ...
    nf,nnz(c.directStateUnresolvedFace),nnz(c.downstreamStateUnresolvedFace), ...
    100*report.oldComputed.areaFraction,100*report.newComputed.areaFraction,nnz(removed), ...
    report.removedByUnknownUpstream.areaM2,nnz(retained));
end
function [source,resolution]=resolveBaseline(root,explicit)
% 模块作用：解析同一v55证据的不同磁盘布局；不查找或回退计算数据。
if ~isempty(explicit)
    if ~(ischar(explicit) || (isstring(explicit) && isscalar(explicit)))
        error('RocketAeroHeat:StateValidity:SourcePath','显式sourcePath须为单个文件路径。');
    end
    source=char(explicit);
    if ~isfile(source)
        error('RocketAeroHeat:StateValidity:MissingExplicitSource', ...
            '显式指定的v55快照不存在：%s。不会改用其他候选文件。',source);
    end
    resolution='explicit v55 baseline file';return;
end
candidates={fullfile(root,'outputs/v55_validation/workflow_workbench.mat'), ...
    fullfile(root,'validation/inherited_v55/v55/workflow_workbench.mat')};
labels={'local v55 output snapshot','packaged inherited v55 validation snapshot'};
for k=1:numel(candidates)
    if isfile(candidates{k}),source=candidates{k};resolution=labels{k};return;end
end
error('RocketAeroHeat:StateValidity:MissingBaseline', ...
    ['找不到v55实际20度对照快照。已检查：\n%s\n%s\n' ...
    '请使用第二参数明确指定v55的workflow_workbench.mat；不会读取fresh v56基线。'], ...
    candidates{1},candidates{2});
end
function mask=reachable(path,sources,n)
% 模块作用：与生产顺序循环独立的有向图可达性检查。
parent=double(path.predecessorFace(:));child=find(parent>0);roots=find(sources);
G=digraph([parent(child);repmat(n+1,numel(roots),1)],[child;roots],[],n+1);
distance=distances(G,n+1);mask=isfinite(distance(1:n))';
end
function f=field(data,key)
id=find(strcmp({data.fields.key},key));assert(isscalar(id));f=data.fields(id);
end
function r=stats(mask,area)
r=struct('faces',nnz(mask),'areaM2',sum(area(mask)),'areaFraction',sum(area(mask))/sum(area));
end
function writeText(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',value);
end
