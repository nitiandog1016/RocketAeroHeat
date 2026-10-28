function report=testAxisymmetricAnnularView(folder,evidenceFile)
%TESTAXISYMMETRICANNULARVIEW 从真实已收敛NS证据核对环隙墙面、身份与曲面热账。
% 模块一：重用冻结推进结果，不重新启动长推进；原运行源码逐项匹配。
% 模块二：独立检查两片开口圆柱拓扑/法向、逐面元owner和原p/q/Tw映射。
% 模块三：单列2*pi*r*dx真面积/热量和48段弦面绘图量，不以缩放q掩盖差异。
% 模块四：篡改稳态、几何、气体或墙面须拒绝；未知温度保持NaN。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v62_validation','annular_view');end
if nargin<2,evidenceFile=fullfile(root,'outputs','v62_validation','advance','axisymmetric_advance_evidence.mat');end
assert(exist(evidenceFile,'file')==2,'AxisFV:AnnularViewTestEvidence','需要真实推进测试成功写出的完整证据。');
frozen=load(evidenceFile,'report','solutions','gas');
assert(all(isfield(frozen,{'report','solutions','gas'}))&& ...
    strcmp(frozen.report.status,'passed-real-annular-NS-time-advance-and-analytic-steady-check')&& ...
    ~frozen.report.physicalSolution&&~frozen.report.eligibleForSurfaceCoverage&& ...
    ~frozen.report.experimentalValidation&&frozen.report.newSolvedSurfaceAreaM2==0&& ...
    iscell(frozen.solutions)&&numel(frozen.solutions)>=2, ...
    'AxisFV:AnnularViewTestEvidence','不能以尝试文件或未通过推进冒充已冻结真实稳态。');
for j=1:numel(frozen.report.sourceBindings)
    item=frozen.report.sourceBindings(j);
    assert(strcmp(fileSHA256(fullfile(root,item.path)),item.sha256), ...
        'AxisFV:AnnularViewTestSource','当前源码与原推进预锁证据不匹配：%s。',item.path);
end
source=frozen.solutions{end};mesh=source.mesh;result=source.result;
assert(result.steadyConverged&&~isequal(result.Q,result.initialQ)&& ...
    result.history.residual(1)>1e4*result.residualNorm);
assert(isequaln(result.gas,axisymmetricFVGas(frozen.gas)));
sourceSHA256=fileSHA256(evidenceFile);
if ~exist(folder,'dir'),mkdir(folder);end
options=struct('evidenceSource',[char(evidenceFile) '；SHA256=' sourceSHA256]);
view=exportAxisymmetricAnnularStudy(mesh,result,fullfile(folder,'web'),options);
assert(exist(view.htmlFile,'file')==2&&exist(view.jsonFile,'file')==2&&exist(view.snapshotFile,'file')==2);
saved=load(view.snapshotFile,'mesh','result','visualMesh','data','audit','report');
assert(isequaln(saved.mesh,mesh)&&isequaln(saved.result,result),'导出不得改变源状态或几何。');
data=jsondecode(fileread(view.jsonFile));audit=saved.audit;w=result.finalInfo.wall;
F=data.faces+1;V=data.vertices;ns=48;nf=size(F,1);wallIndex=data.wallBinding.sourceWallIndex;
assert(nf==2*ns*numel(w.faceID)&&size(V,1)==2*size(mesh.X,1)*ns);
assert(all(wallIndex>=1&wallIndex<=numel(w.faceID))&& ...
    all(accumarray(wallIndex,1)==2*ns));

%% 一、曲面真面积与平三角弦面面积独立复算；面元热密度完全沿用原值
edgeA=V(F(:,2),:)-V(F(:,1),:);edgeB=V(F(:,3),:)-V(F(:,1),:);
areaVector=cross(edgeA,edgeB,2);triangleArea=vecnorm(areaVector,2,2)/2;
assert(all(triangleArea>0)&&max(abs(triangleArea-data.faceArea))<1e-16);
ri=min(mesh.R,[],'all');ro=max(mesh.R,[],'all');lengthX=max(mesh.X,[],'all')-min(mesh.X,[],'all');
exactCurvedArea=2*pi*(ri+ro)*lengthX;exactRatio=sin(pi/ns)/(pi/ns);
assert(abs(sum(triangleArea)/exactCurvedArea-exactRatio)<1e-13);
assert(abs(audit.physicalCurvedWallAreaM2-exactCurvedArea)<1e-14&& ...
    abs(audit.displayFacetedAreaM2-sum(triangleArea))<1e-14&& ...
    abs(audit.relativeAreaDeficit-(1-exactRatio))<1e-13);
assert(~audit.renderedHeatFluxWasRenormalized&&audit.endcapFaces==0);
assert(isequal(data.wallBinding.sourceFaceID,w.faceID(wallIndex))&& ...
    isequal(data.wallBinding.ownerCell,w.ownerCell(wallIndex))&& ...
    isequal(data.wallBinding.sourceFluidNormalXR,w.normal(wallIndex,:)));
assert(max(abs(accumarray(wallIndex,data.wallBinding.curvedAreaPerTriangleM2)-2*pi*w.areaPerRadian))<1e-14);
heat=field(data,'heatFlux');pressure=field(data,'pressure');Tw=field(data,'wallInput');
assert(all(heat.status==1)&&all(pressure.status==1)&&all(Tw.status==2));
assert(max(abs(heat.values*1000-w.heatFluxIntoWall(wallIndex)))<1e-11&& ...
    max(abs(pressure.values*1000-w.pressure(wallIndex)))<1e-11&& ...
    isequal(Tw.values,w.temperatureInput(wallIndex)));
trueHeat=2*pi*sum(w.heatFluxIntoWall.*w.areaPerRadian);
facetHeat=sum(triangleArea.*w.heatFluxIntoWall(wallIndex));
assert(abs(audit.physicalCurvedWallHeatW-trueHeat)<1e-12&& ...
    abs(audit.displayFacetedHeatIntegralW-facetHeat)<1e-12&& ...
    abs(facetHeat/trueHeat-exactRatio)<1e-13&&trueHeat>facetHeat&&facetHeat>0);
assert(abs(audit.physicalMovingWallPowerW-2*pi*sum(w.mechanicalPowerIntoFluidPerRadian))<1e-12);
assert(abs(audit.physicalWallEnergyFluxOutW- ...
    (audit.physicalCurvedWallHeatW-audit.physicalMovingWallPowerW))<1e-12);
assert(contains(heat.scope,'不缩放原q')&&contains(heat.methodNames{1},result.finalInfo.method));

%% 二、只画两个真实圆柱壁；内壁流体外法向-er，外壁+er，两端各有两个开口圆环
centers=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
normal=areaVector./vecnorm(areaVector,2,2);radial=[zeros(nf,1),centers(:,2:3)];
radial=radial./vecnorm(radial,2,2);expectedSign=w.normal(wallIndex,2);
assert(all(sum(normal.*radial,2).*expectedSign>.99)&&all(abs(normal(:,1))<1e-12));
radius=sqrt(V(:,2).^2+V(:,3).^2);triangleRadius=radius(F);
assert(all(max(triangleRadius,[],2)-min(triangleRadius,[],2)<1e-14));
edges=sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2);
[uniqueEdges,~,group]=unique(edges,'rows');count=accumarray(group,1);
assert(all(count<=2)&&nnz(count==1)==4*ns,'只允许四个真实开口圆环，不能有假端盖或裂缝。');
rim=uniqueEdges(count==1,:);rimX=V(rim,1);rimX=reshape(rimX,size(rim));
assert(all(abs(rimX(:,1)-rimX(:,2))<1e-14)&& ...
    all(abs(rimX(:,1)-min(mesh.X,[],'all'))<1e-14| ...
        abs(rimX(:,1)-max(mesh.X,[],'all'))<1e-14));
assert(numel(data.regions)==2&&sum([data.regions.faces])==nf);
for j=1:numel(data.regions)
    assert(isequal(sort(fieldnames(data.regions(j).fields)),sort({data.fields.key}')));
    for k=1:numel(data.fields)
        row=data.regions(j).fields.(data.fields(k).key);
        assert(abs(row.computedArea+row.prescribedArea+row.missingArea-data.regions(j).area)<1e-14);
    end
end

%% 三、数值身份、实际参数与未知量同步进入云图/探针/截面/导出合同
assert(~data.researchIdentity.physicalSolution&&~data.researchIdentity.eligibleForSurfaceCoverage&& ...
    ~data.researchIdentity.experimentalValidation&&data.researchIdentity.newProductionSolvedAreaM2==0&& ...
    ~data.researchIdentity.axialEndsArePhysicalWalls&&data.researchIdentity.steadyConverged);
assert(contains(data.caseNote,'数值校核')&&contains(data.caseNote,'新增生产面积0')&& ...
    contains(data.caseNote,'曲面')&&contains(data.caseNote,'弦面')&& ...
    data.resultCondition.gas.R==frozen.gas.R&&data.resultCondition.gas.mu==frozen.gas.mu&& ...
    data.resultCondition.residualNorm==result.residualNorm);
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    current=field(data,key{1});
    assert(all(current.status==0)&&all(isnan(current.values))&& ...
        all(~cellfun(@isempty,current.missingReason))&&current.computedAreaFraction==0);
end
assert(all(heat.methodIndex==0)&&numel(heat.methodNames)==1&& ...
    isequal(heat.scopeIndex,wallIndex-1)&&numel(heat.scopeNames)==numel(w.faceID));
assert(contains(fileread(view.htmlFile),'● 数值校核'));
assert(strcmp(fileSHA256(fullfile(root,'web','viewer.css')),fileSHA256(fullfile(folder,'web','viewer.css'))));
for name={'section.js','export.js','viewer.js'}
    original=fileread(fullfile(root,'web',name{1}));exported=fileread(fullfile(folder,'web',name{1}));
    assert(strcmp(exported,strrep(original,'工程计算','数值校核'))&&~contains(exported,'工程计算'), ...
        '只允许当前研究页将状态文字改为数值校核，不得改变既有绘图/导出算法。');
end
for j=1:numel(view.sourceBindings)
    item=view.sourceBindings(j);
    assert(strcmp(fileSHA256(fullfile(root,item.path)),item.sha256));
end

%% 四、损坏/过域快照须在写出新目录之前拒绝；不得用当前工况覆盖旧结果
refusals=0;badFolder=[tempname(folder) '_must_not_export'];
bad=result;bad.steadyConverged=false;reject(mesh,bad,'AxisFV:AnnularViewIdentity');
bad=result;bad.stopReason='maximum-steps-reached';reject(mesh,bad,'AxisFV:AnnularViewIdentity');
bad=result;bad.eligibleForSurfaceCoverage=true;reject(mesh,bad,'AxisFV:AnnularViewIdentity');
bad=result;bad.physicalSolution=true;reject(mesh,bad,'AxisFV:AnnularViewIdentity');
bad=result;bad.residualNorm=10*result.options.steadyResidualTolerance;reject(mesh,bad,'AxisFV:AnnularViewConvergence');
bad=result;bad.history.residual(end)=10*result.options.steadyResidualTolerance;reject(mesh,bad,'AxisFV:AnnularViewConvergence');
badMesh=mesh;badMesh.X(1,1)=badMesh.X(1,1)+1e-4;reject(badMesh,result,'AxisFV:AnnularViewBinding');
badMesh=mesh;badMesh.periodicI=false;reject(badMesh,result,'AxisFV:AnnularViewGeometry');
bad=result;bad.geometryBinding.volumePerRadian(1)=2*bad.geometryBinding.volumePerRadian(1);reject(mesh,bad,'AxisFV:AnnularViewBinding');
bad=result;bad.gas.R=1.01*bad.gas.R;bad.gas=rmfield(bad.gas,{'cp','cv','k'});reject(mesh,bad,'AxisFV:AnnularViewState');
bad=result;bad.gas.mu=1.01*bad.gas.mu;bad.gas=rmfield(bad.gas,{'cp','cv','k'});reject(mesh,bad,'AxisFV:AnnularViewReplay');
bad=result;bad.gas.Pr=1.01*bad.gas.Pr;bad.gas=rmfield(bad.gas,{'cp','cv','k'});reject(mesh,bad,'AxisFV:AnnularViewReplay');
bad=result;bad.finalInfo.wall.temperatureInput(1)=bad.finalInfo.wall.temperatureInput(1)+1;reject(mesh,bad,'AxisFV:AnnularViewReplay');
bad=result;bad.finalInfo.wall.heatFluxIntoWall(1)=bad.finalInfo.wall.heatFluxIntoWall(1)+1;
bad.finalInfo.wall.totalEnergyFluxOutPerRadian(1)=bad.finalInfo.wall.totalEnergyFluxOutPerRadian(1)+bad.finalInfo.wall.areaPerRadian(1);
reject(mesh,bad,'AxisFV:AnnularViewReplay');
bad=result;bad.finalInfo.wall.ownerCell(1)=mod(bad.finalInfo.wall.ownerCell(1),mesh.cellCount)+1;reject(mesh,bad,'AxisFV:AnnularViewBinding');
bad=result;bad.finalInfo.wall.centroid(1,1)=bad.finalInfo.wall.centroid(1,1)+1e-4;reject(mesh,bad,'AxisFV:AnnularViewBinding');
bad=result;bad.finalInfo.wall.areaPerRadian(1)=2*bad.finalInfo.wall.areaPerRadian(1);reject(mesh,bad,'AxisFV:AnnularViewBinding');
bad=result;bad.finalInfo.wall.pressure(1)=bad.finalInfo.wall.pressure(1)+1;reject(mesh,bad,'AxisFV:AnnularViewBinding');
bad=result;bad.finalInfo.wall.normal(1,:)=-bad.finalInfo.wall.normal(1,:);reject(mesh,bad,'AxisFV:AnnularViewBinding');
bad=result;bad.finalInfo.wall.heatFluxIntoWall(1)=NaN;reject(mesh,bad,'AxisFV:AnnularViewWall');
bad=result;bad.finalInfo.wall.temperatureInput(1)=NaN;reject(mesh,bad,'AxisFV:AnnularViewWall');
bad=result;bad.finalInfo.wall.totalEnergyFluxOutPerRadian(1)= ...
    bad.finalInfo.wall.totalEnergyFluxOutPerRadian(1)+1;reject(mesh,bad,'AxisFV:AnnularViewWall');
assert(~exist(badFolder,'dir'),'拒绝输入不得生成可误认的结果页面。');
assert(strcmp(sourceSHA256,fileSHA256(evidenceFile)),'只读冻结证据不得改变。');
report=struct('status','passed-bound-annular-surface-display-and-curved-area-audit', ...
    'evidenceFile',char(evidenceFile),'evidenceSHA256',sourceSHA256, ...
    'sourceAdvanceBindings',frozen.report.sourceBindings,'sourceSteadyResidual',result.residualNorm, ...
    'sourceAcceptedSteps',result.acceptedSteps,'view',view,'refusals',refusals, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false, ...
    'newSolvedSurfaceAreaM2',0,'fluidTemperatureIsNotEdgeTemperature',true, ...
    'wallTemperatureIsPrescribed',true,'renderingAlgorithmsUnchanged',true, ...
    'onlyExportedStudyStatusLabelsAdapted',true,'spatialResidualReplayChecked',true);
save(fullfile(folder,'annular_view_evidence.mat'),'report','sourceSHA256');
fid=fopen(fullfile(folder,'annular_view_summary.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report,'ConvertInfAndNaN',true));
fprintf('Annular view passed: %d actual wall faces -> %d triangles, curved %.12g m2 / facet %.12g m2, %d refusals; no production coverage.\n', ...
    numel(w.faceID),nf,audit.physicalCurvedWallAreaM2,audit.displayFacetedAreaM2,refusals);

    function reject(currentMesh,currentResult,identifier)
        try
            exportAxisymmetricAnnularStudy(currentMesh,currentResult,badFolder);
        catch exception
            assert(strcmp(exception.identifier,identifier),'AxisFV:AnnularViewTestRefusal', ...
                '拒绝标识错误：期望%s，实际%s。',identifier,exception.identifier);
            refusals=refusals+1;return;
        end
        error('AxisFV:AnnularViewTestNoRefusal','损坏或超域输入没有被拒绝。');
    end
end

function value=field(data,key)
value=data.fields(strcmp({data.fields.key},key));
assert(isscalar(value),'AxisFV:AnnularViewTestField','字段缺失或重复：%s。',key);
end

function value=fileSHA256(path)
% 模块五：实际文件字节哈希；用于确认原证据及原网页模块未被替换。
fid=fopen(path,'rb');assert(fid>=0,'AxisFV:AnnularViewTestSource','缺少来源文件%s。',path);
cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
