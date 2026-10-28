function report=testTmrAxisymmetricSteadyView(folder,evidenceFile)
%TESTTMRAXISYMMETRICSTEADYVIEW 真实TMR稳态墙面/三维数据/面积热账定向核对。
% 模块一：只读已收敛continue1检查点，通过当前完整空间算子逐值重放导出。
% 模块二：独立复算旋转真面积、锥面参数三角面积和显示三角面积/热积分。
% 模块三：核对未知温度、逐面算法/状态/源面映射及伪时间身份。
% 模块四：篡改气体、边界、墙热、几何或稳态合同均须在生成页面前拒绝。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','tmr_steady_view');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v63_validation','tmr_hlle_continue1','tmr_implicit_budget.mat');
    if ~exist(evidenceFile,'file'),evidenceFile=fullfile(root,'validation','v63','tmr_hlle_continue1','tmr_implicit_budget.mat');end
end
sourceSHA=fileSHA(evidenceFile);source=load(evidenceFile,'mesh','implicit','boundary','report');
mesh=source.mesh;result=source.implicit;boundary=source.boundary;
assert(result.steadyConverged&&result.residualNorm<=2e-6&&source.report.steadyConverged);
sourceFiles={'web/index.html','web/viewer.css','web/viewer.js','web/section.js','web/export.js'};
sourceUIHashes=cellfun(@(p)fileSHA(fullfile(root,p)),sourceFiles,'UniformOutput',false);
opts=struct('azimuthSegments',48,'evidenceSource',evidenceFile);
view=exportTmrAxisymmetricStudy(mesh,result,boundary,fullfile(folder,'web'),opts);
data=jsondecode(fileread(view.jsonFile));w=result.finalInfo.wall;map=data.wallBinding.sourceWallIndex;
V=data.vertices;F=data.faces+1;ns=opts.azimuthSegments;nf=size(F,1);
assert(nf==3840&&size(V,1)==1968&&numel(data.faceArea)==nf&&numel(unique(map))==40);
q=getField(data,'heatFlux');p=getField(data,'pressure');tw=getField(data,'wallInput');
assert(isequal(q.values*1000,w.heatFluxIntoWall(map))||max(abs(q.values*1000-w.heatFluxIntoWall(map)))<1e-10);
assert(max(abs(p.values*1000-w.pressure(map)))<1e-10&&isequal(tw.values,w.temperatureInput(map)));
assert(all(q.status==1)&&all(p.status==1)&&all(tw.status==2));
assert(isequal(data.wallBinding.sourceFaceID,w.faceID(map))&&isequal(data.wallBinding.ownerCell,w.ownerCell(map)));
assert(all(q.methodIndex==0)&&contains(q.methodNames{1},result.finalInfo.method)&& ...
    isequal(q.scopeIndex,map-1)&&numel(q.scopeNames)==40&&contains(q.scope,'不缩放signed q'));

%% 一、独立几何与热积分：参数三角形真面积在锥面上不相等
crossArea=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
facetArea=vecnorm(crossArea,2,2)/2;assert(max(abs(facetArea-data.faceArea))<1e-15);
x=mesh.X(:,1);r=mesh.R(:,1);dx=diff(x);dr=diff(r);L=hypot(dx,dr);
segment=data.wallBinding.meridionalSegment;parameter=data.wallBinding.parameterTriangle;
expectedCurved=zeros(nf,1);expectedFacet=zeros(nf,1);delta=2*pi/ns;
for j=1:nf
    ix=segment(j);
    if parameter(j)==1
        expectedCurved(j)=L(ix)*delta*(r(ix)+2*r(ix+1))/6;
        expectedFacet(j)=r(ix+1)*sin(delta/2)*sqrt(dx(ix)^2+dr(ix)^2*cos(delta/2)^2);
    else
        expectedCurved(j)=L(ix)*delta*(2*r(ix)+r(ix+1))/6;
        expectedFacet(j)=r(ix)*sin(delta/2)*sqrt(dx(ix)^2+dr(ix)^2*cos(delta/2)^2);
    end
end
assert(max(abs(expectedCurved-data.wallBinding.curvedAreaPerTriangleM2))<1e-15&& ...
    max(abs(expectedFacet-facetArea))<1e-15);
sourceCurved=accumarray(map,expectedCurved,[40 1]);assert(max(abs(sourceCurved-2*pi*w.areaPerRadian))<1e-14);
cone=abs(dr(segment))>1e-10;cylinder=~cone;
assert(any(cone)&&any(cylinder));
assert(any(abs(expectedCurved(cone)-2*pi*w.areaPerRadian(map(cone))/(2*ns))>1e-8), ...
    'AxisFV:TmrViewTestArea','锥面参数对角真面积被错误套用等分圆柱公式。');
audit=data.surfaceIntegralAudit;
trueArea=sum(2*pi*w.areaPerRadian);trueHeat=sum(2*pi*w.areaPerRadian.*w.heatFluxIntoWall);
facetHeat=sum(facetArea.*q.values*1000);
assert(abs(audit.physicalCurvedWallAreaM2-trueArea)<1e-13&& ...
    abs(audit.displayFacetedAreaM2-sum(facetArea))<1e-13&& ...
    abs(audit.physicalCurvedWallHeatW-trueHeat)<1e-9&& ...
    abs(audit.displayFacetedHeatIntegralW-facetHeat)<1e-9&& ...
    abs(sum(expectedCurved.*q.values*1000)-trueHeat)<1e-9&&~audit.renderedHeatFluxWasRenormalized);
centers=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;radius=hypot(centers(:,2),centers(:,3));
normals=crossArea./vecnorm(crossArea,2,2);fluid=[w.normal(map,1),w.normal(map,2).*centers(:,2:3)./radius];
assert(all(sum(normals.*fluid,2)>.995),'AxisFV:TmrViewTestNormal','绘图面法向不对应源流体外法向。');
edges=sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2);[uniqueEdges,~,group]=unique(edges,'rows');counts=accumarray(group,1);
assert(all(counts<=2)&&nnz(counts==1)==2*ns&&audit.endcapFaces==0);
rim=uniqueEdges(counts==1,:);rimX=reshape(V(rim,1),size(rim));
assert(all(abs(rimX(:,1)-rimX(:,2))<1e-14)&&all(abs(rimX(:,1)-x(1))<1e-14|abs(rimX(:,1)-x(end))<1e-14));

%% 二、字段身份与未知温度保持原义；共享工作台源文件字节保持不变
assert(~data.researchIdentity.physicalSolution&&~data.researchIdentity.eligibleForSurfaceCoverage&& ...
    ~data.researchIdentity.experimentalValidation&&~data.researchIdentity.wholeRocketValidated&& ...
    data.researchIdentity.newProductionSolvedAreaM2==0&&data.researchIdentity.currentSpatialResidualReplayPassed);
assert(strcmp(data.resultCondition.timeKind,'pseudo-time-steady-search')&& ...
    data.resultCondition.pseudoTimeS==result.pseudoTime&&data.resultCondition.physicalTimeAdvancedS==0&& ...
    ~isfield(data.resultCondition,'timeS')&&~data.researchIdentity.pseudoTimeIsPhysical);
assert(contains(data.caseNote,'非整箭')&&contains(data.caseNote,'非物理时间')&& ...
    contains(data.caseNote,'新增生产面积0')&&contains(fileread(view.htmlFile),'● 数值研究'));
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    f=getField(data,key{1});assert(all(isnan(f.values))&&all(f.status==0)&&all(~cellfun(@isempty,f.missingReason)));
end
for j=1:numel(data.regions)
    for k=1:numel(data.fields)
        f=data.regions(j).fields.(data.fields(k).key);
        assert(abs(f.computedArea+f.prescribedArea+f.missingArea-data.regions(j).area)<1e-13);
    end
end
for j=1:numel(sourceFiles),assert(strcmp(sourceUIHashes{j},fileSHA(fullfile(root,sourceFiles{j}))));end
assert(strcmp(fileSHA(fullfile(root,'web','viewer.css')),fileSHA(fullfile(folder,'web','viewer.css'))));
snapshot=load(view.snapshotFile,'mesh','result','boundary');
assert(isequaln(snapshot.mesh,mesh)&&isequaln(snapshot.result,result)&&isequaln(snapshot.boundary,boundary));

%% 三、拒绝伪收敛、换气体/边界、假壁热或几何错配；不得生成误导页面
refusals=struct('label',{},'identifier',{});badFolder=[tempname(folder) '_must_not_export'];
bad=result;bad.steadyConverged=false;reject('not-steady',mesh,bad,boundary);
bad=result;bad.physicalSolution=true;reject('promoted-physics',mesh,bad,boundary);
bad=result;bad.eligibleForSurfaceCoverage=true;reject('promoted-area',mesh,bad,boundary);
bad=result;bad.stopReason='maximum-steps-reached';reject('budget-terminal',mesh,bad,boundary);
bad=result;bad.options.steadyResidualTolerance=2e-3;reject('loosened-threshold',mesh,bad,boundary);
bad=result;bad.options.steadyRequiredSteps=1;reject('missing-consecutive-steps',mesh,bad,boundary);
bad=result;bad.history.residual(end-1)=3e-6;reject('invalid-history',mesh,bad,boundary);
bad=result;bad.time=bad.pseudoTime;reject('pseudo-time-as-time',mesh,bad,boundary);
bad=result;bad.pseudoTimeIsPhysical=true;reject('pseudo-time-physical',mesh,bad,boundary);
badMesh=mesh;badMesh.X(1,1)=badMesh.X(1,1)+1e-4;reject('nodal-mismatch',badMesh,result,boundary);
badMesh=mesh;badMesh.gradientX(1,1)=badMesh.gradientX(1,1)+.01;reject('gradient-mismatch',badMesh,result,boundary);
bad=result;bad.geometryBinding.volumePerRadian(1)=2*bad.geometryBinding.volumePerRadian(1);reject('volume-binding',mesh,bad,boundary);
bad=result;bad.gas.mu=1.01*bad.gas.mu;bad.gas=rmfield(bad.gas,{'cp','cv','k'});reject('transport-mu',mesh,bad,boundary);
bad=result;bad.gas.Pr=1.01*bad.gas.Pr;bad.gas=rmfield(bad.gas,{'cp','cv','k'});reject('transport-Pr',mesh,bad,boundary);
bad=result;bad.finalInfo.wall.heatFluxIntoWall(1)=bad.finalInfo.wall.heatFluxIntoWall(1)+1;reject('changed-wall-q',mesh,bad,boundary);
bad=result;bad.finalInfo.wall.pressure(1)=bad.finalInfo.wall.pressure(1)+1;reject('changed-wall-p',mesh,bad,boundary);
badBoundary=boundary;badBoundary.jMin.temperature=312;reject('changed-wall-Tw',mesh,result,badBoundary);
badBoundary=boundary;badBoundary.iMax.primitive(1)=1.01*badBoundary.iMax.primitive(1);reject('changed-reservoir',mesh,result,badBoundary);
assert(~exist(badFolder,'dir')&&strcmp(sourceSHA,fileSHA(evidenceFile)));
report=struct('status','passed-actual-steady-TMR-surface-replay-geometry-field-data-audit', ...
    'sourceEvidenceSHA256',sourceSHA,'steadyResidual',result.residualNorm,'sourceAcceptedSteps',result.acceptedSteps, ...
    'triangles',nf,'vertices',size(V,1),'sourceWallFaces',40,'refusals',refusals,'view',view, ...
    'trueCurvedWallAreaM2',trueArea,'facetedWallAreaM2',sum(facetArea), ...
    'trueCurvedWallHeatW',trueHeat,'facetedWallHeatW',facetHeat, ...
    'sharedUISourceUnchanged',true,'browserVisuallyInspected',false,'webGLRuntimeTested',false, ...
    'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'tmr_steady_view_evidence.mat'),'report','sourceSHA','sourceUIHashes');
fid=fopen(fullfile(folder,'tmr_steady_view_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('TMR steady view passed: 40 walls -> %d triangles, curve area%.12g/facet%.12g m2, true heat%.12g/facet%.12g W; %d refusals; no browser/production claim.\n', ...
    nf,trueArea,sum(facetArea),trueHeat,facetHeat,numel(refusals));
    function reject(label,currentMesh,currentResult,currentBoundary)
        try
            exportTmrAxisymmetricStudy(currentMesh,currentResult,currentBoundary,badFolder,opts);
        catch exception
            assert(startsWith(exception.identifier,'AxisFV:TmrView'),'AxisFV:TmrViewTestReject','非预期拒绝%s：%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFV:TmrViewTestNoRefusal','损坏输入%s未被拒绝。',label);
    end
end
function f=getField(data,key)
f=data.fields(strcmp({data.fields.key},key));assert(isscalar(f),'AxisFV:TmrViewTestField','字段缺失或重复。');
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
