function report=testTmrAxisymmetricSAView(folder,evidenceFile)
%TESTTMRAXISYMMETRICSAVIEW 真实TMR稳态墙面/三维数据/面积热账定向核对。
% 模块一：只读真正五量SA稳态检查点，完整RANS重放并保留原残差/五量BE门槛。
% 模块二：独立复算旋转真面积、锥面参数三角面积和显示三角面积/热积分。
% 模块三：SA来源/Prt、给定N/nu=5与墙muT=0；未知温度/源face/伪时间分开。
% 模块四：篡改气体、边界、墙热、几何或稳态合同均须在生成页面前拒绝。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v65_validation','tmr_sa_view');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v65_validation','tmr_sa_continue1','tmr_sa_evidence.mat');
    if ~exist(evidenceFile,'file'),evidenceFile=fullfile(root,'validation','v65','tmr_sa_continue1','tmr_sa_evidence.mat');end
end
assert(exist(evidenceFile,'file')==2,'AxisFV:TmrSAViewTestEvidence','缺少实际五量SA稳态证据，未稳态检查点不得冒充。');
sourceSHA=fileSHA(evidenceFile);source=load(evidenceFile,'mesh','implicit','boundary','report');
mesh=source.mesh;result=source.implicit;boundary=source.boundary;
ni=size(mesh.X,1);wallCount=ni-1;
assert(result.steadyConverged&&result.residualNorm<=2e-6&&source.report.steadyConverged);
sourceFiles={'web/index.html','web/viewer.css','web/viewer.js','web/section.js','web/export.js'};
sourceUIHashes=cellfun(@(p)fileSHA(fullfile(root,p)),sourceFiles,'UniformOutput',false);
opts=struct('azimuthSegments',48,'evidenceSource',evidenceFile);
view=exportTmrAxisymmetricSAStudy(mesh,result,boundary,fullfile(folder,'web'),opts);
data=jsondecode(fileread(view.jsonFile));w=result.finalInfo.wall;map=data.wallBinding.sourceWallIndex;
V=data.vertices;F=data.faces+1;ns=opts.azimuthSegments;nf=size(F,1);
assert(nf==wallCount*2*ns&&size(V,1)==ni*ns&&numel(data.faceArea)==nf&& ...
    isequal(sort(unique(map)),(1:wallCount)'));
q=getField(data,'heatFlux');p=getField(data,'pressure');tw=getField(data,'wallInput');
assert(isequal(q.values*1000,w.heatFluxIntoWall(map))||max(abs(q.values*1000-w.heatFluxIntoWall(map)))<1e-10);
assert(max(abs(p.values*1000-w.pressure(map)))<1e-10&&isequal(tw.values,w.temperatureInput(map)));
assert(all(q.status==1)&&all(p.status==1)&&all(tw.status==2));
assert(isequal(data.wallBinding.sourceFaceID,w.faceID(map))&&isequal(data.wallBinding.ownerCell,w.ownerCell(map)));
assert(all(q.methodIndex==0)&&contains(q.methodNames{1},result.finalInfo.method)&& ...
    isequal(q.scopeIndex,map-1)&&numel(q.scopeNames)==wallCount&&contains(q.scope,'不缩放signed q'));

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
sourceCurved=accumarray(map,expectedCurved,[wallCount 1]);assert(max(abs(sourceCurved-2*pi*w.areaPerRadian))<1e-14);
cone=abs(dr(segment))>1e-10;cylinder=~cone;
assert(any(cone)&&any(cylinder));
assert(any(abs(expectedCurved(cone)-2*pi*w.areaPerRadian(map(cone))/(2*ns))>1e-8), ...
    'AxisFV:TmrSAViewTestArea','锥面参数对角真面积被错误套用等分圆柱公式。');
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
assert(all(sum(normals.*fluid,2)>.995),'AxisFV:TmrSAViewTestNormal','绘图面法向不对应源流体外法向。');
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

%% 三、SA五量与给定湍流边界身份；墙分子系数与参考系数不能混写
assert(size(result.Q,2)==5&&data.researchIdentity.turbulenceTransportImplemented&& ...
    ~data.researchIdentity.turbulenceValidated&&contains(data.researchIdentity.kind,'SA-noft2')&& ...
    strcmp(data.turbulence.modelID,'SA-noft2-nonnegative-AJS2012-v1')&& ...
    data.turbulence.turbulentPrandtl==.9&&data.turbulence.constants.turbulentPrandtl==.9&& ...
    data.turbulence.densityGradientCorrectionIncluded&&~data.turbulence.duplicateScalarAxisSourceAdded&& ...
    ~data.turbulence.extraCompressibilityCorrection&&~data.turbulence.wallFunctionUsed&& ...
    ~data.turbulence.experimentalValidation&&~data.turbulence.eligibleForSurfaceCoverage);
assert(numel(data.turbulence.boundaryInputs)==3&&contains(data.turbulence.inletVariableMeaning,'not muT/mu=5'));
for j=1:numel(data.turbulence.boundaryInputs)
    input=data.turbulence.boundaryInputs(j);bc=boundary.(input.side);gas=result.gas;
    mu=gas.transport.referenceViscosity*(bc.primitive(4)/gas.transport.referenceTemperature)^1.5* ...
        (gas.transport.referenceTemperature+gas.transport.sutherlandTemperature)/(bc.primitive(4)+gas.transport.sutherlandTemperature);
    nu=mu/bc.primitive(1);chi=bc.nuTilde/nu;muRatio=chi^4/(chi^3+7.1^3);
    assert(abs(input.nuTildeToLaminarKinematicViscosity/5-1)<1e-12&& ...
        input.nuTildeInputM2s==bc.nuTilde&&abs(input.laminarKinematicViscosityM2s/nu-1)<1e-14&& ...
        abs(input.turbulentToLaminarViscosity-muRatio)<1e-14&& ...
        abs(input.turbulentToLaminarViscosity-1.2942343413175514)<2e-12&& ...
        contains(input.valueKind,'prescribed')&&contains(input.valueKind,'not measured'));
end
assert(all(data.turbulence.wallNuTilde==0)&&all(data.turbulence.wallTurbulentViscosityPaS==0)&& ...
    all(data.turbulence.wallTurbulentConductivityWmK==0)&& ...
    all(data.wallBinding.sourceTurbulentViscosityPaS==0));
assert(isequal(data.wallBinding.sourceEffectiveViscosityPaS,w.viscosity(map))&& ...
    isequal(data.wallBinding.sourceEffectiveConductivityWmK,w.thermalConductivity(map))&& ...
    isequal(data.wallBinding.sourceViscosityPaS,w.viscosity(map))&& ...
    isequal(data.wallBinding.sourceThermalConductivityWmK,w.thermalConductivity(map)));
assert(data.transport.referenceViscosityPaS==result.gas.mu&& ...
    data.transport.referenceTemperatureK==result.gas.transport.referenceTemperature&& ...
    data.resultCondition.turbulentPrandtl==.9&& ...
    data.resultCondition.globalBEFiveQuantityBalanceRelative==result.globalBalanceRelative&& ...
    data.resultCondition.conservationTolerance==5e-10&&numel(data.resultCondition.quantityOrder)==5);
assert(contains(data.caseNote,'SA-noft2-nonnegative-AJS2012-v1')&&contains(data.caseNote,'Prt=.90')&& ...
    contains(data.caseNote,'并非muT/mu=5')&&contains(data.caseNote,'墙N=0、muT=0')&& ...
    contains(q.methodNames{1},'墙热未乘增强因子')&&contains(fileread(view.htmlFile),'SA-noft2五量数值研究'));
bindings={data.researchIdentity.currentSourceBindings.path};
assert(all(ismember({'regionalFlow/axisymmetric/axisymmetricRANSResidual.m', ...
    'regionalFlow/axisymmetric/axisymmetricSADiscretization.m','regionalFlow/axisymmetric/axisymmetricSAModel.m', ...
    'regionalFlow/axisymmetric/axisymmetricSAWallDistance.m'},bindings)));

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
bad=result;bad.Q=bad.Q(:,1:4);reject('four-quantity-not-SA',mesh,bad,boundary);
bad=result;bad.Q(1,5)=-abs(bad.Q(1,5));reject('negative-SA-state',mesh,bad,boundary);
bad=result;bad.Q(1,5)=1.01*bad.Q(1,5);reject('changed-fifth-state',mesh,bad,boundary);
bad=result;bad.options.conservationTolerance=5e-8;reject('loosened-five-quantity-ledger',mesh,bad,boundary);
bad=result;bad.options.nonlinearRelativeTolerance=1e-6;reject('loosened-nonlinear-tolerance',mesh,bad,boundary);
bad=result;bad.initialConservedIntegral(5)=bad.initialConservedIntegral(5)+1;reject('changed-fifth-initial-integral',mesh,bad,boundary);
bad=result;bad.accumulatedSourceIntegral(5)=bad.accumulatedSourceIntegral(5)+1;reject('changed-fifth-source-ledger',mesh,bad,boundary);
bad=result;bad.history.stepDefect(end,5)=1;reject('changed-fifth-step-ledger',mesh,bad,boundary);
bad=result;bad.finalInfo.sa.modelID='SA-standard-with-ft2';reject('different-SA-model',mesh,bad,boundary);
bad=result;bad.finalInfo.sa.cellModel.constants.turbulentPrandtl=.85;reject('different-Prt',mesh,bad,boundary);
bad=result;bad.finalInfo.sa.faceMuT(w.faceID(1))=1e-5;reject('wall-nonzero-eddy-viscosity',mesh,bad,boundary);
bad=result;bad.finalInfo.sa.effectiveFaceK(w.faceID(1))=1.01*w.thermalConductivity(1);reject('wall-heat-multiplier',mesh,bad,boundary);
bad=result;bad.finalInfo.sa.sourceDensityCorrection(:)=0;reject('omitted-density-gradient-correction',mesh,bad,boundary);
bad=result;bad.finalInfo=rmfield(bad.finalInfo,'sa');reject('missing-SA-transport',mesh,bad,boundary);
bad=result;bad.finalInfo.turbulenceValidated=true;reject('promoted-turbulence-validation',mesh,bad,boundary);
badBoundary=boundary;badBoundary.iMin=rmfield(badBoundary.iMin,'nuTilde');reject('missing-given-SA-inlet',mesh,result,badBoundary);
badBoundary=boundary;badBoundary.iMax.nuTilde=1.01*badBoundary.iMax.nuTilde;reject('changed-given-SA-backflow',mesh,result,badBoundary);
assert(~exist(badFolder,'dir')&&strcmp(sourceSHA,fileSHA(evidenceFile)));
report=struct('status','passed-actual-steady-five-quantity-SA-TMR-surface-audit', ...
    'sourceEvidenceSHA256',sourceSHA,'steadyResidual',result.residualNorm,'sourceAcceptedSteps',result.acceptedSteps, ...
    'triangles',nf,'vertices',size(V,1),'sourceWallFaces',wallCount,'refusals',refusals,'view',view, ...
    'trueCurvedWallAreaM2',trueArea,'facetedWallAreaM2',sum(facetArea), ...
    'trueCurvedWallHeatW',trueHeat,'facetedWallHeatW',facetHeat, ...
    'modelID',data.turbulence.modelID,'fiveQuantityLedgerChecked',true,'givenNuTildeRatio',5, ...
    'wallEddyViscosityZero',true,'molecularAndReferenceTransportSeparated',true, ...
    'sharedUISourceUnchanged',true,'browserVisuallyInspected',false,'webGLRuntimeTested',false, ...
    'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'tmr_sa_view_evidence.mat'),'report','sourceSHA','sourceUIHashes');
fid=fopen(fullfile(folder,'tmr_sa_view_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('TMR five-quantity SA steady view passed: %d walls -> %d triangles, curve area%.12g/facet%.12g m2, true heat%.12g/facet%.12g W; %d refusals; no browser/production claim.\n', ...
    wallCount,nf,trueArea,sum(facetArea),trueHeat,facetHeat,numel(refusals));
    function reject(label,currentMesh,currentResult,currentBoundary)
        try
            exportTmrAxisymmetricSAStudy(currentMesh,currentResult,currentBoundary,badFolder,opts);
        catch exception
            assert(startsWith(exception.identifier,'AxisFV:TmrSAView'),'AxisFV:TmrSAViewTestReject','非预期拒绝%s：%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFV:TmrSAViewTestNoRefusal','损坏输入%s未被拒绝。',label);
    end
end
function f=getField(data,key)
f=data.fields(strcmp({data.fields.key},key));assert(isscalar(f),'AxisFV:TmrSAViewTestField','字段缺失或重复。');
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
