function report=testTmrAxisymmetricWallComparison(folder,evidenceFile)
%TESTTMRAXISYMMETRICWALLCOMPARISON 真实Sutherland稳态壁面实验诊断的定向核对。
% 模块一：无真实MAT时明确等待，不制造假流场、不启动求解器。
% 模块二：独立按实际jMin节点区间核对每点face/owner及当前原p/q逐值映射。
% 模块三：检查24独立压力点/27热流点、重复不加权、SI量及signed热身份。
% 模块四：内存篡改NaN、墙热、边界与原稳态合同均须拒绝；不保存假证据。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'tests'),fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v64_validation','tmr_wall_comparison');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v64_validation','tmr_sutherland','tmr_sutherland_evidence.mat');
    if ~exist(evidenceFile,'file')
        evidenceFile=fullfile(root,'validation','v64','tmr_sutherland','tmr_sutherland_evidence.mat');
    end
end
assert(exist(evidenceFile,'file')==2,'AxisFV:TmrWallTestEvidence', ...
    '需要等待真实Sutherland稳态MAT；本测试不写假数据，也不执行求解。');
sourceSHA=fileSHA(evidenceFile);s=load(evidenceFile,'mesh','gas','boundary','implicit','report');
comparisonSourceSHA=fileSHA(fullfile(root,'tests','compareTmrAxisymmetricWall.m'));
testSourceSHA=fileSHA([mfilename('fullpath') '.m']);
comparison=compareTmrAxisymmetricWall(folder,evidenceFile);b=loadTmrMach7Benchmark;
p=struct2table(comparison.pressure);q=struct2table(comparison.heating);w=s.implicit.finalInfo.wall;
assert(height(p)==27&&height(q)==27&&nnz(p.independentComparisonPoint)==24&& ...
    isequal(p.sourceRow,(1:27)')&&isequal(q.sourceRow,(1:27)')&& ...
    isequal(find(~p.independentComparisonPoint),[25;26;27]));
assert(isequal(p.experimentalRatio,b.pressure.ratio)&&isequal(q.experimentalRatio,b.heating.ratio)&& ...
    isequal(p.experimentalPressurePa,p.experimentalRatio*576)&& ...
    isequal(q.experimentalReportedHeatingMagnitudeWm2,q.experimentalRatio*9300));
assert(isequal(p.projectX_m,-p.sourceX_m)&&isequal(q.projectX_m,-q.sourceX_m));
assert(max(abs(p.projectX_m(p.S_m>=0)+p.S_m(p.S_m>=0)*cosd(20)))<1e-15&& ...
    max(abs(q.projectX_m(q.S_m>=0)+q.S_m(q.S_m>=0)*cosd(20)))<1e-15);

%% 一、独立按真实墙面centroid及节点X范围核对，不调用比较器内部采样函数
x=s.mesh.X(:,1);r=s.mesh.R(:,1);centroid=w.centroid(:,1);
for kind=1:2
    if kind==1,t=p;else,t=q;end
    for j=1:height(t)
        k=find(t.projectX_m(j)<=x(1:end-1)&t.projectX_m(j)>x(2:end));
        if isempty(k)&&t.projectX_m(j)==x(end),k=40;end
        if isempty(k)
            assert(t.missingMask(j)&&isnan(t.faceID(j))&&isnan(t.ownerCell(j)));
            continue;
        end
        index=find(centroid<=x(k)&centroid>=x(k+1));assert(isscalar(index));
        assert(~t.missingMask(j)&&t.sourceSegment(j)==k&&t.sourceWallIndex(j)==index&& ...
            t.faceID(j)==w.faceID(index)&&t.ownerCell(j)==w.ownerCell(index)&& ...
            t.faceUpstreamXM(j)==x(k)&&t.faceDownstreamXM(j)==x(k+1));
        expectedRadius=r(k)+(t.projectX_m(j)-x(k))/(x(k+1)-x(k))*(r(k+1)-r(k));
        assert(t.actualGridRadiusAtXM(j)==expectedRadius);
        if kind==1
            assert(t.computedPressurePa(j)==w.pressure(index)&&t.computedRatio(j)==w.pressure(index)/576);
        else
            assert(t.computedSignedHeatFluxIntoWallWm2(j)==w.heatFluxIntoWall(index)&& ...
                t.computedSignedRatio(j)==w.heatFluxIntoWall(index)/9300);
        end
    end
end
assert(~any(p.missingMask)&&~any(q.missingMask),'AxisFV:TmrWallTestRange','当前真实网格应覆盖全部公开测点。');
assert(isequal(p.sourceWallIndex(24:27),repmat(p.sourceWallIndex(24),4,1))&& ...
    isequal(p.computedPressurePa(24:27),repmat(p.computedPressurePa(24),4,1)));
assert(isequal(p.differencePa,p.computedPressurePa-p.experimentalPressurePa)&& ...
    isequal(q.signedMinusReportedMagnitudeWm2,q.computedSignedHeatFluxIntoWallWm2-q.experimentalReportedHeatingMagnitudeWm2));
mask=p.independentComparisonPoint&~p.missingMask;v=p.relativeDifference(mask);
assert(comparison.summary.pressure.comparedPointCount==24&& ...
    comparison.summary.pressure.meanRelativeDifference==mean(v)&& ...
    comparison.summary.pressure.meanAbsoluteRelativeDifference==mean(abs(v))&& ...
    comparison.summary.pressure.rootMeanSquareRelativeDifference==sqrt(mean(v.^2))&& ...
    comparison.summary.heatSignedVersusReportedPositiveMagnitude.comparedPointCount==27);
assert(isequal(q.positiveHeatingMagnitudeAvailable,q.computedSignedHeatFluxIntoWallWm2>=0)&& ...
    isequaln(q.computedPositiveHeatingMagnitudeWm2(q.positiveHeatingMagnitudeAvailable), ...
    q.computedSignedHeatFluxIntoWallWm2(q.positiveHeatingMagnitudeAvailable))&& ...
    all(isnan(q.computedPositiveHeatingMagnitudeWm2(~q.positiveHeatingMagnitudeAvailable))));

%% 二、报告保持模型差异、来源与非统计实验误差身份，CSV/JSON/MAT可逐点重载
assert(strcmp(comparison.status,'diagnostic-laminar-vs-turbulent-experiment-mismatches')&& ...
    comparison.sourceHashesVerified&&comparison.currentSpatialResidualReplayPassed&&comparison.sourceSnapshotMatchesEvidenceFile&& ...
    strcmp(comparison.evidenceSHA256,sourceSHA)&&comparison.steadyTolerance==2e-6&&comparison.steadyRequiredSteps==3);
assert(~isfield(comparison,'accuracyPass')&&~isfield(comparison.summary,'accuracyPass')&& ...
    ~comparison.experimentalValidation&&~comparison.physicalSolution&&~comparison.eligibleForSurfaceCoverage&& ...
    comparison.newProductionSolvedAreaM2==0&&~comparison.pseudoTimeIsPhysical&&comparison.physicalTimeAdvancedS==0);
assert(comparison.uncertainty.relativeSurfacePressureEstimate==.1&& ...
    comparison.uncertainty.relativeSurfaceHeatEstimate==.1&& ...
    ~comparison.uncertainty.isStatisticalConfidenceInterval&&isempty(comparison.uncertainty.confidenceLevel)&& ...
    ~comparison.uncertainty.ratioUncertaintyIndependenceAssumed);
assert(isequal(p.experimentalEstimatedLowerPa,.9*p.experimentalPressurePa)&& ...
    isequal(q.experimentalEstimatedUpperWm2,1.1*q.experimentalReportedHeatingMagnitudeWm2));
m=comparison.modelDifferences;assert(isequaln(m.gas,s.gas)&&isequaln(m.boundaryInput,s.boundary)&& ...
    contains(m.flowIdentity,'laminar')&&contains(m.flowIdentity,'turbulent experiment')&& ...
    contains(m.openBoundaryIdentity,'not source NMF extrapolation outlet')&&m.nominalFlareAngleDeg==20&& ...
    all(abs(m.actualOriginalGridFlareAnglesDeg-19.999867)<1e-5)&& ...
    all(m.actualMinusNominalAngleDeg~=0)&&~m.upstreamProfileAgreementVerified&&~m.gridAndDomainSensitivityVerified);
tr=axisymmetricFVTransport(s.boundary.iMin.primitive(4),s.gas);v0=s.boundary.iMin.primitive;
assert(m.derivedInletRePerM==v0(1)*norm(v0(2:3))/tr.mu&& ...
    m.relativeInletReDifference==m.derivedInletRePerM/5.706e6-1&& ...
    m.experimentalLocalReference.pressure_Pa==576&&m.experimentalLocalReference.Mach==7.05);
bindings=comparison.sourceBindings;assert(ismember('tests/compareTmrAxisymmetricWall.m',{bindings.path}));
for j=1:numel(bindings),assert(strcmp(bindings(j).sha256,fileSHA(fullfile(root,bindings(j).path))));end
csvP=readtable(comparison.files.pressureCSV);csvQ=readtable(comparison.files.heatingCSV);
assert(height(csvP)==27&&height(csvQ)==27&&isequal(csvP.sourceRow,p.sourceRow)&&isequal(csvQ.faceID,q.faceID));
assert(max(abs(csvP.computedPressurePa-p.computedPressurePa)./p.computedPressurePa)<1e-13&& ...
    max(abs(csvQ.computedSignedHeatFluxIntoWallWm2-q.computedSignedHeatFluxIntoWallWm2)./ ...
    max(1,abs(q.computedSignedHeatFluxIntoWallWm2)))<1e-13);
json=jsondecode(fileread(comparison.files.json));mat=load(comparison.files.mat,'comparison');
assert(isequaln(mat.comparison,comparison)&&strcmp(json.evidenceSHA256,sourceSHA)&& ...
    numel(json.pressure)==27&&numel(json.heating)==27&&nnz([json.pressure.independentComparisonPoint])==24);

%% 三、只改真实证据的内存副本；每个篡改必须在创建输出目录前被拒绝
refusals=struct('label',{},'identifier',{});badFolder=[tempname(folder) '_must_not_compare'];
bad=s.implicit;bad.steadyConverged=false;reject('not-steady',struct('implicit',bad));
bad=s.implicit;bad.options.steadyResidualTolerance=2e-3;reject('loosened-tolerance',struct('implicit',bad));
bad=s.implicit;bad.options.steadyRequiredSteps=1;reject('one-step-steady',struct('implicit',bad));
bad=s.implicit;bad.history.residual(end-1)=3e-6;reject('broken-three-step-history',struct('implicit',bad));
bad=s.implicit;bad.history.residual(end-1)=NaN;reject('nan-history',struct('implicit',bad));
bad=s.implicit;bad.residualNorm=NaN;reject('nan-final-residual',struct('implicit',bad));
bad=s.implicit;bad.Q(1,1)=NaN;reject('nan-state',struct('implicit',bad));
bad=s.implicit;bad.finalInfo.wall.pressure(1)=NaN;reject('nan-pressure',struct('implicit',bad));
bad=s.implicit;bad.finalInfo.wall.heatFluxIntoWall(1)=NaN;reject('nan-heat',struct('implicit',bad));
bad=s.implicit;bad.finalInfo.wall.pressure(1)=bad.finalInfo.wall.pressure(1)+1;reject('stale-pressure',struct('implicit',bad));
bad=s.implicit;bad.finalInfo.wall.heatFluxIntoWall(1)=-bad.finalInfo.wall.heatFluxIntoWall(1);reject('changed-signed-heat',struct('implicit',bad));
bad=s.boundary;bad.jMin.temperature=312;reject('changed-prescribed-wall',struct('boundary',bad));
bad=s.boundary;bad.iMax.primitive(1)=1.01*bad.iMax.primitive(1);reject('changed-open-boundary',struct('boundary',bad));
bad=s.implicit;bad.pseudoTimeIsPhysical=true;reject('physical-pseudo-time',struct('implicit',bad));
bad=s.implicit;bad.eligibleForSurfaceCoverage=true;reject('promoted-coverage',struct('implicit',bad));
bad=s.gas;bad=rmfield(bad,'transport');reject('old-constant-transport',struct('gas',bad));
assert(~exist(badFolder,'dir')&&strcmp(sourceSHA,fileSHA(evidenceFile)));
assert(strcmp(comparisonSourceSHA,fileSHA(fullfile(root,'tests','compareTmrAxisymmetricWall.m')))&& ...
    strcmp(testSourceSHA,fileSHA([mfilename('fullpath') '.m'])),'AxisFV:TmrWallTestSource','比较或测试源码在运行中改变。');
report=struct('status','passed-real-Sutherland-wall-diagnostic-mapping-and-refusals', ...
    'evidenceSHA256',sourceSHA,'comparisonFiles',comparison.files,'sourceBindings',comparison.sourceBindings, ...
    'testSourceSHA256',testSourceSHA,'comparisonSourceSHA256Before',comparisonSourceSHA,'sourceStableAfter',true, ...
    'pressureRawRows',27,'independentPressurePoints',24, ...
    'heatingPoints',27,'wallFaces',40,'refusals',refusals,'summary',comparison.summary, ...
    'solverRun',false,'currentSpatialResidualReplayPassed',true,'experimentalValidation',false, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'tmr_wall_comparison_test.mat'),'report');
fid=fopen(fullfile(folder,'tmr_wall_comparison_test.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report,'ConvertInfAndNaN',true));clear cleanup;
fprintf('TMR wall comparison checks passed: 24 independent pressure / 27 heat, exact current p/q map, %d refusals; no physical-validation pass.\n',numel(refusals));
    function reject(label,overrides)
        try
            compareTmrAxisymmetricWall(badFolder,evidenceFile,overrides);
        catch exception
            assert(startsWith(exception.identifier,'AxisFV:TmrWall'), ...
                'AxisFV:TmrWallTestReject','%s出现非预期异常%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFV:TmrWallTestNoRefusal','篡改%s没有被拒绝。',label);
    end
end
function digest=fileSHA(path)
% 模块用途：只读锁定真实证据和本次测试源码，测试不写入冻结v63或求解MAT。
if exist(path,'file')~=2&&exist([path '.m'],'file')==2,path=[path '.m'];end
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
