function report=testTmrAxisymmetricWallThermalView(folder,evidenceFile)
%TESTTMRAXISYMMETRICWALLTHERMALVIEW 显式斜壁热梯度稳态的三维导出与溯源核对。
% 模块一：先执行既有五量SA查看器全部35项拒绝、真实几何/热账与未知温度检查。
% 模块二：完整audit逐值重放并映射至每个显示面，独立检查切向/位移温差和signed墙热。
% 模块三：中文探针保留gt/dn及守恒均值单点代理限制；墙温输入仍不是求解壁温。
% 模块四：改回legacy、篡改gt/gn或删改来源必须在生成页面前拒绝，不改冻结v65。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v66_validation','tmr_wall_thermal_view');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v66_validation','tmr_wall_thermal','tmr_wall_thermal_evidence.mat');
    if ~exist(evidenceFile,'file')
        evidenceFile=fullfile(root,'validation','v66','tmr_wall_thermal','tmr_wall_thermal_evidence.mat');
    end
end
assert(exist(evidenceFile,'file')==2,'AxisFV:TmrWallThermalViewTestEvidence', ...
    '缺少v66实际新方程稳态，不能用旧Q的单次正热流诊断代替。');
sourceSHA=fileSHA(evidenceFile);source=load(evidenceFile,'mesh','implicit','boundary','report');
mesh=source.mesh;result=source.implicit;boundary=source.boundary;
assert(result.steadyConverged&&source.report.steadyConverged&&result.residualNorm<=2e-6&& ...
    isfield(result.finalInfo,'wallThermalGradient')&& ...
    isfield(boundary.jMin,'thermalGradient')&&strcmp(boundary.jMin.thermalGradient,'dirichlet-tangent'), ...
    'AxisFV:TmrWallThermalViewTestEvidence','必须是真正显式Dirichlet切向策略的五量稳态。');

%% 一、复用既有35项拒绝；不复制或削弱v65查看器的门槛
base=testTmrAxisymmetricSAView(folder,evidenceFile);
assert(numel(base.refusals)==35,'AxisFV:TmrWallThermalViewTestBase','既有35项拒绝没有完整执行。');
data=jsondecode(fileread(base.view.jsonFile));
saved=load(base.view.snapshotFile,'data');
audit=result.finalInfo.wallThermalGradient;w=result.finalInfo.wall;
assert(isequaln(saved.data.wallThermalGradient,audit), ...
    'AxisFV:TmrWallThermalViewTestAudit','保存的完整audit不等于实际稳态的原始梯度audit。');
assertJSONAudit(data.wallThermalGradient,audit);
[rate,replay]=axisymmetricRANSResidual(result.Q,mesh,result.gas,boundary);
assert(isequaln(replay,result.finalInfo)&& ...
    max(abs(rate)./result.residualScale,[],'all')==result.residualNorm);

%% 二、从真实几何而非audit的预存向量独立核对温度和法向通量约束
ids=audit.faceID;o=mesh.faces.owner(ids);normal=mesh.faces.normal(ids,:);
tangent=[-normal(:,2) normal(:,1)];delta=mesh.faces.centroid(ids,:)-mesh.volumeCentroid(o,:);
dn=sum(delta.*normal,2);dt=sum(delta.*tangent,2);g=audit.gradientXR;
assert(numel(ids)==40&&numel(unique(ids))==40&&isequal(sort(ids),sort(w.faceID))&& ...
    isequal(audit.ownerCell,o)&&all(audit.tangentialGradientInputKm==0)&&all(dn>0)&&~audit.cellStateClipped);
assertNumeric(audit.normalXR,normal);assertNumeric(audit.tangentXR,tangent);
assertNumeric(audit.normalDistanceM,dn);assertNumeric(audit.tangentialOffsetM,dt);
P=axisymmetricFVState(result.Q(:,1:4),result.gas);
assertNumeric(audit.cellTemperatureK,P(o,4));
tangentDefect=sum(g.*tangent,2)-audit.tangentialGradientInputKm;
temperatureDefect=sum(g.*delta,2)-(audit.wallTemperatureInputK-P(o,4));
tangentBudget=256*eps*(sum(abs(g.*tangent),2)+abs(audit.tangentialGradientInputKm));
temperatureBudget=256*eps*(sum(abs(g.*delta),2)+abs(audit.wallTemperatureInputK)+abs(P(o,4)));
assert(all(abs(tangentDefect)<=tangentBudget)&&all(abs(temperatureDefect)<=temperatureBudget), ...
    'AxisFV:TmrWallThermalViewTestConstraint','真实墙几何下的切向导数或cell-to-wall温差约束未闭合。');
normalGradient=sum(g.*normal,2);assertNumeric(audit.normalGradientKm,normalGradient);
[present,wallAuditIndex]=ismember(w.faceID,ids);assert(all(present));
qFromAudit=-w.thermalConductivity.*normalGradient(wallAuditIndex);
heatBudget=256*eps*(abs(qFromAudit)+abs(w.heatFluxIntoWall));
assert(all(abs(qFromAudit-w.heatFluxIntoWall)<=heatBudget), ...
    'AxisFV:TmrWallThermalViewTestHeat','实际signed墙热不等于完整audit梯度及分子k(Tw)的Fourier通量。');

%% 三、逐显示面保留源face、gt/dn/gn和中文适用限制
map=data.wallBinding.sourceWallIndex;b=data.wallThermalBinding;index=wallAuditIndex(map);
assert(strcmp(b.policy,'dirichlet-tangent')&&isequal(b.sourceAuditIndex,index)&& ...
    isequal(b.sourceFaceID,w.faceID(map))&&isequal(b.ownerCell,w.ownerCell(map))&& ...
    isequal(data.wallBinding.sourceFaceID,b.sourceFaceID));
assertNumeric(b.normalDistanceM,dn(index));assertNumeric(b.tangentialOffsetM,dt(index));
assertNumeric(b.tangentialGradientInputKm,audit.tangentialGradientInputKm(index));
assertNumeric(b.normalGradientKm,audit.normalGradientKm(index));assertNumeric(b.gradientXR,g(index,:));
assertNumeric(b.ownerTemperatureProxyK,P(o(index),4));
assertNumeric(b.wallTemperatureInputK,audit.wallTemperatureInputK(index));
q=data.fields(strcmp({data.fields.key},'heatFlux'));assert(isscalar(q));
assert(all(abs(q.values*1000-qFromAudit(map))<=heatBudget(map)+128*eps*abs(qFromAudit(map))));
assert(strcmp(data.researchIdentity.wallThermalGradientPolicy,'dirichlet-tangent')&& ...
    data.researchIdentity.wallThermalGradientReplayed&& ...
    strcmp(data.resultCondition.wallThermalGradientPolicy,'dirichlet-tangent')&& ...
    contains(data.caseNote,'Dirichlet')&&contains(data.caseNote,'gt=')&&contains(data.caseNote,'dn=')&& ...
    contains(data.caseNote,'单点代理')&&contains(data.transport.evaluationPolicy,'墙热梯度显式')&& ...
    contains(q.methodNames{1},'单点代理')&&contains(q.methodNames{1},'gn=(Tw-Towner-gt*dt)/dn'));
for j=1:40
    assert(contains(q.scopeNames{j},'显式Dirichlet切向约束')&&contains(q.scopeNames{j},'gt=')&& ...
        contains(q.scopeNames{j},'dn=')&&contains(q.scopeNames{j},'单点代理')&& ...
        contains(q.scopeNames{j},'SA标量墙梯度未改'));
end
bindings=data.researchIdentity.currentSourceBindings;
helper='regionalFlow/axisymmetric/axisymmetricWallThermalGradient.m';
entry=bindings(strcmp({bindings.path},helper));
assert(isscalar(entry)&&strcmp(entry.sha256,fileSHA(fullfile(root,helper))));

%% 四、不能改策略沿用新解，也不能篡改audit再保留原热流
refusals=struct('label',{},'identifier',{});badFolder=[tempname(folder) '_must_not_export'];
badBoundary=boundary;badBoundary.jMin.thermalGradient='legacy';
reject('reuse-new-solution-with-legacy-policy',result,badBoundary);
badBoundary=boundary;badBoundary.jMin=rmfield(badBoundary.jMin,'thermalGradient');
reject('omit-explicit-policy-on-new-solution',result,badBoundary);
bad=result;bad.finalInfo.wallThermalGradient.tangentialGradientInputKm(1)=1;
reject('tampered-audit-tangential-derivative',bad,boundary);
bad=result;bad.finalInfo.wallThermalGradient.normalGradientKm(1)=bad.finalInfo.wallThermalGradient.normalGradientKm(1)+1;
reject('tampered-audit-normal-gradient',bad,boundary);
bad=result;bad.finalInfo.wallThermalGradient.gradientXR(1,1)=bad.finalInfo.wallThermalGradient.gradientXR(1,1)+1;
reject('tampered-audit-gradient-vector',bad,boundary);
bad=result;bad.finalInfo=rmfield(bad.finalInfo,'wallThermalGradient');
reject('missing-complete-gradient-audit',bad,boundary);
bad=result;bad.finalInfo.wallThermalGradient.faceID(1)=bad.finalInfo.wallThermalGradient.faceID(2);
reject('duplicate-audit-source-face',bad,boundary);
badBoundary=boundary;badBoundary.jMin.temperatureTangentialGradient=1;
reject('nonzero-given-tangent-on-isothermal-wall',result,badBoundary);
assert(~exist(badFolder,'dir')&&strcmp(sourceSHA,fileSHA(evidenceFile)));

report=struct('status','passed-actual-steady-explicit-Dirichlet-wall-gradient-view', ...
    'sourceEvidenceSHA256',sourceSHA,'baseSAView',base,'baseRefusalCount',numel(base.refusals), ...
    'additionalRefusals',refusals,'completeAuditSavedAndReplayed',true, ...
    'sourceWallCount',40,'displayTriangleCount',numel(map),'signedHeatBindingChecked',true, ...
    'maximumTangentialDefectKm',max(abs(tangentDefect)), ...
    'maximumTemperatureConstraintDefectK',max(abs(temperatureDefect)), ...
    'maximumSignedHeatDefectWm2',max(abs(qFromAudit-w.heatFluxIntoWall)), ...
    'temperatureIsConservativeMeanPointProxy',true,'SAWallGradientChanged',false, ...
    'sourceBindings',sourceBindings(root),'view',base.view,'physicalSolution',false, ...
    'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0, ...
    'browserVisuallyInspected',false,'webGLRuntimeTested',false);
save(fullfile(folder,'tmr_wall_thermal_view_evidence.mat'),'report','sourceSHA','audit');
fid=fopen(fullfile(folder,'tmr_wall_thermal_view_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('TMR explicit wall thermal view: complete audit/40 source walls/%d triangles; %d base + %d new refusals; no browser or physical validation claim.\n', ...
    numel(map),numel(base.refusals),numel(refusals));
    function reject(label,currentResult,currentBoundary)
        try
            exportTmrAxisymmetricSAStudy(mesh,currentResult,currentBoundary,badFolder, ...
                struct('azimuthSegments',48,'evidenceSource',evidenceFile));
        catch exception
            assert(startsWith(exception.identifier,'AxisFV:TmrSAView'), ...
                'AxisFV:TmrWallThermalViewTestRefusal','非法合同%s得到非预期拒绝%s。',label,exception.identifier);
            assert(~exist(badFolder,'dir'),'AxisFV:TmrWallThermalViewTestWrite','非法来源已写出误导页面。');
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFV:TmrWallThermalViewTestNoRefusal','非法墙热来源%s没有拒绝。',label);
    end
end

function assertJSONAudit(actual,expected)
% 模块用途：JSON必须保留全部原audit字段和数值，不能仅保存通过标志或省去旧梯度。
assert(isequal(sort(fieldnames(actual)),sort(fieldnames(expected))), ...
    'AxisFV:TmrWallThermalViewTestAudit','导出audit字段不完整。');
for key=fieldnames(expected)'
    a=actual.(key{1});b=expected.(key{1});
    if isnumeric(b),assertNumeric(a,b);
    else,assert(isequaln(a,b),'AxisFV:TmrWallThermalViewTestAudit','audit非数值字段%s改变。',key{1});end
end
end

function assertNumeric(actual,expected)
% 模块用途：按各自有单位数值的舍入预算核对，零值不借用无单位max(1)容差。
assert(isnumeric(actual)&&isreal(actual)&&isequal(size(actual),size(expected))&& ...
    all(isfinite(actual),'all')&&all(abs(actual-expected)<=128*eps*(abs(actual)+abs(expected)),'all'), ...
    'AxisFV:TmrWallThermalViewTestNumeric','导出数值或几何约束超过独立舍入预算。');
end

function records=sourceBindings(root)
% 模块用途：绑定新增测试、复用的旧测试、导出器及实际墙梯度helper的完整源码。
paths={'tests/testTmrAxisymmetricWallThermalView.m','tests/testTmrAxisymmetricSAView.m', ...
    'regionalFlow/axisymmetric/exportTmrAxisymmetricSAStudy.m', ...
    'regionalFlow/axisymmetric/axisymmetricWallThermalGradient.m'};records=struct('path',{},'sha256',{});
for j=1:numel(paths),records(j)=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));end
end

function digest=fileSHA(path)
% 模块用途：只读证据/源码哈希；不修改旧版检查点或测试文件。
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
