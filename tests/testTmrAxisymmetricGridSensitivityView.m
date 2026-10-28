function report=testTmrAxisymmetricGridSensitivityView(folder,evidenceFile)
%TESTTMRAXISYMMETRICGRIDSENSITIVITYVIEW 四档官方嵌套网格的真实SA稳态科研页核对。
% 模块一：仅读取已收敛五量检查点，重用完整SA查看测试及既有35项拒绝。
% 模块二：独立核对完整官方原节点、采样步长、曲面积/热账与所有真实墙段。
% 模块三：每个显示面保留Dirichlet墙热原audit和gt/dn/gn来源，未知温度保持null。
% 模块四：同尺寸改节点、漏墙/漏audit及不支持尺寸须在写页面前拒绝。
% 本测试不推进时间、不生成伪稳态，不宣称网格无关、实验或生产区域验收。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'tests'),fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v67_validation','normal_view');end
if nargin<2
    evidenceFile=pick({fullfile(root,'outputs','v67_validation','normal','tmr_sensitivity_evidence.mat'), ...
        fullfile(root,'validation','v67','normal','tmr_sensitivity_evidence.mat')});
end
assert(exist(evidenceFile,'file')==2,'AxisFV:TmrGridViewEvidence','需要实际收敛的五量网格敏感性检查点。');
sourceSHA=fileSHA(evidenceFile);source=load(evidenceFile,'mesh','implicit','boundary','report');
mesh=source.mesh;result=source.implicit;boundary=source.boundary;
assert(result.steadyConverged&&source.report.steadyConverged&&result.residualNorm<=2e-6&& ...
    isfield(result.finalInfo,'wallThermalGradient')&&isfield(boundary.jMin,'thermalGradient')&& ...
    strcmp(boundary.jMin.thermalGradient,'dirichlet-tangent'), ...
    'AxisFV:TmrGridViewEvidence','必须是真实五量显式墙热策略稳态，预算终态不得导出。');

%% 一、完整原SA查看器合同与四档官方原节点身份
base=testTmrAxisymmetricSAView(folder,evidenceFile);
assert(numel(base.refusals)==35,'AxisFV:TmrGridViewBase','必须完整保留既有35项拒绝。');
jsonText=fileread(base.view.jsonFile);data=jsondecode(jsonText);
saved=load(base.view.snapshotFile,'data','mesh','result','boundary');
assert(isequaln(saved.mesh,mesh)&&isequaln(saved.result,result)&&isequaln(saved.boundary,boundary));
g=loadTmrMach7Grid(3);ni=size(mesh.X,1);nj=size(mesh.X,2);wallCount=ni-1;
strideI=(g.ni-1)/(ni-1);strideJ=(g.nj-1)/(nj-1);
assert(ismember(strideI,[1 2])&&ismember(strideJ,[1 2]));
sourceI=(1:strideI:g.ni)';sourceJ=(1:strideJ:g.nj)';
assert(isequal(mesh.X,g.X(sourceI,sourceJ))&&isequal(mesh.R,g.R(sourceI,sourceJ)));
id=data.gridIdentity;
assert(strcmp(id.kind,'official-TMR-level3-nested-original-nodes-v1')&& ...
    id.ni==ni&&id.nj==nj&&id.cellCount==mesh.cellCount&&mesh.cellCount==(ni-1)*(nj-1)&& ...
    isequaln(id.source,g.source)&&id.streamwiseStride==strideI&&id.normalStride==strideJ&& ...
    isequal(id.sourceNodeIndicesI,sourceI)&&isequal(id.sourceNodeIndicesJ,sourceJ)&& ...
    id.officialNodesPreserved&&~id.reinterpolated&&id.boundaryExtentUnchanged&& ...
    id.geometryOnly&&~id.externalCFDSolutionImported);
bindings=data.researchIdentity.currentSourceBindings;
for key={'regionalFlow/axisymmetric/loadTmrMach7Grid.m','regionalFlow/axisymmetric/axisymmetricWallThermalGradient.m'}
    entry=bindings(strcmp({bindings.path},key{1}));
    assert(isscalar(entry)&&strcmp(entry.sha256,fileSHA(fullfile(root,key{1}))));
end
assert(numel(unique(data.wallBinding.sourceWallIndex))==wallCount&& ...
    data.surfaceIntegralAudit.sourceWallFaces==wallCount&&data.surfaceIntegralAudit.endcapFaces==0);

%% 二、完整墙热audit逐值保存，独立用真实几何重算法向及切向约束
audit=result.finalInfo.wallThermalGradient;w=result.finalInfo.wall;
assert(isequaln(saved.data.wallThermalGradient,audit));assertJSONAudit(data.wallThermalGradient,audit);
ids=audit.faceID;o=mesh.faces.owner(ids);normal=mesh.faces.normal(ids,:);
tangent=[-normal(:,2) normal(:,1)];delta=mesh.faces.centroid(ids,:)-mesh.volumeCentroid(o,:);
dn=sum(delta.*normal,2);dt=sum(delta.*tangent,2);gradient=audit.gradientXR;
assert(numel(ids)==wallCount&&numel(unique(ids))==wallCount&& ...
    isequal(sort(ids),find(mesh.faces.boundarySide==3))&&isequal(audit.ownerCell,o)&& ...
    all(audit.tangentialGradientInputKm==0)&&all(dn>0)&&~audit.cellStateClipped);
assertNumeric(audit.normalXR,normal);assertNumeric(audit.tangentXR,tangent);
assertNumeric(audit.normalDistanceM,dn);assertNumeric(audit.tangentialOffsetM,dt);
P=axisymmetricFVState(result.Q(:,1:4),result.gas);assertNumeric(audit.cellTemperatureK,P(o,4));
tangentDefect=sum(gradient.*tangent,2)-audit.tangentialGradientInputKm;
temperatureDefect=sum(gradient.*delta,2)-(audit.wallTemperatureInputK-P(o,4));
tangentBudget=256*eps*(sum(abs(gradient.*tangent),2)+abs(audit.tangentialGradientInputKm));
temperatureBudget=256*eps*(sum(abs(gradient.*delta),2)+abs(audit.wallTemperatureInputK)+abs(P(o,4)));
assert(all(abs(tangentDefect)<=tangentBudget)&&all(abs(temperatureDefect)<=temperatureBudget));
normalGradient=sum(gradient.*normal,2);assertNumeric(audit.normalGradientKm,normalGradient);
[present,wallAuditIndex]=ismember(w.faceID,ids);assert(all(present));
qFromAudit=-w.thermalConductivity.*normalGradient(wallAuditIndex);
heatBudget=256*eps*(abs(qFromAudit)+abs(w.heatFluxIntoWall));
assert(all(abs(qFromAudit-w.heatFluxIntoWall)<=heatBudget));

%% 三、每个显示面可追溯至真实墙面热梯度；未知Te/Tr/壁温解不得借值填充
map=data.wallBinding.sourceWallIndex;b=data.wallThermalBinding;index=wallAuditIndex(map);
assert(strcmp(b.policy,'dirichlet-tangent')&&isequal(b.sourceAuditIndex,index)&& ...
    isequal(b.sourceFaceID,w.faceID(map))&&isequal(b.ownerCell,w.ownerCell(map))&& ...
    isequal(data.wallBinding.sourceFaceID,b.sourceFaceID));
assertNumeric(b.normalDistanceM,dn(index));assertNumeric(b.tangentialOffsetM,dt(index));
assertNumeric(b.tangentialGradientInputKm,audit.tangentialGradientInputKm(index));
assertNumeric(b.normalGradientKm,audit.normalGradientKm(index));assertNumeric(b.gradientXR,gradient(index,:));
assertNumeric(b.ownerTemperatureProxyK,P(o(index),4));
assertNumeric(b.wallTemperatureInputK,audit.wallTemperatureInputK(index));
q=getField(data,'heatFlux');p=getField(data,'pressure');tw=getField(data,'wallInput');
assert(all(abs(q.values*1000-qFromAudit(map))<=heatBudget(map)+128*eps*abs(qFromAudit(map)))&& ...
    all(abs(p.values*1000-w.pressure(map))<=256*eps*abs(w.pressure(map)))&& ...
    isequal(tw.values,w.temperatureInput(map))&&all(q.status==1)&&all(p.status==1)&&all(tw.status==2));
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    field=getField(data,key{1});assert(all(isnan(field.values))&&all(field.status==0));
end
assert(contains(jsonText,'null')&&~contains(jsonText,':NaN')&&~contains(jsonText,':Infinity'));
for j=1:wallCount
    assert(contains(q.scopeNames{j},'显式Dirichlet切向约束')&&contains(q.scopeNames{j},'gt=')&& ...
        contains(q.scopeNames{j},'dn=')&&contains(q.scopeNames{j},'单点代理')&& ...
        contains(q.scopeNames{j},'SA标量墙梯度未改'));
end
assert(data.researchIdentity.wallThermalGradientReplayed&& ...
    strcmp(data.resultCondition.wallThermalGradientPolicy,'dirichlet-tangent')&& ...
    contains(q.methodNames{1},'gn=(Tw-Towner-gt*dt)/dn')&& ...
    ~data.researchIdentity.physicalSolution&&~data.researchIdentity.eligibleForSurfaceCoverage&& ...
    ~data.researchIdentity.experimentalValidation&&data.researchIdentity.newProductionSolvedAreaM2==0);

%% 四、不能以相同尺寸、缺失墙段或不完整audit绕过原节点及来源合同
refusals=struct('label',{},'identifier',{});badFolder=[tempname(folder) '_grid_must_not_export'];
badX=mesh.X;badX(2,2)=badX(2,2)+1e-8;
badMesh=buildAxisymmetricFVGeometry(badX,mesh.R);reject('same-shape-rebuilt-nonofficial-node',badMesh,result);
bad=result;bad.finalInfo.wall.faceID=bad.finalInfo.wall.faceID(2:end);
reject('omitted-source-wall-face',mesh,bad);
bad=result;bad.finalInfo.wallThermalGradient.faceID=bad.finalInfo.wallThermalGradient.faceID(2:end);
reject('omitted-thermal-audit-face',mesh,bad);
badMesh=buildAxisymmetricFVGeometry(mesh.X(:,1:end-1),mesh.R(:,1:end-1));
reject('unsupported-node-shape',badMesh,result);
assert(numel(refusals)==4&&~exist(badFolder,'dir')&&strcmp(sourceSHA,fileSHA(evidenceFile)));
report=struct('status','passed-official-nested-TMR-SA-grid-sensitivity-view', ...
    'sourceEvidenceSHA256',sourceSHA,'gridIdentity',id,'baseSAView',base,'baseRefusalCount',numel(base.refusals), ...
    'additionalRefusals',refusals,'sourceWallCount',wallCount,'displayTriangleCount',numel(map), ...
    'completeWallThermalAuditSavedAndReplayed',true,'perDisplayFaceThermalProvenanceChecked',true, ...
    'maximumTangentialDefectKm',max(abs(tangentDefect)), ...
    'maximumTemperatureConstraintDefectK',max(abs(temperatureDefect)), ...
    'maximumSignedHeatDefectWm2',max(abs(qFromAudit-w.heatFluxIntoWall)), ...
    'sourceBindings',sourceBindings(root),'view',base.view,'physicalSolution',false, ...
    'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0, ...
    'gridIndependenceDemonstrated',false,'browserVisuallyInspected',false,'webGLRuntimeTested',false);
save(fullfile(folder,'tmr_grid_sensitivity_view_evidence.mat'),'report','sourceSHA','audit');
fid=fopen(fullfile(folder,'tmr_grid_sensitivity_view_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('TMR official nested SA view: %d x %d nodes, %d walls, %d triangles; 35+4 refusals; no physical/experimental/grid-independence claim.\n', ...
    ni,nj,wallCount,numel(map));
    function reject(label,currentMesh,currentResult)
        try
            exportTmrAxisymmetricSAStudy(currentMesh,currentResult,boundary,badFolder, ...
                struct('azimuthSegments',48,'evidenceSource',evidenceFile));
        catch exception
            assert(startsWith(exception.identifier,'AxisFV:TmrSAView'), ...
                'AxisFV:TmrGridViewRefusal','非法合同%s得到非预期拒绝%s。',label,exception.identifier);
            assert(~exist(badFolder,'dir'),'AxisFV:TmrGridViewWrite','非法数据已生成页面。');
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFV:TmrGridViewNoRefusal','非法网格/来源%s没有拒绝。',label);
    end
end

function path=pick(candidates)
for j=1:numel(candidates),if exist(candidates{j},'file')==2,path=candidates{j};return;end,end
error('AxisFV:TmrGridViewEvidence','缺少实际稳态证据：%s。',strjoin(candidates,' / '));
end
function field=getField(data,key)
field=data.fields(strcmp({data.fields.key},key));assert(isscalar(field));
end
function assertJSONAudit(actual,expected)
assert(isequal(sort(fieldnames(actual)),sort(fieldnames(expected))));
for key=fieldnames(expected)'
    a=actual.(key{1});b=expected.(key{1});
    if isnumeric(b),assertNumeric(a,b);else,assert(isequaln(a,b));end
end
end
function assertNumeric(actual,expected)
% 模块用途：各有单位量只使用自身浮点舍入预算，不用无单位max(1)放宽零值。
assert(isnumeric(actual)&&isreal(actual)&&isequal(size(actual),size(expected))&& ...
    all(isfinite(actual),'all')&&all(abs(actual-expected)<=128*eps*(abs(actual)+abs(expected)),'all'), ...
    'AxisFV:TmrGridViewNumeric','导出数值或几何约束超出独立舍入预算。');
end
function records=sourceBindings(root)
paths={'tests/testTmrAxisymmetricGridSensitivityView.m','tests/testTmrAxisymmetricSAView.m', ...
    'regionalFlow/axisymmetric/exportTmrAxisymmetricSAStudy.m','regionalFlow/axisymmetric/loadTmrMach7Grid.m', ...
    'regionalFlow/axisymmetric/axisymmetricWallThermalGradient.m'};records=struct('path',{},'sha256',{});
for j=1:numel(paths),records(j)=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
