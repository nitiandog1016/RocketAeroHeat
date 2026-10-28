function report=exportRocketExactFluxStudy(folder,evidenceRoot,resultFolder)
%EXPORTROCKETEXACTFLUXSTUDY 将真实20度同态通量对照及新方法修正接回统一科研工作台。
% 模块一：旧HLLC末态、新精确通量同态基线及新方法末态分别标识，不拼接残差历史。
% 模块二：原面压力/热流/邻壁气体温度/算子温度斜率共享色标，逐变量保留方法。
% 模块三：给定Tw与未知Te/Tr/求解Tw分开；所有有限研究数值仍不具生产合格资格。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<2
    evidenceRoot=fullfile(root,'outputs','v82_validation');
    if ~exist(fullfile(evidenceRoot,'integrate01','exact_integration_evidence.mat'),'file'),evidenceRoot=fullfile(root,'validation','v82');end
end
if nargin<3,resultFolder='exact_krylov01';end
assert(ismember(resultFolder,{'exact_krylov01','primitive_krylov01'}),'ExactFluxView:Source','请明确选择已保存的真实修正组。');
assert(~exist(folder,'dir'),'ExactFluxView:Overwrite','请使用新的科研导出目录。');
d=load(publishedInput(root,80,fullfile('domain02','shell_domain_evidence.mat')),'domain','fs');
cut=d.domain;fs=d.fs;mesh=cut.inputGeometry;n=size(mesh.F,1);
s=load(publishedInput(root,81,fullfile('arithmetic02','steady_krylov_evidence.mat')),'result','gas','report');
gas=s.gas;scale=s.result.residualScale;infos={pick(s.result.info)};rates={s.result.rate};
sourceSummary=struct('status','frozen-HLLC-source-state-not-recomputed-here', ...
    'acceptedNonlinearIterations',0,'candidateTrials',0,'residualEvaluations',0,'elapsedSeconds',0, ...
    'maximumResidual',s.result.maximumResidual,'rmsMerit',s.result.rmsMerit,'steadyConverged',false, ...
    'physicalTimeAdvanced',0,'BETimeStepsClaimed',0,'sourceHLLCAcceptedIterations',6,'interiorEulerFlux','hllc', ...
    'iterationCounterMeaning','frozen-source-state-not-recomputed-HLLC-iterations','frozenSourceReport',s.report);
summaries={sourceSummary};histories={struct([])};clear s;
s=load(fullfile(evidenceRoot,'integrate01','exact_integration_evidence.mat'),'info','rate','report');
infos{2}=pick(s.info);rates{2}=s.rate;integrationReport=s.report;
summary=s.report;summary.acceptedNonlinearIterations=0;summary.candidateTrials=0;summary.residualEvaluations=1;
summary.steadyConverged=false;summary.physicalTimeAdvanced=0;summary.BETimeStepsClaimed=0;
summary.interiorEulerFlux='exact-godunov';summary.sourceHLLCAcceptedIterations=6;
summary.iterationCounterMeaning='new-equation-baseline-zero-corrections';summaries{2}=summary;histories{2}=struct([]);clear s;
s=load(fullfile(evidenceRoot,resultFolder,'steady_krylov_evidence.mat'),'result','report');
infos{3}=pick(s.result.info);rates{3}=s.result.rate;summaries{3}=s.report;histories{3}=s.result.history;clear s;
labels={'HLLC源状态 · 旧方法六次末态','精确Godunov · 相同Q新方程基线','精确Godunov · 新方程有界修正'};
ids={'hllc-source','exact-baseline','exact-krylov'};policies={'hllc','exact-godunov','exact-godunov'};
if strcmp(resultFolder,'primitive_krylov01')
    labels{3}='精确Godunov · primitive切向有界修正';ids{3}='exact-primitive-krylov';
    summaries{3}.resultEvidenceFolder=resultFolder;
end
mkdir(folder);exportRocketTestGeometryView({mesh},fullfile(folder,'geometry'),{'原19752面真实箭体；X鼻尖正向'});
exportRocketTestGeometryView({mesh,struct('V',cut.nodeLayers(:,:,2),'F',mesh.F),struct('V',cut.nodeLayers(:,:,end),'F',mesh.F)}, ...
    fullfile(folder,'domain'),{'原墙','第一近壁节点层','最外计算边界；不是激波'});
common=struct('pressure',[Inf -Inf],'heatFlux',[Inf -Inf],'wallAdjacentTemperature',[Inf -Inf],'temperatureRate',[Inf -Inf]);
datasets=cell(3,1);entries=struct([]);cv=gas.R/(gas.gamma-1);owner=cut.walls.ownerCell;
for j=1:3
    info=infos{j};r=rates{j};P=info.primitive;u=P(:,2:4);rho=P(:,1);
    temperatureRate=(r(:,5)-sum(u.*r(:,2:4),2)+(.5*sum(u.^2,2)-cv*P(:,5)).*r(:,1))./(rho*cv);
    mapped=mapCutCellWallFields(cut,info,false);data=buildAerothermalViewerData(mesh);state=ones(n,1);
    fluxName='Davis波界HLLC内部对流';if j>1,fluxName='精确理想气体Godunov内部对流，xi=0波系采样';end
    scheme=[fluxName '；Newton应力/Fourier导热；内部黏性功采用算术平均速度代理'];
    scope='M5/30km/20°；常gamma/冻结mu层流研究，未稳态/网格外域或实验合格；Tw给定300K；新增生产合格面积0';
    data.fields=field('pressure','研究墙压','kPa',mapped.pressure.value/1000,state,['当前CV内能EOS恢复墙压；' scheme],scope);
    data.fields(end+1)=field('heatFlux','研究入壁热流（未收敛）','kW/m²',mapped.heatFluxIntoWall.value/1000,state,['真实Dirichlet墙梯度+Fourier；' scheme],scope);
    data.fields(end+1)=field('wallAdjacentTemperature','邻壁单元气体温度','K',P(owner,5),state,['当前CV内能恢复，不是Te/Tr/Tw；' scheme],scope);
    data.fields(end+1)=field('temperatureRate','邻壁气体温度变化率（算子诊断）','K/s',temperatureRate(owner),state, ...
        ['完整五量残差转换的气体温度斜率，不是轨迹热史；' scheme],scope);
    data.fields(end+1)=field('wallInput','给定壁温','K',mapped.wallTemperatureInput.value,2*state,'显式静止无滑移墙温输入，非求解壁温',scope);
    data.fields(end+1)=field('freestreamNormalProjection','给定来流沿实体外法向投影','m/s',-cut.walls.normal*fs.flowVelocityVector(:),state, ...
        '给定U∞点乘原实体外法向，仅几何诊断，不是当地速度或压力',scope);
    data.fields(end).computedStateLabel='几何诊断';
    for z={{'edgeTemperature','外缘温度','未执行合格的当地边界层外缘提取'}, ...
           {'recoveryTemperature','恢复温度','没有适用的当地外缘/恢复温度解'}, ...
           {'wallSolved','求解壁温','Tw300K为输入，未执行独立壁体能量求解'}}
        spec=z{1};data.fields(end+1)=field(spec{1},spec{2},'K',nan(n,1),zeros(n,1),['未求解：' spec{3}],spec{3}); %#ok<AGROW>
    end
    for key=fieldnames(common)',k=key{1};f=data.fields(strcmp({data.fields.key},k));common.(k)=[min(common.(k)(1),min(f.values)),max(common.(k)(2),max(f.values))];end
    summary=summaries{j};caption=sprintf('M5 / 30km / α20° / Tw给定300K；%s；本方程最大残差%.7g s^-1；未收敛；物理时间推进0',labels{j},summary.maximumResidual);
    data.defaultField='heatFlux';data.caseNote=caption;data.coordinateNote='SI / m · X鼻尖正向 · 零攻角来流-X · 灰色为当前场量未知';data.flowDirection=fs.unitFlowDirection;
    data.region=repmat({'真实原三角面研究域'},n,1);data.regions=struct([]);data.heatingGaps=struct([]);data.missingReason=repmat({''},n,1);
    data.resultCondition=struct('caption',caption,'Mach',5,'altitudeM',30000,'angleOfAttackDeg',20, ...
        'flowVelocityVector',fs.flowVelocityVector,'wallTemperatureInputK',300,'physicalTimeAdvanced',0, ...
        'BETimeStepsClaimed',0,'steadyConverged',false);
    data.researchIdentity=struct('caseID',ids{j},'interiorEulerFlux',policies{j},'viscousWorkVelocity','arithmetic', ...
        'discretizationChangedFromHLLCSource',j>1,'HLLCHistoryAppended',false,'productionQualified',false, ...
        'physicalSolution',false,'steadyConverged',false,'newProductionSolvedAreaM2',0,'fullThreeDimensionalMomentum',true);
    datasets{j}=data;
    entry=struct('id',ids{j},'label',labels{j},'interiorEulerFlux',policies{j}, ...
        'resultsURL',['cases/' ids{j} '/index.html'],'resultsJSONURL',['cases/' ids{j} '/results.json'], ...
        'summary',summary,'history',histories{j},'signedHeatPowerW',mapped.summedKnownHeatPower, ...
        'adjacentTemperatureRangeK',[min(P(owner,5)),max(P(owner,5))], ...
        'physicalTimeAdvanced',0,'steadyConverged',false,'newProductionSolvedAreaM2',0);
    if isempty(entries),entries=entry;else,entries(end+1)=entry;end %#ok<AGROW>
end
for j=1:3
    data=datasets{j};data.comparisonColorRanges=common;out=fullfile(folder,'cases',ids{j});mkdir(out);
    for name={'index.html','viewer.css','viewer.js','section.js','export.js'},copyfile(fullfile(root,'web',name{1}),fullfile(out,name{1}));end
    json=jsonencode(data,'ConvertInfAndNaN',true);write(fullfile(out,'results.json'),json);write(fullfile(out,'case-data.js'),['window.ROCKET_CASE=' json ';' newline]);
end
payload=struct('kind','actual-20deg-explicit-exact-Godunov-NS-study-v1','cases',entries,'commonRanges',common, ...
    'sourceFaces',n,'cellCount',numel(cut.cells.volume),'layerCount',cut.quality.layerCount,'gas',gas, ...
    'geometryURL','geometry/index.html','domainURL','domain/index.html','integrationReport',integrationReport, ...
    'readOnlyFrozenStudy',true,'productionQualified',false,'newProductionSolvedAreaM2',0,'browserVisualOrWebGLRuntimeVerified',false);
json=jsonencode(payload,'ConvertInfAndNaN',true);write(fullfile(folder,'exact-study-data.json'),json);write(fullfile(folder,'exact-study-data.js'),['window.ROCKET_EXACT_STUDY=' json ';' newline]);
for name={'rocket-exact.html','rocket-exact.js','rocket-sensitivity.css'}
    dest=name{1};if endsWith(dest,'.html'),dest='index.html';end;copyfile(fullfile(root,'web',name{1}),fullfile(folder,dest));
end
report=struct('status','exported-three-actual-20deg-original-facet-explicit-exact-flux-cases', ...
    'cases',3,'htmlFile',fullfile(folder,'index.html'),'sourceFaces',n,'productionQualified',false, ...
    'browserVisualOrWebGLRuntimeVerified',false,'newProductionSolvedAreaM2',0);
write(fullfile(folder,'study_summary.json'),jsonencode(report));
end
function out=pick(info)
keys={'primitive','wallPressure','wallHeatFluxIntoWall','wallTemperatureInput','wallActiveInFluidDomain'};
out=struct;for name=keys,key=name{1};out.(key)=info.(key);end
end
function f=field(key,label,unit,values,status,method,scope)
n=numel(values);reasons=repmat({scope},n,1);reasons(status>0)={''};
f=struct('key',key,'label',label,'unit',unit,'values',values,'status',status,'methodNames',{{method}}, ...
    'methodIndex',zeros(n,1),'scope',scope,'scopeNames',{{scope}},'scopeIndex',zeros(n,1), ...
    'missingReason',{reasons},'stateNote',scope,'computedStateLabel','未收敛研究数值');
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
function path=publishedInput(root,version,relative)
% 模块四：发布包保留内部已测检查点；无本机outputs时直接读取validation中的原字节。
path=fullfile(root,'outputs',sprintf('v%d_validation',version),relative);
if ~exist(path,'file'),path=fullfile(root,'validation',sprintf('v%d',version),relative);end
assert(exist(path,'file')==2,'ExactFluxView:MissingEvidence','缺少本地内部研究输入：%s',path);
end
