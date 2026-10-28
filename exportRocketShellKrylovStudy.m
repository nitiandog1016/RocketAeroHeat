function report=exportRocketShellKrylovStudy(folder,evidenceRoot)
%EXPORTROCKETSHELLKRYLOVSTUDY 将同态功策略与三种实际迭代结果接回统一科研工作台。
% 模块一：原V/F、同一内部初猜、五量实际残差及不同迭代身份明确保留。
% 模块二：压力/热流/邻壁气体温度/算子温度斜率独立展示，未知Te/Tr/求解Tw不补值。
% 模块三：共同色标、原面探针/任意截线/CSV；研究页不冒充正式稳态环境。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<2,evidenceRoot=fullfile(root,'outputs','v81_validation');end
assert(~exist(folder,'dir'),'ShellKrylovView:Overwrite','请使用新的研究导出目录。');
base=fullfile(root,'outputs','v80_validation');if ~exist(fullfile(base,'advance01','advance_evidence.mat'),'file'),base=fullfile(root,'validation','v80');end
d=load(fullfile(base,'domain02','shell_domain_evidence.mat'),'domain','fs');cut=d.domain;fs=d.fs;mesh=cut.inputGeometry;n=size(mesh.F,1);
s=load(fullfile(base,'advance01','advance_evidence.mat'),'result','gas');gas=s.gas;scale=s.result.residualScale;
infos={pick(s.result.finalSpatialInfo)};rates={s.result.residual};labels={'同一初态 · 旧重构功速度'};methods={'initial-legacy'};
summaries={struct('acceptedNonlinearIterations',0,'candidateTrials',0,'residualEvaluations',0,'elapsedSeconds',0, ...
    'maximumResidual',s.result.residualNorm,'rmsMerit',sqrt(mean((s.result.residual./s.result.residualScale).^2,'all')), ...
    'steadyConverged',false,'physicalTimeAdvanced',0,'BETimeStepsClaimed',0)};histories={struct([])};clear s;
w=load(fullfile(evidenceRoot,'work01','work_velocity_evidence.mat'),'info','rate','report');
infos{2}=pick(w.info);rates{2}=w.rate;labels{2}='同一初态 · 算术功速度';methods{2}='initial-arithmetic';
summaries{2}=summaries{1};summaries{2}.maximumResidual=max(abs(w.rate)./scale,[],'all');
summaries{2}.rmsMerit=sqrt(mean((w.rate./scale).^2,'all'));histories{2}=struct([]);workReport=w.report;clear w;
arithmeticFolder='arithmetic01';if exist(fullfile(evidenceRoot,'arithmetic02','steady_krylov_evidence.mat'),'file'),arithmeticFolder='arithmetic02';end
for spec={{'actual01','steady_line_evidence.mat','旧功 · 单独法向线','legacy-line'}, ...
          {'krylov01','steady_krylov_evidence.mat','旧功 · 完整Krylov','legacy-krylov'}, ...
          {arithmeticFolder,'steady_krylov_evidence.mat','算术功 · 完整Krylov','arithmetic-krylov'}}
    z=spec{1};path=fullfile(evidenceRoot,z{1},z{2});s=load(path,'result','report');
    if isfield(s.report,'totalAcceptedNonlinearIterations'),extra=load(path,'displayHistory');s.displayHistory=extra.displayHistory;clear extra;end
    summary=s.report;history=s.result.history;accepted=s.result.acceptedNonlinearIterations;
    if isfield(s,'displayHistory'),history=s.displayHistory;accepted=s.report.totalAcceptedNonlinearIterations;summary.totalResidualEvaluations=history.residualEvaluations(end);end
    infos{end+1}=pick(s.result.info);rates{end+1}=s.result.rate;summaries{end+1}=summary;histories{end+1}=history;
    labels{end+1}=sprintf('%s · 累计%d次修正',z{3},accepted);methods{end+1}=z{4};clear s; %#ok<AGROW>
end
assert(numel(infos)==5&&n==19752);mkdir(folder);
exportRocketTestGeometryView({mesh},fullfile(folder,'geometry'),{'真实源箭体原19752面，X鼻尖正向'});
exportRocketTestGeometryView({mesh,struct('V',cut.nodeLayers(:,:,2),'F',mesh.F),struct('V',cut.nodeLayers(:,:,end),'F',mesh.F)}, ...
    fullfile(folder,'domain'),{'原墙','第一近壁节点层','最外计算边界；不是激波'});
common=struct('pressure',[Inf -Inf],'heatFlux',[Inf -Inf],'wallAdjacentTemperature',[Inf -Inf],'temperatureRate',[Inf -Inf]);
dataSet=cell(5,1);entries=struct([]);cv=gas.R/(gas.gamma-1);owner=cut.walls.ownerCell;
for j=1:5
    info=infos{j};r=rates{j};P=info.primitive;u=P(:,2:4);rho=P(:,1);
    temperatureRate=(r(:,5)-sum(u.*r(:,2:4),2)+(.5*sum(u.^2,2)-cv*P(:,5)).*r(:,1))./(rho*cv);
    map=mapCutCellWallFields(cut,info,false);data=buildAerothermalViewerData(mesh);state=ones(n,1);
    scope='M5/30km/20°；常gamma/冻结mu层流研究，尚未稳态/网格外域或实验合格；给定Tw300K；新增生产面积0';
    scheme='旧WLS重构内面黏性功速度';if info.arithmeticWork,scheme='显式一阶算术平均内面黏性功速度代理；不是裁剪旧迹';end
    data.fields=field('pressure','研究墙压','kPa',map.pressure.value/1000,state,['当前五量状态EOS恢复原墙压；' scheme],scope);
    data.fields(end+1)=field('heatFlux','研究入壁热流（未收敛）','kW/m²',map.heatFluxIntoWall.value/1000,state,['真实Dirichlet墙梯度+Fourier；流体能量采用' scheme],scope);
    data.fields(end+1)=field('wallAdjacentTemperature','邻壁单元气体温度','K',P(owner,5),state,['当前真实CV内能恢复，不是Te/Tr/Tw；' scheme],scope);
    data.fields(end+1)=field('temperatureRate','邻壁气体温度变化率（算子诊断）','K/s',temperatureRate(owner),state, ...
        ['由完整五量空间残差转换的局部气体温度斜率，不是轨迹热史；' scheme],scope);
    data.fields(end+1)=field('wallInput','给定壁温','K',map.wallTemperatureInput.value,2*state,'显式无滑移等温墙输入，非求解壁温',scope);
    data.fields(end+1)=field('freestreamNormalProjection','给定来流沿实体外法向投影','m/s',-cut.walls.normal*fs.flowVelocityVector(:),state,'给定U∞点乘原实体外法向，仅几何诊断，不是当地速度或压力',scope);
    data.fields(end).computedStateLabel='几何诊断';
    for z={{'edgeTemperature','外缘温度','没有经过验证的边界层外缘提取'}, ...
           {'recoveryTemperature','恢复温度','没有适用的当地外缘/边界层恢复温度解'}, ...
           {'wallSolved','求解壁温','Tw300K是输入，没有执行独立壁体能量求解'}}
        spec=z{1};data.fields(end+1)=field(spec{1},spec{2},'K',nan(n,1),zeros(n,1),['未求解：' spec{3}],spec{3}); %#ok<AGROW>
    end
    for key=fieldnames(common)',k=key{1};f=data.fields(strcmp({data.fields.key},k));common.(k)=[min(common.(k)(1),min(f.values)),max(common.(k)(2),max(f.values))];end
    summary=summaries{j};caption=sprintf('M5 / 30km / α20° / Tw给定300K；%s；原最大残差%.7g s^-1；未收敛；无物理时间推进',labels{j},summary.maximumResidual);
    data.defaultField='heatFlux';data.caseNote=caption;data.coordinateNote='SI / m · X鼻尖正向 · 零攻角来流-X · 灰色为当前场量未知';data.flowDirection=fs.unitFlowDirection;
    data.region=repmat({'原三角面研究域'},n,1);data.regions=struct([]);data.heatingGaps=struct([]);data.missingReason=repmat({''},n,1);
    data.resultCondition=struct('caption',caption,'Mach',5,'altitudeM',30000,'angleOfAttackDeg',20, ...
        'flowVelocityVector',fs.flowVelocityVector,'wallTemperatureInputK',300,'physicalTimeAdvanced',0, ...
        'BETimeStepsClaimed',0,'steadyConverged',false);
    data.researchIdentity=struct('caseID',methods{j},'arithmeticWork',info.arithmeticWork,'productionQualified',false, ...
        'physicalSolution',false,'steadyConverged',false,'newProductionSolvedAreaM2',0,'fullThreeDimensionalMomentum',true);
    dataSet{j}=data;
    entry=struct('id',methods{j},'label',labels{j},'arithmeticWork',info.arithmeticWork, ...
        'resultsURL',['cases/' methods{j} '/index.html'],'resultsJSONURL',['cases/' methods{j} '/results.json'], ...
        'summary',summary,'history',histories{j},'signedHeatPowerW',map.summedKnownHeatPower, ...
        'adjacentTemperatureRangeK',[min(P(owner,5)),max(P(owner,5))],'physicalTimeAdvanced',0, ...
        'steadyConverged',false,'newProductionSolvedAreaM2',0);
    if isempty(entries),entries=entry;else,entries(end+1)=entry;end %#ok<AGROW>
end
for j=1:5
    data=dataSet{j};data.comparisonColorRanges=common;out=fullfile(folder,'cases',methods{j});mkdir(out);
    for name={'index.html','viewer.css','viewer.js','section.js','export.js'},copyfile(fullfile(root,'web',name{1}),fullfile(out,name{1}));end
    json=jsonencode(data,'ConvertInfAndNaN',true);write(fullfile(out,'results.json'),json);write(fullfile(out,'case-data.js'),['window.ROCKET_CASE=' json ';' newline]);
end
payload=struct('kind','rocket-real-20deg-viscous-work-and-full-Krylov-study-v1','cases',entries, ...
    'commonRanges',common,'sourceFaces',n,'cellCount',numel(cut.cells.volume),'layerCount',cut.quality.layerCount, ...
    'workReport',workReport,'gas',gas,'geometryURL','geometry/index.html','domainURL','domain/index.html', ...
    'readOnlyFrozenStudy',true,'productionQualified',false,'newProductionSolvedAreaM2',0,'browserVisualOrWebGLRuntimeVerified',false);
json=jsonencode(payload,'ConvertInfAndNaN',true);write(fullfile(folder,'krylov-study-data.json'),json);write(fullfile(folder,'krylov-study-data.js'),['window.ROCKET_KRYLOV_STUDY=' json ';' newline]);
for name={'rocket-krylov.html','rocket-krylov.js','rocket-sensitivity.css'}
    dest=name{1};if endsWith(dest,'.html'),dest='index.html';end;copyfile(fullfile(root,'web',name{1}),fullfile(folder,dest));
end
report=struct('status','exported-five-real-original-facet-work-and-Krylov-research-cases','cases',5, ...
    'htmlFile',fullfile(folder,'index.html'),'sourceFaces',n,'productionQualified',false, ...
    'browserVisualOrWebGLRuntimeVerified',false,'newProductionSolvedAreaM2',0);write(fullfile(folder,'study_summary.json'),jsonencode(report));
end
function out=pick(info)
keys={'primitive','wallPressure','wallHeatFluxIntoWall','wallTemperatureInput','wallActiveInFluidDomain'};
out=struct;for name=keys,key=name{1};out.(key)=info.(key);end
out.arithmeticWork=isfield(info,'viscousWorkAudit')&&strcmp(info.viscousWorkAudit.policy,'arithmetic');
end
function f=field(key,label,unit,values,status,method,scope)
n=numel(values);reasons=repmat({scope},n,1);reasons(status>0)={''};
f=struct('key',key,'label',label,'unit',unit,'values',values,'status',status,'methodNames',{{method}}, ...
    'methodIndex',zeros(n,1),'scope',scope,'scopeNames',{{scope}},'scopeIndex',zeros(n,1), ...
    'missingReason',{reasons},'stateNote',scope,'computedStateLabel','未收敛研究数值');
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
