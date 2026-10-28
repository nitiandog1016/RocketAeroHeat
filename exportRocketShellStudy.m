function report=exportRocketShellStudy(folder,evidenceRoot)
%EXPORTROCKETSHELLSTUDY 将真实20度三维算子/三步启动接回同窗科研工作台。
% 模块一：只读同一原三角域及已测试证据，保留初猜、给定壁温、未收敛身份。
% 模块二：原19752面映射、共同色标、任意X截线和逐场量来源，不新增生产覆盖。
% 模块三：独立展示原墙、近壁层和外边界真实节点；不画虚构流场或激波。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<2,evidenceRoot=fullfile(root,'outputs','v80_validation');end
assert(~exist(folder,'dir'),'RocketShellView:Overwrite','请使用新的研究导出目录。');
d=load(fullfile(evidenceRoot,'domain02','shell_domain_evidence.mat'),'domain','fs','report');
old=load(fullfile(evidenceRoot,'initial01','initial_residual_evidence.mat'),'info','report');
t=load(fullfile(evidenceRoot,'thermal02','wall_thermal_evidence.mat'),'infoNew','report','gas');
a=load(fullfile(evidenceRoot,'advance01','advance_evidence.mat'),'result','report');
cut=d.domain;fs=d.fs;mesh=cut.inputGeometry;n=size(mesh.F,1);
assert(n==19752&&fs.alpha==20&&a.result.acceptedSteps==3&&~a.result.steadyConverged&& ...
    t.report.newInitialHeatWrongDirectionFaces==0,'RocketShellView:Evidence','证据或数值身份不符。');
mkdir(folder);
exportRocketTestGeometryView({mesh},fullfile(folder,'geometry'),{'原19752三角面；未修改既有表面重网格'});
meshes=cell(3,1);layers=[1 2 size(cut.nodeLayers,3)];
for j=1:3,meshes{j}=struct('V',cut.nodeLayers(:,:,layers(j)),'F',mesh.F);end
exportRocketTestGeometryView(meshes,fullfile(folder,'domain'), ...
    {'原实体墙面','第一近壁节点层：真实节点，不是流场','最外边界节点层：真实节点，不是激波'});
infos={old.info,t.infoNew,a.result.finalSpatialInfo};labels={'旧梯度初态 · 保留失败对照','Dirichlet修正初态','Dirichlet三步数值启动'};
identifiers={'legacy-initial','dirichlet-initial','dirichlet-startup'};entries=struct([]);
common=struct('pressure',[Inf -Inf],'heatFlux',[Inf -Inf],'wallAdjacentTemperature',[Inf -Inf]);
for j=1:3
    info=infos{j};common.pressure=[min(common.pressure(1),min(info.wallPressure)/1000),max(common.pressure(2),max(info.wallPressure)/1000)];
    common.heatFlux=[min(common.heatFlux(1),min(info.wallHeatFluxIntoWall)/1000),max(common.heatFlux(2),max(info.wallHeatFluxIntoWall)/1000)];
    T=info.primitive(cut.walls.ownerCell,5);common.wallAdjacentTemperature=[min(common.wallAdjacentTemperature(1),min(T)),max(common.wallAdjacentTemperature(2),max(T))];
end
for j=1:3
    info=infos{j};mapped=mapCutCellWallFields(cut,info,false);data=buildAerothermalViewerData(mesh);
    scope=['真实三维常gamma/冻结来流mu/Pr0.72层流研究；无湍流/真实空气/网格外域精度资格；' ...
        '当前状态未收敛，不代表飞行时刻的气动热环境，不计生产合格面积'];
    state=ones(n,1);pressureState=state;if j<3,pressureState=2*state;end
    pressureMethod='原三维有限体积当前密度与内能恢复压力；未收敛';
    if j<3,pressureMethod='给定均匀来流初猜通过EOS恢复压力；不是当地飞行壁压解';end
    gradientMethod='显式Dirichlet墙温两切向约束 + 真实正法向距离 + Fourier热流';
    if j==1,gradientMethod='旧完整位移梯度修正 + Fourier；冷初态184面热向错误，保留失败对照';end
    data.fields=field('pressure','研究墙压','kPa',mapped.pressure.value/1000,pressureState,pressureMethod,scope);
    data.fields(end+1)=field('heatFlux','研究入壁热流（未收敛）','kW/m²',mapped.heatFluxIntoWall.value/1000,state,gradientMethod,scope);
    data.fields(end+1)=field('wallAdjacentTemperature','邻壁单元气体温度','K',info.primitive(cut.walls.ownerCell,5),pressureState, ...
        '真实第一层CV当前内能恢复气体温度；不是外缘、恢复温度或壁温',scope);
    data.fields(end+1)=field('wallInput','给定壁温','K',mapped.wallTemperatureInput.value,2*state,'显式等温无滑移壁边界输入',scope);
    projected=-cut.walls.normal*fs.flowVelocityVector(:);
    data.fields(end+1)=field('freestreamNormalProjection','给定来流沿实体外法向投影','m/s',projected,state, ...
        '给定完整U∞点乘原三角面实体外法向；纯几何诊断，不是当地速度或压力',scope);
    data.fields(end).computedStateLabel='几何诊断';
    for spec={{'edgeTemperature','外缘温度','未提取且未验证边界层外缘'}, ...
              {'recoveryTemperature','恢复温度','没有受支持的当地外缘和边界层解'}, ...
              {'wallSolved','求解壁温','本研究仅给定300K壁温，没有执行壁体能量平衡求解'}}
        s=spec{1};data.fields(end+1)=field(s{1},s{2},'K',nan(n,1),zeros(n,1),['未求解：' s{3}],s{3}); %#ok<AGROW>
    end
    time=0;if j==3,time=a.result.time;end
    caption=sprintf(['M5 / 30km / α20° / Tw给定300K；%s；数值启动积分t=%.9g s，' ...
        '不是飞行时间；未收敛；新增合格热流面积0'],labels{j},time);
    data.caseNote=caption;data.defaultField='heatFlux';data.coordinateNote='SI / m · X鼻尖正向 · 零攻角U∞=-X · 灰色为该场量未知';
    data.flowDirection=fs.unitFlowDirection;data.comparisonColorRanges=common;
    data.region=repmat({'原三角面研究域'},n,1);data.regions=struct([]);data.heatingGaps=struct([]);
    data.missingReason=repmat({''},n,1);
    data.resultCondition=struct('caption',caption,'Mach',5,'altitudeM',30000,'angleOfAttackDeg',20, ...
        'flowVelocityVector',fs.flowVelocityVector,'wallTemperatureInputK',300,'startupIntegrationTimeSeconds',time, ...
        'startupTimeIsFlightTime',false,'steadyConverged',false);
    data.researchIdentity=struct('caseID',identifiers{j},'productionQualified',false,'physicalSolution',false, ...
        'newProductionSolvedAreaM2',0,'wallTemperatureKind','prescribed','steadyConverged',false, ...
        'legacyWallGradientFailureComparison',j==1,'fullThreeDimensionalMomentum',true);
    out=fullfile(folder,'cases',identifiers{j});mkdir(out);
    for name={'index.html','viewer.css','viewer.js','section.js','export.js'},copyfile(fullfile(root,'web',name{1}),fullfile(out,name{1}));end
    json=jsonencode(data,'ConvertInfAndNaN',true);write(fullfile(out,'results.json'),json);write(fullfile(out,'case-data.js'),['window.ROCKET_CASE=' json ';' newline]);
    entry=struct('id',identifiers{j},'label',labels{j},'resultsURL',['cases/' identifiers{j} '/index.html'], ...
        'resultsJSONURL',['cases/' identifiers{j} '/results.json'],'startupIntegrationTimeSeconds',time, ...
        'pressureRangePa',[min(mapped.pressure.value),max(mapped.pressure.value)], ...
        'heatRangeWm2',[min(mapped.heatFluxIntoWall.value),max(mapped.heatFluxIntoWall.value)], ...
        'signedHeatPowerW',mapped.summedKnownHeatPower,'steadyConverged',false,'newProductionSolvedAreaM2',0);
    if isempty(entries),entries=entry;else,entries(end+1)=entry;end %#ok<AGROW>
end
history=a.result.history;rows=struct([]);
for j=1:numel(history.time)
    row=struct('acceptedStep',j-1,'startupIntegrationTimeSeconds',history.time(j),'dt',history.acceptedDt(j), ...
        'retries',history.retries(j),'residualPerSecond',history.residualNorm(j),'globalBalanceRelative',history.conservationError(j));
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
payload=struct('kind','rocket-real-20deg-three-dimensional-startup-v1','cases',entries,'history',rows, ...
    'sourceFaces',n,'cellCount',numel(cut.cells.volume),'layerCount',size(cut.nodeLayers,3)-1, ...
    'wallNormalDistanceRangeM',d.report.wallNormalDistanceRangeM,'domainQuality',d.report, ...
    'thermalReport',t.report,'advanceReport',a.report,'gas',t.gas,'commonRanges',common, ...
    'geometryURL','geometry/index.html','domainURL','domain/index.html','readOnlyFrozenStudy',true, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'browserVisualOrWebGLRuntimeVerified',false);
json=jsonencode(payload,'ConvertInfAndNaN',true);write(fullfile(folder,'shell-study-data.json'),json);write(fullfile(folder,'shell-study-data.js'),['window.ROCKET_SHELL_STUDY=' json ';' newline]);
for name={'rocket-shell.html','rocket-shell.js','rocket-sensitivity.css'}
    dest=name{1};if endsWith(dest,'.html'),dest='index.html';end;copyfile(fullfile(root,'web',name{1}),fullfile(folder,dest));
end
report=struct('status','exported-real-20deg-original-facet-three-way-numerical-startup-study', ...
    'htmlFile',fullfile(folder,'index.html'),'cases',3,'sourceFaces',n,'steadyConverged',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'browserVisualOrWebGLRuntimeVerified',false);
write(fullfile(folder,'study_summary.json'),jsonencode(report));
end
function f=field(key,label,unit,values,status,method,scope)
% 模块用途：逐变量记录代数/数值状态，研究数值明确不标成生产工程计算。
n=numel(values);reasons=repmat({scope},n,1);reasons(status>0)={''};
f=struct('key',key,'label',label,'unit',unit,'values',values,'status',status, ...
    'methodNames',{{method}},'methodIndex',zeros(n,1),'scope',scope,'scopeNames',{{scope}}, ...
    'scopeIndex',zeros(n,1),'missingReason',{reasons},'stateNote',scope,'computedStateLabel','未收敛研究数值');
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
