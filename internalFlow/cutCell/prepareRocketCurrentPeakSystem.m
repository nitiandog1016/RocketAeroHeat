function preparation=prepareRocketCurrentPeakSystem(folder,stateFile,geometryFile,geometryVariable)
%PREPAREROCKETCURRENTPEAKSYSTEM 从实际接受检查点构造当前方向系统，不执行失败暖启动。
% 模块一：显式读取result末态，完整NS逐值回放；最后接受方向只作当前状态权重参考。
% 模块二：重建当前缓存、完整C、局部五量逆、全部像列与峰值度量；物理网格/NS不变。
% 模块三：保存各独立产物及SHA清单，后续求解从零或同一清单的检查点开始。
folder=canonical(folder);stateFile=canonical(stateFile);geometryFile=canonical(geometryFile);
assert(~exist(folder,'dir')&&isfile(stateFile)&&isfile(geometryFile)&&ischar(geometryVariable), ...
    'PeakSystem:Input','需要新目录、真实状态文件及明确几何变量。');
s=load(stateFile,'result','gas','boundary');r=s.result;
assert(all(isfield(r,{'Q','rate','info','residualScale','trace'}))&&~isempty(r.trace), ...
    'PeakSystem:State','必须提供含真实末态和最后接受方向的result检查点。');
Q=r.Q;scale=r.residualScale;gas=s.gas;boundary=s.boundary;reference=r.trace{end}.linearUpdateU;
assert(isreal(reference)&&isequal(size(reference),size(Q))&&all(isfinite(reference),'all'), ...
    'PeakSystem:Reference','最后接受方向只能作为有限的同尺寸数值参考。');
d=load(geometryFile,geometryVariable);cut=d.(geometryVariable);
[rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);
replayReport=assertRocketNSReplay(rate,info,r.rate,r.info,cut.cells.volume);clear s r d;
mkdir(folder);sourceBinding=struct('sourceFile',stateFile,'sourceFileSHA256',rocketStreamingFileSHA256(stateFile),'stateLocation','result');
geometryBinding=struct('file',geometryFile,'variable',geometryVariable,'sha256',rocketStreamingFileSHA256(geometryFile));
sourceIdentity=rocketSelectedDirectionSourceIdentity;methodIdentity=rocketPeakSystemSourceIdentity;rhs=rate./scale;timer=tic;
cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);
cacheFile=fullfile(folder,'current_peak_cache_evidence.mat');save(cacheFile,'cache','Q','rate','scale','gas','boundary','sourceBinding','sourceIdentity','geometryBinding','-v7');
graph=buildRocketShellGraphPreconditioner(Q,cache.geometry,gas,boundary,scale,1e-7,1);
operator=buildRocketWlsNormalPreconditioner(graph,cache,512*1024^2);clear graph cache info cut;
operatorFile=fullfile(folder,'current_peak_operator_evidence.mat');save(operatorFile,'operator','Q','rhs','sourceBinding','sourceIdentity','-v7');
packed=packRocketSteadyGraph(operator,struct('ordering','symrcm'));clear operator;
packedFile=fullfile(folder,'current_peak_packed_evidence.mat');save(packedFile,'packed','operatorFile','sourceBinding','sourceIdentity','-v7');
partition=buildRocketStrongPatches(packed,8);factor=factorRocketCoupledPatches(packed,partition,struct('failureEvidenceFile',fullfile(folder,'factor_failure.mat')));clear partition;
factorFile=fullfile(folder,'current_peak_factor_evidence.mat');save(factorFile,'factor','sourceBinding','sourceIdentity','-v7');
[columnScale,columnNorm,columnConstruction]=buildRocketParentImageScale(packed,factor);
scaleFile=fullfile(folder,'current_peak_scale_evidence.mat');save(scaleFile,'columnScale','columnNorm','columnConstruction','sourceBinding','sourceIdentity','-v7');
baseMetric=buildRocketPrimitiveUpdateMetric(Q,scale,gas);[metric,weighting]=reweightRocketPrimitiveUpdateMetric(baseMetric,reference);
referenceMetric=applyRocketPrimitiveUpdateMetric(metric,reference);metricNormalization=max(1,norm(referenceMetric,'fro'));
metricFile=fullfile(folder,'current_peak_metric_evidence.mat');save(metricFile,'baseMetric','metric','weighting','reference','rhs','metricNormalization','sourceBinding','sourceIdentity','-v7');
storage=whos('packed','factor','columnScale','columnNorm','baseMetric','metric','reference','referenceMetric');cacheBytes=sum([storage.bytes]);
assert(cacheBytes<=1024*1024^2,'PeakSystem:Budget','声明的方向输入缓存超1024MiB上界。');
artifacts=struct('cache',receipt(cacheFile),'operator',receipt(operatorFile),'packed',receipt(packedFile), ...
    'factor',receipt(factorFile),'scale',receipt(scaleFile),'metric',receipt(metricFile));
preparation=struct('kind','current-state-peak-direction-system-v1','sourceBinding',sourceBinding,'geometryBinding',geometryBinding, ...
    'sourceIdentity',sourceIdentity,'methodIdentity',methodIdentity,'artifacts',artifacts,'cellCount',size(Q,1), ...
    'referencePolicy','last-accepted-direction-for-current-metric-only','referenceReusedAsInitialGuess',false, ...
    'oldDerivativeOrQualificationInherited',false,'currentStateCompleteNSExactlyReplayed',true,'replayReport',replayReport, ...
    'metricNormalization',metricNormalization,'weighting',weighting,'factorReport',factor.factorizationReport, ...
    'columnConstruction',columnConstruction,'actualDirectionInputCacheArrayBytes',cacheBytes,'processPeakMemoryClaimed',false, ...
    'constructionSeconds',toc(timer),'newActualNSEvaluations',1,'newFlowAcceptances',0,'physicalSolution',false, ...
    'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
fid=fopen(fullfile(folder,'peak_system.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(preparation));
fprintf('PEAK_SYSTEM_PREPARED cells=%d NS=1 referenceOnly=true cacheBytes=%.0f\n',size(Q,1),cacheBytes);
end
function file=canonical(file),assert(ischar(file));file=char(java.io.File(file).getCanonicalPath());end
function out=receipt(file),out=struct('file',file,'sha256',rocketStreamingFileSHA256(file));end
