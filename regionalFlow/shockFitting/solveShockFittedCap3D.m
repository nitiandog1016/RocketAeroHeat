function study=solveShockFittedCap3D(upstream,gas,cap,options)
%SOLVESHOCKFITTEDCAP3D 固定真实球鼻壁面、由下游相容推进共同激波的三维候选区域求解。
% 模块一：从给定初始外壳的静止RH初始化；不是预先给定的最终压力或激波形状。
% 模块二：每步重算下游相容、共享节点速度和ALE声学CFL；侧面只允许全特征出流。
% 模块三：新几何/状态正性失败时只缩短时间步并记录，不裁温/压力或修改边界模型。
% 模块四：保存实际五量、几何、拒绝和未稳态指标；未通过区域验证不授予生产资格。
% 模块五：推进前把调用者保留数组、一步输入及内部工作数组合计预算；不是进程峰值。
defaults=struct('CFL',.18,'maximumSteps',200,'endTimeS',2e-5,'maximumArrayBytes',128*1024^2, ...
    'checkpointFile','','checkpointEvery',25,'resumeFile','');
if nargin<4,options=struct;end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'ShockFit:CapOptions','未知曲面推进选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumSteps,{'double'},{'scalar','real','finite','integer','positive'});
validateattributes(options.endTimeS,{'double'},{'scalar','real','finite','positive'});
validateattributes(options.CFL,{'double'},{'scalar','real','finite','positive','<=',1});
validateattributes(options.checkpointEvery,{'double'},{'scalar','real','finite','integer','positive'});
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:ShockSpeedGas','曲面推进当前只常gamma，NASA7未实现。');
shell=cap.shell;mesh=shell.mesh;N=numel(shell.originalVertexIDs);X=shell.nodes;
b=mesh.boundaryFaces;outerMask=shell.outerBoundaryMask;outerFaces=b(outerMask);sideMask=shell.sideBoundaryMask;
g=shockFittedALEGeometry(mesh,X,X,1);jump=movingShockJump3D(upstream,gas,-g.normal(outerFaces,:),zeros(numel(outerFaces),1));
faceLabels=shell.cellOriginalFaceIDs(mesh.owner(outerFaces));[found,source]=ismember(shell.cellOriginalFaceIDs,faceLabels);
assert(all(found),'ShockFit:InitialColumns','单元没有对应的初始激波入口。');
Q=jump.state.Q(source,:);initialQ=Q;initialNodes=X;time=0;history=zeros(0,16);rejections=cell(0,1);last=struct;
conservedIntegrals=sum(g.volume0.*Q,1);boundaryFluxRates=zeros(0,5);
sourceIdentity=identity();
if ~isempty(options.resumeFile)
    old=load(options.resumeFile,'checkpoint');prior=old.checkpoint;
    assert(isequaln(prior.upstream,upstream)&&isequaln(prior.gas,gas)&&isequaln(prior.cap,cap)&& ...
        isequaln(prior.sourceIdentity,sourceIdentity),'ShockFit:Resume','检查点输入、几何或推进方程版本不同。');
    Q=prior.Q;X=prior.nodes;time=prior.timeS;history=prior.history;rejections=prior.rejections;
    initialQ=prior.initialQ;initialNodes=prior.initialNodes;
    conservedIntegrals=prior.conservedIntegrals;boundaryFluxRates=prior.boundaryFluxRates;
end
boundary=struct('kind',ones(numel(b),1),'exteriorQ',nan(numel(b),5),'upstream',upstream);
boundary.kind(shell.wallBoundaryMask)=2;boundary.kind(outerMask)=3;boundary.kind(sideMask)=4;
freestream=shockFittedState3D([upstream.density upstream.velocityVector' upstream.temperature],gas,'primitive');
U=norm(upstream.velocityVector);scales=[upstream.density,repmat(upstream.density*U,1,3),upstream.density*U^2];
stopReason='maximum-steps-without-steady-qualification';steady=false;timer=tic;
proposed=[];dt=NaN;shock=struct;motion=struct;estimate=struct;
lastInputQ=Q;lastNodes0=X;lastShock=struct;lastMotion=struct;lastEstimate=struct;maximumDeclaredRegionArrayBytes=0;
stepWorkspaceUpper=8*(1000*mesh.faceCount+200*mesh.cellCount+120*mesh.faceCount+24*mesh.nodeCount);
try
for iteration=1:options.maximumSteps
    if time>=options.endTimeS,stopReason='requested-transient-time-reached-not-steady';break;end
    variables=whos;inputs=whos('Q','mesh','X','boundary','gas');
    regionUpper=sum([variables.bytes])+sum([inputs.bytes])+numel(X)*8+stepWorkspaceUpper;
    maximumDeclaredRegionArrayBytes=max(maximumDeclaredRegionArrayBytes,regionUpper);
    assert(regionUpper<=options.maximumArrayBytes,'ShockFit:RegionMemoryBudget', ...
        '保留数组、一步输入与内部临时数组的合计预算超限，不启动新步。');
    geom=shockFittedALEGeometry(mesh,X,X,1);
    shock=shockFittedBoundarySpeed3D(Q(mesh.owner(outerFaces),:),geom.normal(outerFaces,:),upstream,gas);
    motion=shockFittedShellMotion3D(shell,X,shock.normalSpeedMS);
    estimate=shockFittedTimeStep3D(Q,mesh,X,motion.nodeVelocityMS,gas,options.CFL,boundary);
    dt=min(estimate.dt,options.endTimeS-time);accepted=false;
    for reduction=0:12
        proposed=X+dt*motion.nodeVelocityMS;
        assert(isequal(proposed(1:N,:),initialNodes(1:N,:)),'ShockFit:BodyMotion','真实壁面不能随候选激波移动。');
        try
            candidate=advanceShockFittedALE3D(Q,mesh,X,proposed,dt,boundary,gas,options.maximumArrayBytes);
            accepted=true;break;
        catch ex
            allowed={'ShockFit:CellVolume','ShockFit:Temperature','CutExactEuler:Temperature', ...
                'CutExactEuler:State','CutExactEuler:InternalEnergy','ShockPatch:Temperature'};
            if ~ismember(ex.identifier,allowed),rethrow(ex);end
            rejections{end+1,1}=struct('timeS',time,'dt',dt,'reduction',reduction,'identifier',ex.identifier,'message',ex.message); %#ok<AGROW>
            dt=dt/2;
        end
    end
    assert(accepted,'ShockFit:TimeStepRejected','有界时间步缩短后仍无合法几何/状态，拒绝继续。');
    rate=(candidate.Q-Q)/dt;relativeRate=max(abs(rate)./scales,[],'all')*cap.radiusM/U;
    globalChange=sum(candidate.geometry.volume1.*candidate.Q-candidate.geometry.volume0.*Q,1);
    balance=max(abs(globalChange+dt*candidate.boundaryFluxRate)./candidate.conservationScale);
    gap=vecnorm(proposed(shell.layers*N+(1:N),:)-proposed(1:N,:),2,2);
    projection=max(abs(motion.projectionResidualMS));idealSpeed=max(abs(shock.normalSpeedMS));
    exactGrid=candidate.geometry.gridNormalSpeedMS(outerFaces);
    history(end+1,:)=[time+dt,dt,relativeRate,idealSpeed/U,projection/U,max(abs(exactGrid))/U, ...
        min(candidate.state.temperature),max(candidate.state.temperature),min(gap),max(gap), ...
        balance,max(candidate.geometry.gclRelative),min(candidate.supersonicOutflowAudit.relativeNormalMach), ...
        candidate.budget.declaredTotalArrayUpperBytes,reduction,max(abs(candidate.shockInteriorPressureRelativeMismatch(outerMask)))]; %#ok<AGROW>
    lastInputQ=Q;lastNodes0=X;lastShock=shock;lastMotion=motion;lastEstimate=estimate;
    conservedIntegrals(end+1,:)=sum(candidate.geometry.volume1.*candidate.Q,1); %#ok<AGROW>
    boundaryFluxRates(end+1,:)=candidate.boundaryFluxRate; %#ok<AGROW>
    Q=candidate.Q;X=proposed;time=time+dt;last=candidate;
    if mod(iteration,25)==0||iteration==1
        fprintf('FITTED_CAP step=%d t=%.9g dt=%.6g residual=%.6g shockSpeed=%.6g projection=%.6g minOutM=%.6g gap=%.6g--%.6g\n', ...
            size(history,1),time,dt,relativeRate,idealSpeed/U,projection/U,history(end,13),min(gap),max(gap));
    end
    if ~isempty(options.checkpointFile)&&mod(iteration,options.checkpointEvery)==0,saveCheckpoint();end
    if relativeRate<5e-4&&idealSpeed/U<1e-4&&projection/U<1e-4
        steady=true;stopReason='candidate-discrete-steady-gates-met-still-needs-domain-and-mesh-validation';break;
    end
end
catch ex
    if ~isempty(options.checkpointFile)
        saveCheckpoint();
        failure=struct('identifier',ex.identifier,'message',ex.message,'timeS',time,'dt',dt, ...
            'lastAcceptedQ',Q,'lastAcceptedNodes',X,'trialNodes',proposed,'shock',shock,'motion',motion,'CFL',estimate); %#ok<NASGU>
        save([options.checkpointFile '.failure.mat'],'failure','-v7');
    end
    rethrow(ex);
end
if ~isempty(options.checkpointFile),saveCheckpoint();end
study=struct('cap',cap,'upstream',upstream,'gas',gas,'initialQ',initialQ,'initialNodes',initialNodes, ...
    'Q',Q,'nodes',X,'timeS',time,'history',history,'rejections',{rejections},'lastStep',last, ...
    'lastInputQ',lastInputQ,'lastNodes0',lastNodes0,'lastShockCompatibility',lastShock, ...
    'lastMotion',lastMotion,'lastTimeStepEstimate',lastEstimate, ...
    'maximumDeclaredRegionArrayBytes',maximumDeclaredRegionArrayBytes,'processPeakMemoryMeasured',false, ...
    'conservedIntegrals',conservedIntegrals,'boundaryFluxRates',boundaryFluxRates, ...
    'boundary',boundary,'sourceIdentity',sourceIdentity,'candidateDiscreteSteady',steady,'stopReason',stopReason, ...
    'elapsedSeconds',toc(timer),'method','下游精确Riemann右压缩波相容/共享射线节点投影/三维ALE五量推进/超声速出口', ...
    'initialStateIsSolvedFlow',false,'candidatePressureAvailable',~isempty(fieldnames(last)), ...
    'bodyVerticesFixed',isequal(X(1:N,:),initialNodes(1:N,:)), ...
    'shockShapeSolved',false,'heatFluxSolved',false,'wallTemperatureSolved',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'upstreamSoundSpeedMS',freestream.soundSpeed);

    function saveCheckpoint()
        checkpoint=struct('upstream',upstream,'gas',gas,'cap',cap,'sourceIdentity',sourceIdentity, ...
            'Q',Q,'nodes',X,'timeS',time,'history',history,'rejections',{rejections}, ...
            'initialQ',initialQ,'initialNodes',initialNodes); %#ok<NASGU>
        checkpoint.conservedIntegrals=conservedIntegrals;checkpoint.boundaryFluxRates=boundaryFluxRates;
        temporary=[options.checkpointFile '.writing.mat'];save(temporary,'checkpoint','-v7');movefile(temporary,options.checkpointFile,'f');
    end
end
function result=identity()
result=struct;
for key={'solveShockFittedCap3D','shockFittedBoundarySpeed3D','shockFittedShellMotion3D','shockFittedTimeStep3D', ...
        'advanceShockFittedALE3D','shockFittedALEGeometry','shockFittedALEFlux3D','shockFittedState3D', ...
        'shockFittedSlipWallFlux3D','movingShockJump3D','shockLayerNormalJumpBatch','shockLayerThermalBatch', ...
        'surfaceFlowThermodynamicState','cutCellExactEulerFlux'}
    result.(key{1})=rocketStreamingFileSHA256(which(key{1}));
end
end
