function report=testShockFittedBoundaryFlux(folder,legacyResumeFile)
%TESTSHOCKFITTEDBOUNDARYFLUX 用已知激波验证实际移动面采样及候选通量续算。
% 模块一：已知平面压缩激波，分别在激波上、下游和上游移动采样面。
% 模块二：直接Euler通量和RH恒等式作为解析依据，同时检查刚体/Galilean变换。
% 模块三：旧六步结果逐位回放；两种运动配合新通量均检查整算/续算一致。
% 模块四：候选通量不能替代原压力相容差；NASA7、非压缩分支和预算越界拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupShockFittedFlow;
assert(~exist(folder,'dir'),'ShockFit:Overwrite','边界通量验证必须使用新目录。');mkdir(folder);
gas=createSurfaceFlowGas('caloricallyPerfect',287,1.4);
baseUp=struct('pressure',2000,'temperature',200,'density',2000/(287*200), ...
    'velocityVector',[-2.5*sqrt(1.4*287*200);0;0],'method','解析M2.5平面驻定压缩激波上游');
baseJump=movingShockJump3D(baseUp,gas,[-1 0 0],0);
box=buildShockFittedBoxMesh([2 2 1],[.02 .012 .01]);mesh=box.mesh;
angle=.37;R=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
frames={struct('rotation',eye(3),'boost',[0;0;0]), ...
    struct('rotation',R,'boost',[0;0;0]),struct('rotation',R,'boost',[27;-18;33])};
speeds=[0 -40 40];cases=cell(3,3);dt=1e-8;
U=norm(baseUp.velocityVector);qscale=baseUp.density*[1 repmat(U,1,3) U^2];
maximumAnalyticFluxError=0;maximumUniformError=0;maximumFrameError=0;
for frame=1:3
    rotation=frames{frame}.rotation;boost=frames{frame}.boost;
    up=baseUp;up.velocityVector=rotation*baseUp.velocityVector+boost;
    down=shockFittedState3D([baseJump.state.density ...
        baseJump.state.velocity*rotation'+boost' baseJump.state.temperature],gas,'primitive');
    first=shockFittedState3D([up.density up.velocityVector' up.temperature],gas,'primitive');
    X0=box.nodes*rotation'+[.03 -.02 .01];normalOut=[1 0 0]*rotation';
    for k=1:3
        speed=speeds(k);velocity=(rotation*[speed;0;0]+boost)';
        X1=X0+dt*velocity;
        bd=struct('kind',2*ones(numel(mesh.boundaryFaces),1), ...
            'exteriorQ',repmat(down.Q,numel(mesh.boundaryFaces),1),'upstream',up);
        bd.kind(box.leftBoundaryMask)=1;bd.kind(box.rightBoundaryMask)=5;
        Q=repmat(down.Q,mesh.cellCount,1);
        step=advanceShockFittedALE3D(Q,mesh,X0,X1,dt,bd,gas);
        original=bd;original.kind(box.rightBoundaryMask)=3;
        rh=advanceShockFittedALE3D(Q,mesh,X0,X1,dt,original,gas);
        rows=find(box.rightBoundaryMask);faces=mesh.boundaryFaces(rows);
        normals=step.geometry.normal(faces,:);D=step.geometry.gridNormalSpeedMS(faces);
        Dstar=repmat(normalOut*boost,numel(rows),1);
        % 已知波系的下游/上游均为常态；不调用待测Riemann求解器构造期望通量。
        if speed<=0,sampled=down;else,sampled=first;end
        expected=eulerALE(sampled.Q,sampled.pressure,normals,D);
        upFlux=eulerALE(first.Q,first.pressure,normals,D);
        downFlux=eulerALE(down.Q,down.pressure,normals,D);
        physicalScale=baseUp.density*[U repmat(U^2,1,3) U^3];
        fluxError=max(abs(step.faceALEFlux(faces,:)-expected)./physicalScale,[],'all');
        assert(fluxError<2e-11,'ShockFit:BoundaryAnalyticFlux','实际移动面通量未匹配解析采样态。');
        knownDefect=(Dstar-D).*(down.Q-first.Q);
        identityError=max(abs(downFlux-upFlux-knownDefect)./physicalScale,[],'all');
        assert(identityError<2e-11,'ShockFit:BoundaryRHIdentity','已知移动RH通量差恒等式不闭合。');
        actualScale=up.density*[max(norm(up.velocityVector),first.soundSpeed) ...
            repmat(max(norm(up.velocityVector),first.soundSpeed)^2,1,3) ...
            max(norm(up.velocityVector),first.soundSpeed)^3];
        expectedDefect=(expected-upFlux)./actualScale;
        assert(max(abs(step.shockBoundaryFluxRelativeMismatch(rows,:)-expectedDefect),[],'all')<2e-11);
        assert(max(abs(step.shockRiemannStarPressurePa(rows)/down.pressure-1))<2e-11);
        stateError=max(abs(step.Q-Q)./qscale,[],'all');
        oldStateError=max(abs(rh.Q-Q)./qscale,[],'all');
        if speed<=0
            assert(stateError<2e-11,'ShockFit:BoundaryUniform','位于均匀下游的移动盒域没有保持常态。');
            maximumUniformError=max(maximumUniformError,stateError);
        end
        if speed<0
            assert(oldStateError>1e-5,'ShockFit:BoundaryOracle','负例未分辨强制上游RH通量造成的更新。');
        elseif speed>0
            % 此处RH通量缺陷可为零，但原单元压力不相容；两指标不能相互代替。
            assert(max(abs(step.shockBoundaryFluxRelativeMismatch(rows,:)),[],'all')<2e-11&& ...
                max(abs(step.shockInteriorPressureRelativeMismatch(rows)))>1e-3);
        end
        nodeVelocity=repmat(velocity,size(X0,1),1);
        estimate=shockFittedTimeStep3D(Q,mesh,X0,nodeVelocity,gas,.18,bd);
        oldEstimate=shockFittedTimeStep3D(Q,mesh,X0,nodeVelocity,gas,.18,original);
        assert(isequal(estimate.faceSpectralRadiusMS,oldEstimate.faceSpectralRadiusMS)&& ...
            isequal(estimate.cellSpectralRateM3PerS,oldEstimate.cellSpectralRateM3PerS));
        frameError=0;
        if frame>1
            transformed=transformQ(cases{1,k}.step.Q,rotation,boost);
            frameError=max(abs(step.Q-transformed)./qscale,[],'all');
            assert(frameError<2e-10,'ShockFit:BoundaryFrame','实际推进的旋转/Galilean变换不符。');
        end
        cases{frame,k}=struct('step',step,'originalRHStep',rh,'Q',Q,'upstream',up, ...
            'downstream',down,'boundary',bd,'nodes0',X0,'nodes1',X1,'dt',dt, ...
            'rotation',rotation,'boost',boost,'baseGridSpeedMS',speed,'knownShockSpeedMS',Dstar, ...
            'expectedFlux',expected,'expectedDefect',expectedDefect,'CFL',estimate, ...
            'analyticFluxError',fluxError,'knownRHIdentityError',identityError, ...
            'uniformStateRelativeError',stateError,'originalRHUniformError',oldStateError, ...
            'frameError',frameError);
        maximumAnalyticFluxError=max(maximumAnalyticFluxError,fluxError);
        maximumFrameError=max(maximumFrameError,frameError);
    end
end

%% 模块五：显式拒绝；旧方法仍保持原数值，新增数组须计入预算。
base=cases{1,1};refusals=struct;
hotGas=createSurfaceFlowGas('thermallyPerfect',287,1.4);
refusals.NASA7=reject(@()advanceShockFittedALE3D(base.Q,mesh,base.nodes0,base.nodes1,dt, ...
    base.boundary,hotGas),'ShockFit:FluxGas');
lowState=shockFittedState3D([baseUp.density baseUp.velocityVector' baseUp.temperature/2],gas,'primitive');
bad=base.boundary;bad.exteriorQ=repmat(lowState.Q,numel(mesh.boundaryFaces),1);
refusals.noncompression=reject(@()advanceShockFittedALE3D(repmat(lowState.Q,mesh.cellCount,1),mesh, ...
    base.nodes0,base.nodes1,dt,bad,gas),'ShockFit:ShockBoundaryWave');
refusals.additionalBudget=reject(@()advanceShockFittedALE3D(base.Q,mesh,base.nodes0,base.nodes1,dt, ...
    base.boundary,gas,base.originalRHStep.budget.declaredTotalArrayUpperBytes),'ShockFit:MemoryBudget');

%% 模块六：旧六步逐位回放与两种运动配合候选通量的六步/三加三续算。
old=load(legacyResumeFile,'complete');reference=old.complete;
up=reference.upstream;legacyGas=reference.gas;cap=reference.cap;
options=struct('maximumSteps',6,'endTimeS',1e-4);
legacy=solveShockFittedCap3D(up,legacyGas,cap,options);
assert(isequal(legacy.Q,reference.Q)&&isequal(legacy.nodes,reference.nodes)&& ...
    isequal(legacy.history,reference.history),'ShockFit:BoundaryLegacyReplay','原默认六步数值发生变化。');
motions={'face-riemann-L2','nodal-characteristic'};dynamics=cell(2,1);
for j=1:2
    options.motionMethod=motions{j};options.shockFluxMethod='upstream-riemann';
    complete=solveShockFittedCap3D(up,legacyGas,cap,options);
    partOptions=options;partOptions.maximumSteps=3;
    partOptions.checkpointFile=fullfile(folder,sprintf('method%d_first.mat',j));
    part=solveShockFittedCap3D(up,legacyGas,cap,partOptions);
    next=partOptions;next.resumeFile=partOptions.checkpointFile;
    next.checkpointFile=fullfile(folder,sprintf('method%d_resumed.mat',j));
    resumed=solveShockFittedCap3D(up,legacyGas,cap,next);
    assert(isequal(complete.Q,resumed.Q)&&isequal(complete.nodes,resumed.nodes)&& ...
        isequal(complete.history,resumed.history)&& ...
        isequaln(complete.shockTracePressureMismatchHistory,resumed.shockTracePressureMismatchHistory)&& ...
        isequal(complete.shockBoundaryFluxMismatchHistory,resumed.shockBoundaryFluxMismatchHistory)&& ...
        isequal(complete.conservedIntegrals,resumed.conservedIntegrals)&& ...
        isequal(complete.boundaryFluxRates,resumed.boundaryFluxRates), ...
        'ShockFit:BoundaryResumeReplay','候选通量分段接续不一致。');
    different=next;different.shockFluxMethod='RH';
    refused=reject(@()solveShockFittedCap3D(up,legacyGas,cap,different),'ShockFit:Resume');
    dynamics{j}=struct('complete',complete,'part',part,'resumed',resumed,'fluxSwitchRejection',refused);
end
report=struct('status','tested-actual-moving-boundary-sampling-and-candidate-resume-not-region-qualification', ...
    'analyticCases',9,'maximumAnalyticFluxError',maximumAnalyticFluxError, ...
    'maximumDownstreamUniformError',maximumUniformError,'maximumFrameError',maximumFrameError, ...
    'legacySixStepBitwise',true,'candidateResumeBitwiseBothMotions',true, ...
    'zeroFluxDefectAloneInsufficient',true,'productionQualified',false,'newProductionSolvedAreaM2',0, ...
    'engineProcessID',feature('getpid'));
save(fullfile(folder,'boundary_flux_evidence.mat'),'report','cases','box','gas','baseUp','baseJump', ...
    'frames','speeds','refusals','legacy','dynamics','legacyResumeFile','-v7');
fid=fopen(fullfile(folder,'boundary_flux_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));disp(report);
end

function flux=eulerALE(Q,p,normal,D)
u=Q(2:4)/Q(1);un=normal*u';
flux=[Q(1)*un un.*Q(2:4)+p*normal (Q(5)+p)*un]-D.*Q;
end
function out=transformQ(Q,R,boost)
momentum=Q(:,2:4)*R';
out=[Q(:,1),momentum+Q(:,1).*boost',Q(:,5)+momentum*boost+.5*sum(boost.^2)*Q(:,1)];
end
function id=reject(action,expected)
id='';try,action();catch ex,id=ex.identifier;end
assert(strcmp(id,expected),'ShockFit:RefusalTest','拒绝身份不符：实际%s，需要%s。',id,expected);
end
