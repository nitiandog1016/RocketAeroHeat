function report=testShockFittedTraceReconstruction(folder,capFile)
%TESTSHOCKFITTEDTRACERECONSTRUCTION 校核单侧线性重构、坐标不变性与真实保存端点。
% 模块一：制造的是单元中心原始变量样本，检验插值线性一致性，不冒充PDE制造解。
% 模块二：常值逐位保持，尺度/旋转/Galilean变换、原温域限制与严格拒绝路径。
% 模块三：同一保存球鼻状态比较常值/重构的Riemann速度与共享节点残差，不推进流场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupShockFittedFlow;
assert(~exist(folder,'dir'),'ShockFit:Overwrite','重构测试使用全新目录。');mkdir(folder);
box=buildShockFittedBoxMesh([4 3 3],[.04 .012 .009]);mesh=box.mesh;X=box.nodes;
gas=createSurfaceFlowGas('caloricallyPerfect',287.05287,1.4);faceIDs=mesh.boundaryFaces;
base=[1.1 150 20 -10 1e5];G=[.1 12 -8;2 -1 .4;-3 .5 1.2;8 -7 4;2e4 -3e4 1e4];
W=base+box.cellCentroidM*G';T=W(:,5)./(gas.R*W(:,1));s=shockFittedState3D([W(:,1:4) T],gas,'primitive');
linear=reconstructShockFittedTrace3D(s.Q,mesh,X,gas,faceIDs);
expected=base+box.boundaryCentroidM*G';scale=[1 300 300 300 1e5];
linearError=max(abs(linear.primitiveRP-expected)./scale,[],'all');
assert(linearError<2e-12&&all(linear.slopeFactor==1),'ShockFit:TraceLinear','线性样本未按原三维梯度恢复。');
uniform=shockFittedState3D([1.2 200 30 -20 300],gas,'primitive');uniformQ=repmat(uniform.Q,mesh.cellCount,1);
constant=reconstructShockFittedTrace3D(uniformQ,mesh,X,gas,faceIDs);
assert(isequal(constant.Q,uniformQ(mesh.owner(faceIDs),:)),'ShockFit:TraceConstant','常状态没有逐位保持。');

axis=[1;2;-3];axis=axis/norm(axis);K=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];
R=eye(3)+sin(.37)*K+(1-cos(.37))*K*K;shift=[.37 -.29 .12];lengthScale=2.3;boost=[19 -7 10];
transformedX=lengthScale*X*R'+shift;transformedW=W;transformedW(:,2:4)=W(:,2:4)*R'+boost;
transformedState=shockFittedState3D([transformedW(:,1:4) T],gas,'primitive');
transformed=reconstructShockFittedTrace3D(transformedState.Q,mesh,transformedX,gas,faceIDs);
expectedFrame=expected;expectedFrame(:,2:4)=expected(:,2:4)*R'+boost;
frameError=max(abs(transformed.primitiveRP-expectedFrame)./scale,[],'all');
assert(frameError<2e-12,'ShockFit:TraceFrame','重构不满足声明的尺度/旋转/Galilean一致性。');

hotT=1000+25500*box.cellCentroidM(:,1);
hot=shockFittedState3D([ones(mesh.cellCount,1)*1.1 repmat([300 20 0],mesh.cellCount,1) hotT],gas,'primitive');
limited=reconstructShockFittedTrace3D(hot.Q,mesh,X,gas,faceIDs);
assert(limited.limitedFaceCount>0&&all(limited.slopeFactor>=0&limited.slopeFactor<=1));
limitedState=shockFittedState3D(limited.Q,gas,'conservative');
assert(all(limitedState.temperature>=80&limitedState.temperature<=2000));
refusals=struct;
refusals.unlimitedInvalidTrace=reject(@()reconstructShockFittedTrace3D(hot.Q,mesh,X,gas,faceIDs, ...
    struct('slopeLimiter','none')),'ShockFit:TraceAdmissibility');
refusals.budget=reject(@()reconstructShockFittedTrace3D(s.Q,mesh,X,gas,faceIDs, ...
    struct('maximumArrayBytes',1)),'ShockFit:TraceMemoryBudget');
refusals.interiorFace=reject(@()reconstructShockFittedTrace3D(s.Q,mesh,X,gas,find(mesh.neighbour>0,1)), ...
    'ShockFit:TraceFaces');
nasa=createSurfaceFlowGas('thermallyPerfect',gas.R,1.4);
refusals.NASA7=reject(@()reconstructShockFittedTrace3D(s.Q,mesh,X,nasa,faceIDs),'ShockFit:TraceGas');
small=buildShockFittedBoxMesh([1 1 1],[1 1 1]);smallQ=repmat(uniform.Q,small.mesh.cellCount,1);
refusals.rank=reject(@()reconstructShockFittedTrace3D(smallQ,small.mesh,small.nodes,gas,small.mesh.boundaryFaces, ...
    struct('minimumNeighbors',3,'initialRings',1)),'ShockFit:TraceStencilRank');

saved=load(capFile,'study');study=saved.study;cap=study.cap;shell=cap.shell;realMesh=shell.mesh;
outer=realMesh.boundaryFaces(shell.outerBoundaryMask);geom=shockFittedALEGeometry(realMesh,study.nodes,study.nodes,1);
realTrace=reconstructShockFittedTrace3D(study.Q,realMesh,study.nodes,study.gas,outer);
cellShock=shockFittedBoundarySpeed3D(study.Q(realMesh.owner(outer),:),geom.normal(outer,:),study.upstream,study.gas);
traceShock=shockFittedBoundarySpeed3D(realTrace.Q,geom.normal(outer,:),study.upstream,study.gas);
cellMotion=shockFittedShellMotion3D(shell,study.nodes,cellShock.normalSpeedMS);
traceMotion=shockFittedShellMotion3D(shell,study.nodes,traceShock.normalSpeedMS);
cellJump=movingShockJump3D(study.upstream,study.gas,-geom.normal(outer,:),-cellMotion.faceNormalSpeedMS);
traceJump=movingShockJump3D(study.upstream,study.gas,-geom.normal(outer,:),-traceMotion.faceNormalSpeedMS);
cellState=shockFittedState3D(study.Q(realMesh.owner(outer),:),study.gas,'conservative');
traceState=shockFittedState3D(realTrace.Q,study.gas,'conservative');
probe=struct('sameSavedPhysicalState',true,'flowAdvanced',false,'limitedFaces',realTrace.limitedFaceCount, ...
    'cellTargetMaximumSpeedMS',max(abs(cellShock.normalSpeedMS)), ...
    'traceTargetMaximumSpeedMS',max(abs(traceShock.normalSpeedMS)), ...
    'cellMaximumProjectionResidualMS',max(abs(cellMotion.projectionResidualMS)), ...
    'traceMaximumProjectionResidualMS',max(abs(traceMotion.projectionResidualMS)), ...
    'cellProjectedPressureMismatch',max(abs(cellState.pressure./cellJump.state.pressure-1)), ...
    'traceProjectedPressureMismatch',max(abs(traceState.pressure./traceJump.state.pressure-1)), ...
    'sourceFileSHA256',rocketStreamingFileSHA256(capFile));
report=struct('status','tested-one-sided-interface-reconstruction-not-a-new-flow-solution', ...
    'linearPrimitivePointSampleError',linearError,'transformedSampleError',frameError, ...
    'constantQBitwise',true,'limitedManufacturedFaces',limited.limitedFaceCount,'realProbe',probe, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'trace_reconstruction_evidence.mat'),'report','box','gas','base','G','s','linear','expected', ...
    'R','shift','lengthScale','boost','transformed','expectedFrame','uniformQ','constant', ...
    'hot','limited','refusals','capFile','realTrace','cellShock','traceShock','cellMotion','traceMotion','-v7');
fid=fopen(fullfile(folder,'trace_reconstruction_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));disp(report);
end
function id=reject(action,expected)
id='';try,action();catch ex,id=ex.identifier;end
assert(strcmp(id,expected),'ShockFit:TraceRefusalTest','实际拒绝%s，需要%s。',id,expected);
end
