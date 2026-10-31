function report=testShockFittedNodalMotion(folder,capFile)
%TESTSHOCKFITTEDNODALMOTION 验证共享节点重构、特征速度和原逐面误差保留。
% 模块一：三维仿射原始变量点样本、常值及坐标变换校核，不冒充PDE精度证明。
% 模块二：同cp移动RH合成状态恢复已知速度，检验旋转/Galilean和严格拒绝。
% 模块三：在同一真实球鼻保存状态上作节点运动探针，不宣称推进或新增面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupShockFittedFlow;
assert(~exist(folder,'dir'),'ShockFit:Overwrite','节点运动测试使用新目录。');mkdir(folder);
box=buildShockFittedBoxMesh([4 3 3],[.04 .012 .009]);mesh=box.mesh;X=box.nodes;
gas=createSurfaceFlowGas('caloricallyPerfect',287.05287,1.4);nodeIDs=unique(mesh.faces(mesh.boundaryFaces,:));
base=[1.1 150 20 -10 1e5];G=[.1 12 -8;2 -1 .4;-3 .5 1.2;8 -7 4;2e4 -3e4 1e4];
W=base+box.cellCentroidM*G';T=W(:,5)./(gas.R*W(:,1));state=shockFittedState3D([W(:,1:4) T],gas,'primitive');
linear=reconstructShockFittedNodalTrace3D(state.Q,mesh,X,gas,nodeIDs);
expected=base+X(nodeIDs,:)*G';scales=[1 300 300 300 1e5];
linearError=max(abs(linear.primitiveRP-expected)./scales,[],'all');
assert(linearError<2e-12&&all(linear.slopeFactor==1),'ShockFit:NodalLinearTest','共享节点仿射重构未通过。');
uniform=shockFittedState3D([1.2 200 30 -20 300],gas,'primitive');
constant=reconstructShockFittedNodalTrace3D(repmat(uniform.Q,mesh.cellCount,1),mesh,X,gas,nodeIDs);
assert(isequal(constant.Q,repmat(uniform.Q,numel(nodeIDs),1)),'ShockFit:NodalConstantTest','节点常值没有逐位保持。');
axis=[1;2;-3];axis=axis/norm(axis);K=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];
R=eye(3)+sin(.37)*K+(1-cos(.37))*K*K;shift=[.37 -.29 .12];lengthScale=2.3;boost=[19 -7 10];
transformedX=lengthScale*X*R'+shift;transformedW=W;transformedW(:,2:4)=W(:,2:4)*R'+boost;
transformedState=shockFittedState3D([transformedW(:,1:4) T],gas,'primitive');
transformed=reconstructShockFittedNodalTrace3D(transformedState.Q,mesh,transformedX,gas,nodeIDs);
expectedFrame=expected;expectedFrame(:,2:4)=expected(:,2:4)*R'+boost;
frameError=max(abs(transformed.primitiveRP-expectedFrame)./scales,[],'all');
assert(frameError<2e-12,'ShockFit:NodalFrameTest','节点重构坐标变换误差超限。');
hotT=1000+25500*box.cellCentroidM(:,1);
hot=shockFittedState3D([ones(mesh.cellCount,1)*1.1 repmat([300 20 0],mesh.cellCount,1) hotT],gas,'primitive');
limited=reconstructShockFittedNodalTrace3D(hot.Q,mesh,X,gas,nodeIDs);
assert(limited.limitedNodeCount>0);hotTrace=shockFittedState3D(limited.Q,gas,'conservative');
assert(all(hotTrace.temperature>=80&hotTrace.temperature<=2000));

saved=load(capFile,'study');study=saved.study;up=study.upstream;realGas=study.gas;
normals=[1 0 0;.9 .3 .1;.8 -.2 .4;.7 .5 -.3];normals=normals./vecnorm(normals,2,2);
sound=sqrt(realGas.gamma*realGas.R*up.temperature);knownMach=[1.2;1.8;2.5;3.5];knownSpeed=normals*up.velocityVector+sound*knownMach;
jump=movingShockJump3D(up,realGas,-normals,-knownSpeed);
characteristic=shockFittedCharacteristicSpeed3D(jump.state.Q,normals,up,realGas);
speedError=max(abs(characteristic.normalSpeedMS-knownSpeed))/sound;
assert(speedError<2e-11,'ShockFit:CharacteristicRecoveryTest','原移动RH速度没有恢复。');
rotatedUp=up;rotatedUp.velocityVector=R*up.velocityVector+boost';rotatedNormals=normals*R';
rotatedState=shockFittedState3D([jump.state.density jump.state.velocity*R'+boost jump.state.temperature],realGas,'primitive');
rotated=shockFittedCharacteristicSpeed3D(rotatedState.Q,rotatedNormals,rotatedUp,realGas);
rotatedExpected=knownSpeed+rotatedNormals*boost';
characteristicFrameError=max(abs(rotated.normalSpeedMS-rotatedExpected))/sound;
assert(characteristicFrameError<2e-11,'ShockFit:CharacteristicFrameTest','特征速度坐标变换未通过。');
refusals=struct;
refusals.budget=reject(@()reconstructShockFittedNodalTrace3D(state.Q,mesh,X,gas,nodeIDs,1),'ShockFit:NodalTraceBudget');
nasa=createSurfaceFlowGas('thermallyPerfect',gas.R,1.4);
refusals.NASA7=reject(@()shockFittedCharacteristicSpeed3D(jump.state.Q,normals,up,nasa),'ShockFit:CharacteristicGas');
noncompressive=shockFittedState3D([up.density (up.velocityVector-sound*[1;0;0])' up.temperature],realGas,'primitive');
refusals.noncompressive=reject(@()shockFittedCharacteristicSpeed3D(noncompressive.Q,[1 0 0],up,realGas), ...
    'ShockFit:CharacteristicWave');
small=buildShockFittedBoxMesh([1 1 1],[1 1 1]);smallNodes=unique(small.mesh.faces(small.mesh.boundaryFaces,:));
refusals.rank=reject(@()reconstructShockFittedNodalTrace3D(repmat(uniform.Q,small.mesh.cellCount,1), ...
    small.mesh,small.nodes,gas,smallNodes),'ShockFit:NodalTraceRank');

cap=study.cap;shell=cap.shell;realMesh=shell.mesh;
nodalMotion=shockFittedNodalMotion3D(shell,study.nodes,study.Q,up,realGas);
outer=realMesh.boundaryFaces(shell.outerBoundaryMask);geo=shockFittedALEGeometry(realMesh,study.nodes,study.nodes,1);
N=numel(shell.originalVertexIDs);dt=1e-8;next=study.nodes+dt*nodalMotion.nodeVelocityMS;
moving=shockFittedALEGeometry(realMesh,study.nodes,next,dt);
assert(isequal(next(1:N,:),study.nodes(1:N,:))&&max(moving.gclRelative)<2e-11);
actualJump=movingShockJump3D(up,realGas,-geo.normal(outer,:),-nodalMotion.faceNormalSpeedMS);
cellState=shockFittedState3D(study.Q(realMesh.owner(outer),:),realGas,'conservative');
probe=struct('sameSavedPhysicalState',true,'flowAdvanced',false, ...
    'limitedNodes',nodalMotion.nodalTrace.limitedNodeCount, ...
    'maximumNodalSpeedMS',max(abs(nodalMotion.nodalNormalSpeedMS)), ...
    'referenceCellFaceTargetResidualMS',max(abs(nodalMotion.projectionResidualMS)), ...
    'originalCellPressureMismatch',max(abs(cellState.pressure./actualJump.state.pressure-1)), ...
    'reconstructedFacePressureMismatch',max(abs(nodalMotion.reconstructedFacePrimitiveRP(:,5)./actualJump.state.pressure-1)), ...
    'maximumGCLRelative',max(moving.gclRelative),'sourceFileSHA256',rocketStreamingFileSHA256(capFile));
report=struct('status','tested-shared-node-traces-and-local-characteristic-motion-not-a-flow-solution', ...
    'linearPrimitivePointSampleError',linearError,'transformedSampleError',frameError,'constantQBitwise',true, ...
    'knownShockSpeedErrorBySoundSpeed',speedError,'characteristicFrameErrorBySoundSpeed',characteristicFrameError, ...
    'realProbe',probe,'productionQualified',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'nodal_motion_evidence.mat'),'report','box','gas','nodeIDs','base','G','state','linear','expected', ...
    'R','shift','lengthScale','boost','transformed','expectedFrame','constant','uniform','hot','limited', ...
    'up','realGas','normals','knownMach','knownSpeed','jump','characteristic','rotatedUp','rotatedNormals', ...
    'rotatedState','rotated','refusals','capFile','nodalMotion','moving','-v7');
fid=fopen(fullfile(folder,'nodal_motion_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));disp(report);
end
function id=reject(action,expected)
id='';try,action();catch ex,id=ex.identifier;end
assert(strcmp(id,expected),'ShockFit:NodalRefusalTest','实际拒绝%s，需要%s。',id,expected);
end
