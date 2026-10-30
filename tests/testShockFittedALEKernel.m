function report=testShockFittedALEKernel(folder,sourceFile)
%TESTSHOCKFITTEDALEKERNEL 移动激波、三维ALE通量/GCL、真实原面几何和活塞双边界基准。
% 模块一：两物性RH、实验室能量、旋转/Galilean、壁面Riemann和严格拒绝路径。
% 模块二：内求激波速度的移动活塞及三维共享面实际推进，分别保存两物性适用范围。
% 模块三：真实20°原箭鼻部只检查体网格/GCL/滑移接口，不冒充求得真实激波或壁面场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupShockFittedFlow;
assert(~exist(folder,'dir'),'ShockFit:Overwrite','必须写入新验证目录。');mkdir(folder);
source=load(sourceFile,'fs','geometry','seed');fs=source.fs;geometry=source.geometry;
up=struct('pressure',double(fs.pressure),'temperature',double(fs.temperature),'density',double(fs.density), ...
    'velocityVector',double(fs.flowVelocityVector(:)),'method','原19752面M5/30km/alpha20均匀上游，仅作明确输入');
modes={'caloricallyPerfect','thermallyPerfect'};gases=cell(2,1);jumpCases=cell(2,1);pistons=cell(2,1);
e1=up.velocityVector/norm(up.velocityVector);e2=[0;1;0];frame=[e1 e2 cross(e1,e2)];
coordinates=[1 0 0;.96 .28 0;.91 -.31 sqrt(1-.91^2-.31^2)];normals=coordinates*frame';
D=[0;75;-90];boost=[43;-29;61];angle=.37;rotation=[cos(angle) 0 sin(angle);0 1 0;-sin(angle) 0 cos(angle)];
for mode=1:2
    gas=createSurfaceFlowGas(modes{mode},double(fs.R),1.4);gases{mode}=gas;
    jump=movingShockJump3D(up,gas,normals,D);
    shifted=up;shifted.velocityVector=up.velocityVector+boost;
    galilean=movingShockJump3D(shifted,gas,normals,D+normals*boost);
    turned=up;turned.velocityVector=rotation*up.velocityVector;
    rotated=movingShockJump3D(turned,gas,normals*rotation',D);
    assert(max(abs(galilean.state.temperature-jump.state.temperature))<2e-9&& ...
        max(abs(galilean.state.pressure./jump.state.pressure-1))<2e-11);
    assert(max(abs(galilean.state.velocity-jump.state.velocity-boost'),[],'all')<2e-9);
    assert(max(abs(rotated.state.velocity-jump.state.velocity*rotation'),[],'all')<2e-9);
    assert(max(abs(rotated.state.pressure./jump.state.pressure-1))<2e-11);
    assert(abs(jump.laboratoryTotalEnthalpyChange(2))>1e3, ...
        'ShockFit:MovingEnergy','测试没有分辨移动激波的实验室系总焓变化。');
    inverse=shockFittedState3D(jump.state.Q,gas,'conservative');
    assert(max(abs(inverse.temperature-jump.state.temperature))<2e-9);
    jumpCases{mode}=struct('base',jump,'galilean',galilean,'rotated',rotated,'boost',boost,'rotation',rotation);
    pistons{mode}=solveShockFittedPiston3D(up,gas,normals(2,:)',dot(up.velocityVector,normals(2,:)')-450, ...
        [2 2 2],linspace(0,2e-6,9)');
    if mode==1
        a=sqrt(gas.gamma*gas.R*up.temperature);delta=450;term=(gas.gamma+1)*delta/4;
        exact=term+sqrt(a^2+term^2);
        assert(abs(pistons{mode}.upstreamRelativeShockSpeedMS-exact)<2e-9);
    else
        assert(~pistons{mode}.ALEAdvanced&&isempty(pistons{mode}.history));
    end
    fprintf('MOVING_SHOCK mode=%s maxRH=%.6g pistonRelativeSpeed=%.12g ALEAdvanced=%d\n', ...
        modes{mode},max(abs(jump.conservationRelative),[],'all'),pistons{mode}.upstreamRelativeShockSpeedMS,pistons{mode}.ALEAdvanced);
end

%% 二、不同左右态的三维ALE通量变换、镜像壁面压力和共享面非线性几何运动
gas=gases{1};left=shockFittedState3D([.4 60 -25 35 400;.7 -100 50 -10 650],gas,'primitive');
right=shockFittedState3D([.3 -20 40 15 550;.5 -30 -20 80 350],gas,'primitive');
n=[1 0 0;1 2 -2]/3;n(1,:)=[1 0 0];gridSpeed=[35;-60];
[flux,fluxAudit]=shockFittedALEFlux3D(left.Q,right.Q,n,gridSpeed,gas);
leftBoost=shockFittedState3D([left.density left.velocity+boost' left.temperature],gas,'primitive');
rightBoost=shockFittedState3D([right.density right.velocity+boost' right.temperature],gas,'primitive');
shiftFlux=shockFittedALEFlux3D(leftBoost.Q,rightBoost.Q,n,gridSpeed+n*boost,gas);
expected=[flux(:,1),flux(:,2:4)+flux(:,1).*boost', ...
    flux(:,5)+flux(:,2:4)*boost+.5*sum(boost.^2)*flux(:,1)];
assert(max(abs(shiftFlux-expected)./max(1,abs(expected)),[],'all')<2e-11);
rotFlux=shockFittedALEFlux3D(rotateQ(left.Q,rotation),rotateQ(right.Q,rotation),n*rotation',gridSpeed,gas);
rotExpected=[flux(:,1),flux(:,2:4)*rotation',flux(:,5)];
assert(max(abs(rotFlux-rotExpected)./max(1,abs(rotExpected)),[],'all')<2e-11);
wallState=shockFittedState3D([1 55 80 10 300;1 -5 20 -40 300;1 25 0 30 300],gas,'primitive');
wallNormals=repmat([1 0 0],3,1);wallD=25*ones(3,1);
[wallFlux,wallAudit]=shockFittedSlipWallFlux3D(wallState.Q,wallNormals,wallD,gas);
assert(wallAudit.pressurePa(1)>wallState.pressure(1)&&wallAudit.pressurePa(2)<wallState.pressure(2)&& ...
    abs(wallAudit.pressurePa(3)/wallState.pressure(3)-1)<2e-12);
assert(all(wallFlux(:,1)==0)&&isequal(wallFlux(:,5),wallAudit.pressurePa.*wallD));
base=pistons{1};X0=base.nodes0;affine=[1.04 .03 -.02;-.01 .98 .015;.02 -.01 1.03];
X1=X0*affine'+[.001 -.002 .003];motion=shockFittedALEGeometry(base.mesh,X0,X1,.01);
velocity=(X1-X0)/.01;f=base.mesh.faces;faceVelocity=(velocity(f(:,1),:)+velocity(f(:,2),:)+velocity(f(:,3),:))/3;
wrongSwept=.01*sum(motion.areaVector1.*faceVelocity,2);
wrongGCL=abs(motion.volume1-motion.volume0-base.mesh.incidence*wrongSwept)./max(motion.volume0,motion.volume1);
assert(max(wrongGCL)>1e-5,'ShockFit:GCLOracle','测试没有暴露端点面积替代时间积分的错误。');
geometricScale=2.3;scaled=shockFittedALEGeometry(base.mesh,X0*geometricScale,X1*geometricScale,.01*geometricScale);
% 逐单元比值，避免矩阵右除把尺度测试退化成一个投影。
assert(max(abs(scaled.volume1./motion.volume1/geometricScale^3-1))<2e-12);
assert(max(abs(scaled.gridNormalSpeedMS-motion.gridNormalSpeedMS))<2e-10);
transformed=shockFittedALEGeometry(base.mesh,X0*rotation'+[.2 -.3 .4],X1*rotation'+[.2 -.3 .4],.01);
assert(max(abs(transformed.volume1./motion.volume1-1))<2e-11);
fluxChecks=struct('left',left,'right',right,'normal',n,'gridNormalSpeed',gridSpeed,'flux',flux,'audit',fluxAudit, ...
    'shiftFlux',shiftFlux,'expectedShiftFlux',expected,'rotatedFlux',rotFlux,'expectedRotatedFlux',rotExpected, ...
    'wallState',wallState,'wallNormals',wallNormals,'wallNormalSpeed',wallD,'wallFlux',wallFlux,'wallAudit',wallAudit);
motionChecks=struct('nodes0',X0,'nodes1',X1,'mesh',base.mesh,'dt',.01,'geometry',motion, ...
    'wrongEndpointAreaGCLRelative',wrongGCL,'scaled',scaled,'scale',geometricScale,'transformed',transformed,'rotation',rotation);

%% 三、真实20°原箭鼻部非轴对称面网格；几何均匀态校核不作真实流场
V=double(geometry.V);F=double(geometry.F);centers=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
ids=find(centers(:,1)>max(V(:,1))-.002);shell=buildShockFittedSurfaceShell(V,F,double(ids),.008,2);
N=numel(shell.originalVertexIDs);X0=shell.nodes;X1=X0;
for layer=1:shell.layers
    indices=layer*N+(1:N);coefficient=sin((1:N)'*.37).*cos((1:N)'*.19);
    X1(indices,:)=X1(indices,:)+layer/shell.layers*1e-6*coefficient.*shell.vertexNormals;
end
ref=shockFittedState3D([up.density up.velocityVector' up.temperature],gas,'primitive');
boundary=struct('kind',ones(numel(shell.mesh.boundaryFaces),1), ...
    'exteriorQ',repmat(ref.Q,numel(shell.mesh.boundaryFaces),1),'upstream',up);
uniform=advanceShockFittedALE3D(repmat(ref.Q,shell.mesh.cellCount,1),shell.mesh,X0,X1,1e-8,boundary,gas);
uniformError=max(abs(uniform.Q-ref.Q)./max(1,abs(ref.Q)),[],'all');
assert(uniformError<2e-10&&shell.bodyVerticesExact&&~shell.shockGeometrySolved);
wallFaces=shell.mesh.boundaryFaces(shell.wallBoundaryMask);bodyNormals=uniform.geometry.normal(wallFaces,:);
wallVelocity=repmat(up.velocityVector',numel(wallFaces),1);
tangent=wallVelocity-sum(wallVelocity.*bodyNormals,2).*bodyNormals;
slipState=shockFittedState3D([repmat(up.density,numel(wallFaces),1),tangent,repmat(up.temperature,numel(wallFaces),1)],gas,'primitive');
[slipFlux,slipAudit]=shockFittedSlipWallFlux3D(slipState.Q,bodyNormals,zeros(numel(wallFaces),1),gas);
assert(max(abs(slipAudit.pressurePa/up.pressure-1))<2e-11&&all(slipFlux(:,[1 5])==0,'all'));
bodyCheck=struct('shell',shell,'nodes1',X1,'uniformStep',uniform,'uniformStateRelativeError',uniformError, ...
    'slipInteriorState',slipState,'slipFlux',slipFlux,'slipAudit',slipAudit, ...
    'sourceMach',fs.Mach,'sourceAlphaDeg',fs.alpha,'testRole','真实几何的共享面/GCL/均匀态与局部滑移恒等式；不是该M5/20度箭体外流', ...
    'pressureOnBodySolved',false,'shockShapeSolved',false,'heatFluxSolved',false,'productionQualified',false);

%% 四、错误气体、错误法向、温域越界、内途中途翻转和预算拒绝；均不可静默修补
rejections=struct;
rejections.NASA7Flux=reject(@()shockFittedALEFlux3D(left.Q,right.Q,n,gridSpeed,gases{2}),'ShockFit:FluxGas');
rejections.normal=reject(@()movingShockJump3D(up,gas,normals*2,D),'ShockFit:Normal');
rejections.subsonic=reject(@()movingShockJump3D(up,gas,normals(1,:),dot(up.velocityVector,normals(1,:)')),'ShockPatch:NormalMach');
rejections.temperature=reject(@()shockFittedState3D([1 0 0 0 2501],gas,'primitive'),'ShockPatch:Temperature');
tet=[0 0 0;1 0 0;0 1 0;0 0 1];mesh=buildShockFittedTetraTopology(tet,[1 2 3 4]);
rejections.intermediateInversion=reject(@()shockFittedALEGeometry(mesh,tet,tet*diag([-4 -2 1]),1),'ShockFit:CellVolume');
rejections.memory=reject(@()advanceShockFittedALE3D(base.initialQ,base.mesh,base.nodes0,base.nodes1,base.times(end),base.boundary,gas,1000),'ShockFit:MemoryBudget');
report=struct('status','tested-moving-shock-ALE-GCL-wall-piston-and-real-mesh-kernel-with-explicit-scope', ...
    'previousGoalTurn','progress-storage-maintenance-delivered','completeKernelValidation',true, ...
    'movingRHGasModes',{modes},'ALEGasModes',{{'caloricallyPerfect'}}, ...
    'originalBodyFaces',size(F,1),'originalNoseFaces',numel(ids),'realTetrahedra',shell.mesh.cellCount, ...
    'uniformStateRelativeError',uniformError,'maximumGCLRelative',max(uniform.geometry.gclRelative), ...
    'pistonMaximumStateRelativeError',max(pistons{1}.history(:,2)), ...
    'pistonMaximumConservationRelative',max(pistons{1}.history(:,7)), ...
    'largestDeclaredArrayBudgetBytes',uniform.budget.declaredTotalArrayUpperBytes, ...
    'originalSourceFile',sourceFile,'originalSourceSHA256',rocketStreamingFileSHA256(sourceFile), ...
    'rocketShockShapeSolved',false,'rocketPressureSolved',false,'heatFluxSolved',false,'wallTemperatureSolved',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'shock_fitted_kernel_evidence.mat'),'report','up','geometry','gases','jumpCases','pistons', ...
    'fluxChecks','motionChecks','bodyCheck','rejections','-v7');
fid=fopen(fullfile(folder,'shock_fitted_kernel_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
disp(report);
end

function Q=rotateQ(Q,R)
Q(:,2:4)=Q(:,2:4)*R';
end
function id=reject(action,expected)
id='';try,action();catch ex,id=ex.identifier;end
assert(strcmp(id,expected),'ShockFit:RefusalTest','拒绝身份不符：实际%s，需要%s。',id,expected);
end
