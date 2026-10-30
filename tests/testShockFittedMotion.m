function report=testShockFittedMotion(folder)
%TESTSHOCKFITTEDMOTION 检查共享节点全局投影，而非错误的独立节点平均。
% 模块一：可表示的非均匀速度精确回代；不可表示速度与独立稀疏QR比较。
% 模块二：真实壁面不动、内层随动与GCL检查，不授予物理区域资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupShockFittedFlow;
assert(~exist(folder,'dir'),'ShockFit:Overwrite','投影测试使用新目录。');mkdir(folder);
cap=buildShockFittedSphericalCap(.05,60,4,12,.008,3);shell=cap.shell;mesh=shell.mesh;
N=numel(shell.originalVertexIDs);outer=mesh.boundaryFaces(shell.outerBoundaryMask);
g=shockFittedALEGeometry(mesh,shell.nodes,shell.nodes,1);normal=g.normal(outer,:);area=g.averageAreaM2(outer);
faces=mesh.faces(outer,:)-shell.layers*N;d=shell.nodes(1:N,:)/cap.radiusM;
C=sparse(numel(outer),N);localNumerator=zeros(N,1);localDenominator=zeros(N,1);
target=10+3*sin((1:numel(outer))');
for j=1:3
    a=sum(d(faces(:,j),:).*normal,2);
    C=C+sparse((1:numel(outer))',faces(:,j),a/3,numel(outer),N);
    localNumerator=localNumerator+accumarray(faces(:,j),area.*a.*target,[N 1]);
    localDenominator=localDenominator+accumarray(faces(:,j),area.*a.^2,[N 1]);
end
known=8+2*cos((1:N)'*.31);representable=shockFittedShellMotion3D(shell,shell.nodes,C*known);
assert(max(abs(representable.faceNormalSpeedMS-C*known))<1e-9);
motion=shockFittedShellMotion3D(shell,shell.nodes,target);
A=spdiags(sqrt(area),0,numel(area),numel(area))*C;rhs=sqrt(area).*target;
reference=A\rhs;
assert(max(abs(motion.radialNodeSpeedMS-reference))<1e-8&&motion.normalEquationStationarityRelative<2e-10);
old=localNumerator./localDenominator;
assert(norm(A*motion.radialNodeSpeedMS-rhs)<=norm(A*old-rhs)*(1+2e-12));
X=shell.nodes+2e-7*motion.nodeVelocityMS;geometry=shockFittedALEGeometry(mesh,shell.nodes,X,2e-7);
assert(isequal(X(1:N,:),shell.nodes(1:N,:))&&max(geometry.gclRelative)<2e-11);
report=struct('status','tested-shared-node-face-centroid-global-weighted-least-squares', ...
    'representableMaximumFaceErrorMS',max(abs(representable.faceNormalSpeedMS-C*known)), ...
    'QRMaximumNodeDifferenceMS',max(abs(motion.radialNodeSpeedMS-reference)), ...
    'normalEquationStationarityRelative',motion.normalEquationStationarityRelative, ...
    'weightedResidualRmsMS',motion.weightedProjectionRmsMS,'oldLocalAverageResidualRmsMS',norm(A*old-rhs)/sqrt(sum(area)), ...
    'maximumGCLRelative',max(geometry.gclRelative),'bodyVerticesFixed',true, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'shock_motion_evidence.mat'),'report','cap','C','known','target','reference','old','motion','representable','geometry','-v7');
fid=fopen(fullfile(folder,'shock_motion_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));disp(report);
end
