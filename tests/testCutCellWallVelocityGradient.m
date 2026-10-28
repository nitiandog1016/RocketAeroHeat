function report=testCutCellWallVelocityGradient(folder)
%TESTSTATICWALLVELOCITYDRAFT 三维静止墙梯度草稿的仿射/旋转/牵引与拒绝校核。
% 模块一：已知静止平壁仿射速度真解，允许不同真实法向和大的切向偏移。
% 模块二：三维坐标旋转、单面/批量一致和Newton牵引非负机械耗散。
% 模块三：非法/未支持数据必须明确拒绝；仅单位校核，不是完整NS或区域热流解。
addpath(fileparts(mfilename('fullpath')));
assert(~exist(folder,'dir'),'WallVelocityDraftTest:Overwrite','不得覆盖已有试算。');mkdir(folder);
normal=[1 0 0;1 2 -1;-.4 .2 1;0 -1 0];normal=normal./vecnorm(normal,2,2);N=size(normal,1);
wallCenter=[.2 .3 .4;-.5 .1 .2;.7 -.2 .4;-.1 .4 .3];distance=[.02;1e-5;.13;.003];
offset=[.6 -.3 .4;.2 .4 .1;-.1 .5 -.3;.4 .7 .2];offset=offset-sum(offset.*normal,2).*normal;
cellCenter=wallCenter-distance.*normal+offset;dn=sum((wallCenter-cellCenter).*normal,2);
coefficient=[100 -30 20;-40 70 15;80 -20 -60;13 25 -10];Ucell=-dn.*coefficient;Uwall=zeros(N,3);
boundary=struct('type','wall','velocity',[0 0 0],'velocityGradient','dirichlet-tangent');ids=(1:N)';
[gradient,audit]=cutCellWallVelocityGradient(cellCenter,wallCenter,normal,Ucell,Uwall,boundary,ids);
reference=zeros(3,3,N);for j=1:N,reference(:,:,j)=coefficient(j,:)'*normal(j,:);end
affineError=max(abs(gradient-reference),[],'all')/max(1,max(abs(reference),[],'all'));assert(affineError<2e-12);
for j=1:N
    scalar=cutCellWallVelocityGradient(cellCenter(j,:),wallCenter(j,:),normal(j,:),Ucell(j,:),Uwall(j,:),boundary,j);
    assert(isequal(scalar,gradient(:,:,j)));
end
axis=[.3 .5 -.2]';axis=axis/norm(axis);K=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];angle=.68;
rotation=cos(angle)*eye(3)+(1-cos(angle))*(axis*axis')+sin(angle)*K;
rotated=cutCellWallVelocityGradient(cellCenter*rotation',wallCenter*rotation',normal*rotation',Ucell*rotation',Uwall,boundary,ids);
rotatedReference=gradient;for j=1:N,rotatedReference(:,:,j)=rotation*gradient(:,:,j)*rotation';end
rotationError=max(abs(rotated-rotatedReference),[],'all')/max(1,max(abs(rotatedReference),[],'all'));assert(rotationError<2e-9);
mu=1.8e-5;traction=zeros(N,3);expected=traction;
for j=1:N
    G=gradient(:,:,j);stress=mu*(G+G'-(2/3)*trace(G)*eye(3));traction(j,:)=(stress*normal(j,:)')';
    expected(j,:)=-mu/dn(j)*(Ucell(j,:)+dot(Ucell(j,:),normal(j,:))*normal(j,:)/3);
end
tractionError=max(abs(traction-expected),[],'all')/max(1,max(abs(expected),[],'all'));assert(tractionError<2e-13);
assert(all(-sum(Ucell.*traction,2)>=0));
legacy=rmfield(boundary,'velocityGradient');assert(strcmp(cutCellWallVelocityPolicy(legacy),'legacy'));
legacy.velocityGradient='legacy';assert(strcmp(cutCellWallVelocityPolicy(legacy),'legacy'));
unknown=boundary;unknown.velocityGradient='wrong';nonzero=boundary;nonzero.velocity=[1 0 0];
functionWall=boundary;functionWall.velocity=@(~,~)[0 0 0];extra=boundary;extra.velocitySurfaceGradient=zeros(3);
bad={@()cutCellWallVelocityPolicy(unknown),@()cutCellWallVelocityPolicy(nonzero), ...
     @()cutCellWallVelocityPolicy(functionWall),@()cutCellWallVelocityPolicy(extra), ...
     @()cutCellWallVelocityGradient(cellCenter,wallCenter,2*normal,Ucell,Uwall,boundary,ids), ...
     @()cutCellWallVelocityGradient(wallCenter,wallCenter,normal,Ucell,Uwall,boundary,ids), ...
     @()cutCellWallVelocityGradient(cellCenter,wallCenter,normal,Ucell,ones(N,3),boundary,ids), ...
     @()cutCellWallVelocityGradient(cellCenter,wallCenter,normal,Ucell,Uwall,boundary,ones(N,1))};
for j=1:numel(bad)
    rejected=false;try,bad{j}();catch ex,rejected=startsWith(ex.identifier,'CutNS:WallVelocity');end;assert(rejected);
end
report=struct('status','passed-draft-3D-static-wall-affine-rotation-traction-and-explicit-rejections', ...
    'affineRelativeDifference',affineError,'rotationRelativeDifference',rotationError, ...
    'tractionRelativeDifference',tractionError,'explicitRejections',numel(bad),'singleBatchIdentical',true, ...
    'integratedIntoNS',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'wall_velocity_unit_evidence.mat'),'report','gradient','audit','reference','rotation','rotated','traction','expected','-v7');
fid=fopen(fullfile(folder,'wall_velocity_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
