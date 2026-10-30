function study=solveShockFittedPiston3D(upstream,gas,normal,pistonNormalVelocity,gridCounts,times)
%SOLVESHOCKFITTEDPISTON3D 内部求激波速度的移动活塞平面压缩基准及三维ALE推进。
% 模块一：由给定活塞法向速度与移动RH共同求未知激波速度，不指定激波后压力。
% 模块二：三维四面体域同时受移动活塞和移动激波约束，侧面为同一均匀态储库。
% 模块三：保存实际推进五量、壁压、激波相容差与能量账；这是解析基准，不是箭体解。
validateattributes(normal,{'double'},{'real','finite','size',[3 1]});
assert(abs(norm(normal)-1)<1e-12,'ShockFit:Normal','基准法向须为单位向量。');
validateattributes(pistonNormalVelocity,{'double'},{'real','finite','scalar'});
validateattributes(gridCounts,{'double'},{'real','finite','size',[1 3],'integer','positive'});
validateattributes(times,{'double'},{'real','finite','column','nonnegative'});
assert(numel(times)>=2&&times(1)==0&&all(diff(times)>0),'ShockFit:Time','需要从0开始严格递增的时间。');
initial=shockFittedState3D([upstream.density upstream.velocityVector' upstream.temperature],gas,'primitive');
un=dot(upstream.velocityVector,normal);delta=un-pistonNormalVelocity;
assert(delta>0,'ShockFit:Piston','本基准只接受压缩活塞；膨胀不是同一激波模型。');
lo=initial.soundSpeed*(1+2e-5);hi=initial.soundSpeed*1.5;
assert(mismatch(lo)<0,'ShockFit:PistonWeakLimit','活塞压缩低于本RH法向Mach可分辨范围。');
while mismatch(hi)<0
    hi=hi*1.25;assert(isfinite(hi),'ShockFit:PistonRoot','活塞激波根无法括区。');
end
for j=1:64
    mid=(lo+hi)/2;if mismatch(mid)>0,hi=mid;else,lo=mid;end
end
relativeSpeed=(lo+hi)/2;D=un-relativeSpeed;
jump=movingShockJump3D(upstream,gas,normal',D);
wallResidual=jump.state.velocity*normal-pistonNormalVelocity;
assert(abs(wallResidual)<2e-9,'ShockFit:PistonRoot','求得激波不能满足给定活塞无穿透条件。');
ex=-normal;trial=[0;1;0];if abs(ex'*trial)>.9,trial=[0;0;1];end
ey=trial-ex*(ex'*trial);ey=ey/norm(ey);frame=[ex ey cross(ex,ey)];
[xi,eta,zeta]=ndgrid(linspace(0,1,gridCounts(1)+1),linspace(-.015,.015,gridCounts(2)+1), ...
    linspace(-.012,.012,gridCounts(3)+1));local=[xi(:) eta(:) zeta(:)];
tetra=zeros(6*prod(gridCounts),4);n=0;index=reshape(1:numel(xi),size(xi));
for k=1:gridCounts(3),for j=1:gridCounts(2),for i=1:gridCounts(1)
    a=index(i,j,k);b=index(i+1,j,k);c=index(i+1,j+1,k);d=index(i,j+1,k);
    e=index(i,j,k+1);f=index(i+1,j,k+1);g=index(i+1,j+1,k+1);h=index(i,j+1,k+1);
    tetra(n+(1:6),:)=[a b c g;a c d g;a d h g;a h e g;a e f g;a f b g];n=n+6;
end,end,end
nodes0=positions(0);tetra=positive(tetra,nodes0);mesh=buildShockFittedTetraTopology(nodes0,tetra);
boundaryFaces=mesh.faces(mesh.boundaryFaces,:);xface=mean(reshape(local(boundaryFaces(:),1),size(boundaryFaces)),2);
kinds=ones(size(mesh.boundaryFaces));kinds(xface==0)=2;kinds(xface==1)=3;
boundary=struct('kind',kinds,'exteriorQ',repmat(jump.state.Q,numel(kinds),1),'upstream',upstream);
assert(any(kinds==2)&&any(kinds==3)&&any(kinds==1),'ShockFit:PistonGeometry','活塞、激波或侧储库面缺失。');
Nc=mesh.cellCount;Q=repmat(jump.state.Q,Nc,1);history=zeros(numel(times)-1,10);last=struct;
if strcmp(gas.thermalPropertyMode,'caloricallyPerfect')
    for j=2:numel(times)
        last=advanceShockFittedALE3D(Q,mesh,positions(times(j-1)),positions(times(j)),times(j)-times(j-1),boundary,gas);
        reference=repmat(jump.state.Q,Nc,1);
        qError=max(abs(last.Q-reference)./max(1,abs(reference)),[],'all');
        pError=max(abs(last.state.pressure/jump.state.pressure-1));
        wallError=max(abs(last.wallPressurePa(kinds==2)/jump.state.pressure-1));
        shockError=max(abs(last.shockInteriorPressureRelativeMismatch(kinds==3)));
        history(j-1,:)=[times(j),qError,pError,wallError,shockError,max(last.geometry.gclRelative), ...
            max(abs(last.conservationResidual)./last.conservationScale),min(last.state.temperature),max(last.state.temperature), ...
            last.budget.declaredTotalArrayUpperBytes];
        Q=last.Q;
    end
    assert(max(history(:,2:7),[],'all')<2e-10,'ShockFit:PistonUniform','真实ALE推进未保持解析均匀激波后状态。');
    advanced=true;
else
    % NASA7只求移动RH/活塞相容，不把未实现的NASA7 ALE通量标为已运行。
    advanced=false;history=[];
end
study=struct('upstream',upstream,'gas',gas,'normal',normal,'pistonNormalVelocityMS',pistonNormalVelocity, ...
    'shockNormalSpeedMS',D,'upstreamRelativeShockSpeedMS',relativeSpeed,'jump',jump,'wallVelocityResidualMS',wallResidual, ...
    'gridCounts',gridCounts,'times',times,'frame',frame,'mesh',mesh,'nodes0',nodes0,'nodes1',positions(times(end)), ...
    'boundary',boundary,'initialQ',repmat(jump.state.Q,Nc,1),'Q',Q,'lastStep',last,'history',history,'ALEAdvanced',advanced, ...
    'method','给定活塞速度+同cp移动RH标量根求激波速度；三维共享面ALE五量推进核对', ...
    'benchmarkShockSpeedSolved',true,'rocketShockShapeSolved',false,'bodyMatchedToRocket',false, ...
    'heatFluxSolved',false,'wallTemperatureSolved',false,'productionQualified',false,'newProductionSolvedAreaM2',0);

    function f=mismatch(w)
        compressed=shockLayerNormalJumpBatch(upstream,gas,w);
        f=w-compressed.normalSpeed-delta;
    end
    function X=positions(t)
        wall=-pistonNormalVelocity*t;shock=.02-D*t;
        coord=local;coord(:,1)=wall+local(:,1)*(shock-wall);X=coord*frame';
    end
end

function T=positive(T,X)
a=X(T(:,1),:);b=X(T(:,2),:);c=X(T(:,3),:);d=X(T(:,4),:);
negative=sum((b-a).*cross(c-a,d-a,2),2)<0;T(negative,[3 4])=T(negative,[4 3]);
end
