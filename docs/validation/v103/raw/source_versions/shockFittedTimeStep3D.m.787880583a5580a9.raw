function estimate=shockFittedTimeStep3D(Q,mesh,nodes,nodeVelocity,gas,CFL,boundary)
%SHOCKFITTEDTIMESTEP3D 按真实移动面声学谱半径给出全五量显式时间步预算。
% 模块一：内部面同时包含初始特征和精确Riemann波扇，不漏掉强压缩波的较快传播。
% 模块二：每单元累加面积×波速，dt=CFL*min(V/ΣAλ)，并保留实际控制单元。
% 模块三：边界显式传入后核对储库/镜像墙/上游激波/超声速外流的波速。
% 模块四：这是当前几何的显式谱预算，不替代整步几何/状态正性及温域检查。
if nargin<7,boundary=[];end
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:CFLGas', ...
    '当前精确波扇CFL仅常gamma，不回退NASA7。');
validateattributes(CFL,{'double'},{'real','finite','positive','scalar','<=',1});
assert(isequal(size(nodeVelocity),size(nodes))&&all(isfinite(nodeVelocity),'all'), ...
    'ShockFit:MeshVelocity','需要每个共享节点唯一的三维网格速度。');
s=shockFittedState3D(Q,gas,'conservative');g=shockFittedALEGeometry(mesh,nodes,nodes,1);f=mesh.faces;
w=(nodeVelocity(f(:,1),:)+nodeVelocity(f(:,2),:)+nodeVelocity(f(:,3),:))/3;
D=sum(w.*g.normal,2);owner=mesh.owner;
radius=abs(sum(s.velocity(owner,:).*g.normal,2)-D)+s.soundSpeed(owner);
inside=find(mesh.neighbour>0);next=mesh.neighbour(inside);
other=abs(sum(s.velocity(next,:).*g.normal(inside,:),2)-D(inside))+s.soundSpeed(next);
radius(inside)=max(radius(inside),other);
waves=nan(mesh.faceCount,2);pressureStar=nan(mesh.faceCount,1);
if ~isempty(inside),includeWaves(inside,Q(next,:));end
boundaryWavesIncluded=~isempty(boundary);
if boundaryWavesIncluded
    assert(isequal(size(boundary.kind),[numel(mesh.boundaryFaces) 1])&& ...
        all(ismember(boundary.kind,[1 2 3 4 5])),'ShockFit:Boundary','CFL需要完整明确的边界类型。');
    for kind=1:5
        pick=find(boundary.kind==kind);if isempty(pick),continue;end
        ids=mesh.boundaryFaces(pick);right=Q(owner(ids),:);
        if kind==1
            right=boundary.exteriorQ(pick,:);
        elseif kind==2
            state=s.velocity(owner(ids),:);relative=sum(state.*g.normal(ids,:),2)-D(ids);
            mirror=shockFittedState3D([s.density(owner(ids)),state-2*relative.*g.normal(ids,:), ...
                s.temperature(owner(ids))],gas,'primitive');right=mirror.Q;
        elseif kind==3||kind==5
            up=shockFittedState3D([boundary.upstream.density boundary.upstream.velocityVector' ...
                boundary.upstream.temperature],gas,'primitive');right=repmat(up.Q,numel(ids),1);
        end
        includeWaves(ids,right);
    end
end
rates=abs(mesh.incidence)*(g.averageAreaM2.*radius);
assert(all(rates>0&isfinite(rates)),'ShockFit:CFL','声学谱率退化。');
[minimum,cellID]=min(g.volume0./rates);
estimate=struct('dt',CFL*minimum,'CFL',CFL,'limitingCell',cellID,'faceSpectralRadiusMS',radius, ...
    'faceGridNormalSpeedMS',D,'cellVolumeM3',g.volume0,'cellSpectralRateM3PerS',rates, ...
    'relativeRiemannWaveBoundsMS',waves,'riemannPressureStarPa',pressureStar, ...
    'boundaryRiemannWavesIncluded',boundaryWavesIncluded, ...
    'method','真实共享面ALE初始特征和精确Riemann波扇包络/单元体积的显式CFL预算', ...
    'positivityProvenByCFLAlone',false,'physicalQualification',false);

    function includeWaves(ids,rightQ)
        left=Q(owner(ids),:);normal=g.normal(ids,:);
        [~,audit]=shockFittedALEFlux3D(left,rightQ,normal,D(ids),gas);
        rightState=shockFittedState3D(rightQ,gas,'conservative');star=audit.riemann.pressureStarPa;
        pL=s.pressure(owner(ids));pR=rightState.pressure;gm=gas.gamma;
        factorL=ones(numel(ids),1);factorR=factorL;
        shockL=star>pL;shockR=star>pR;
        factorL(shockL)=sqrt((gm+1)/(2*gm)*star(shockL)./pL(shockL)+(gm-1)/(2*gm));
        factorR(shockR)=sqrt((gm+1)/(2*gm)*star(shockR)./pR(shockR)+(gm-1)/(2*gm));
        unL=sum(s.velocity(owner(ids),:).*normal,2)-D(ids);
        unR=sum(rightState.velocity.*normal,2)-D(ids);
        bounds=[unL-s.soundSpeed(owner(ids)).*factorL,unR+rightState.soundSpeed.*factorR];
        initialRight=abs(unR)+rightState.soundSpeed;
        radius(ids)=max([radius(ids),initialRight,abs(bounds),abs(audit.riemann.normalVelocityStarMs)],[],2);
        waves(ids,:)=bounds;pressureStar(ids)=star;
    end
end
