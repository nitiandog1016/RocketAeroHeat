function [full,graph,detail]=splitRocketGraphDirection(cache,operator,directionU)
%SPLITROCKETGRAPHDIRECTION 分解完整选中分支J与已装配Jgraph在同一方向上的作用。
% 模块一：沿当前Q/几何/边界缓存计算Euler、完整WLS应力、导热和两项黏性功。
% 模块二：单独重建图预条件器的冻结两点扩散/功速度及外储库Rusanov近似。
% 模块三：所有分量采用J而非-J；时间对角从B=I/tau-Jgraph中明确消去。
% 本函数仅用于有界诊断，不修改Q，不构造新预条件器，不授予物理或覆盖资格。
Q=cache.Q;cut=cache.geometry;gas=cache.gas;scale=cache.scale;f=cut.openFaces;w=cut.walls;
n=size(Q,1);nf=numel(f.area);nw=numel(w.area);v=cut.cells.volume(:);
assert(isequal(size(directionU),size(Q))&&all(isfinite(directionU),'all'));
assert(operator.sourceFaces*operator.layers==n&&isequal(operator.residualScale,scale)&& ...
    strcmp(operator.internalViscousLinearization,'frozen-two-point')&& ...
    strcmp(operator.normalEulerLinearization,'exact-branch')&&strcmp(operator.sideEulerLinearization,'exact-branch'), ...
    'GraphSplit:Operator','本诊断要求当前明确的内部精确Euler/冻结两点图预条件器。');
assert(all(f.cellMinus>0)&&all(w.ownerCell>0)&&all(~cache.state.wall), ...
    'GraphSplit:Topology','本诊断限定完整原面层域及非wall外边界，不静默跳过其他拓扑。');
dq=directionU.*scale;rho=Q(:,1);u=Q(:,2:4)./rho;cv=gas.R/(gas.gamma-1);
T=(Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);du=(dq(:,2:4)-u.*dq(:,1))./rho;
dT=(dq(:,5)-sum(u.*dq(:,2:4),2)+(.5*sum(u.^2,2)-cv*T).*dq(:,1))./(rho*cv);
dP=[dq(:,1) du dT];given=zeros(nf,5);ids=cache.outsideFaceIDs;
given(ids,:)=reshape(pagemtimes(cache.boundaryJacobian,permute(dP(cache.outsideOwner,:),[2 3 1])),5,[])';
bd=struct('primitiveDerivative',given,'incoming',cache.state.incoming,'wall',cache.state.wall,'boundarySensitivityIncluded',true);
[viscous,vd]=cutCellWlsViscousDirection(Q,directionU,cache.baseInfo,cut,gas,cache.boundary,scale,cache.gradient,cache.state,bd);
[euler,ed]=cutCellSelectedEulerDirection(Q,directionU,cache.baseInfo,cut,gas,cache.boundary,scale,cache.state,bd,cache.waves);

%% 模块四：内部法向/侧向、特征外边界与固壁Euler贡献分别保留。
inside=f.cellPlus>0;normal=inside&f.kind==1;side=inside&f.kind==3;outside=~inside;
assert(all(normal|side|outside),'GraphSplit:Topology','存在没有显式分类的共享面。');
full=struct;full.interiorNormalEuler=eulerPart(normal,false);full.interiorSideEuler=eulerPart(side,false);
full.outerEuler=eulerPart(outside,false);full.solidWallEuler=eulerPart(false(nf,1),true);
full.momentumStress=vd.momentumU;full.heatConduction=vd.heatU;full.workGradient=vd.workGradientU;full.workVelocity=vd.workVelocityU;
fullTotal=total(full);assert(relative(fullTotal,euler+viscous)<2e-13,'GraphSplit:FullSum','完整分量之和不等于原模块作用。');
graph=full;graph.momentumStress=zeros(n,5);graph.heatConduction=zeros(n,5);graph.workGradient=zeros(n,5);graph.workVelocity=zeros(n,5);

%% 模块五：内部Euler沿用完全相同的波系；外储库冻结，不传播特征状态导数。
out=find(outside);m=double(f.cellMinus(out));normalOut=f.normal(out,:);
[A,~,P]=cutCellConservativeDerivatives(Q(m,:),normalOut,gas);
wave=abs(sum(P(:,2:4).*normalOut,2))+sqrt(gas.gamma*gas.R*P(:,5));
outerFlux=.5*reshape(pagemtimes(A,permute(dq(m,:),[2 3 1])),5,[])'+.5*wave.*dq(m,:);
face=zeros(nf,5);face(out,:)=outerFlux;graph.outerEuler=-divergence(face,zeros(nw,5))./scale;

%% 模块六：两点应力/导热与冻结功速度，保留所有真实共享面及原固壁。
stress=zeros(nf,3);heat=zeros(nf,1);work=zeros(nf,1);k=gas.mu*gas.gamma*cv/gas.Pr;
for first=1:2048:nf
    ids=first:min(first+2047,nf);m=double(f.cellMinus(ids));p=double(f.cellPlus(ids));in=p>0;
    delta=f.centroid(ids,:)-cut.cells.centroid(m,:);delta(in,:)=cut.cells.centroid(p(in),:)-cut.cells.centroid(m(in),:);
    distance=sum(delta.*f.normal(ids,:),2);assert(all(distance>0),'GraphSplit:Distance','两点近似真实法向距非正。');
    dv=-du(m,:);dt=-dT(m);velocity=u(m,:);
    dv(in,:)=du(p(in),:)-du(m(in),:);dt(in)=dT(p(in))-dT(m(in));velocity(in,:)=.5*(u(m(in),:)+u(p(in),:));
    norm=f.normal(ids,:);traction=gas.mu./distance.*(dv+sum(dv.*norm,2).*norm/3);
    stress(ids,:)=traction;heat(ids)=k*dt./distance;work(ids)=sum(velocity.*traction,2);
end
owner=double(w.ownerCell);distance=sum((w.centroid-cut.cells.centroid(owner,:)).*w.normal,2);assert(all(distance>0));
wallStress=-gas.mu./distance.*(du(owner,:)+sum(du(owner,:).*w.normal,2).*w.normal/3);wallHeat=-k*dT(owner)./distance;
graph.momentumStress(:,2:4)=divergence(stress,wallStress);graph.momentumStress=graph.momentumStress./scale;
graph.heatConduction(:,5)=divergence(heat,wallHeat);graph.heatConduction=graph.heatConduction./scale;
graph.workGradient(:,5)=divergence(work,zeros(nw,1));graph.workGradient=graph.workGradient./scale;
graphTotal=total(graph);assembled=directionU/operator.pseudoTimeStepSeconds-applyRocketShellGraphOperator(operator,directionU);
assembledError=relative(graphTotal,assembled);assert(assembledError<2e-11,'GraphSplit:Assembly','分项重建不等于保存B去除时间对角后的作用。');
detail=struct('fullTotal',fullTotal,'graphTotal',graphTotal,'assembledGraphTotal',assembled,'graphAssemblyRelativeDifference',assembledError, ...
    'boundaryDirection',bd,'eulerDetail',ed,'viscousDetail',vd,'fullComponentSumRelativeDifference',relative(fullTotal,euler+viscous), ...
    'completeSelectedActionApplications',1,'newActualNSEvaluations',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0);
    function value=eulerPart(mask,wall)
        face=zeros(nf,5);face(mask,:)=ed.faceEulerFluxDerivative(mask,:);wallFlux=zeros(nw,5);
        if wall,wallFlux=ed.wallEulerFluxDerivative;end
        value=-divergence(face,wallFlux)./scale;
    end
    function out=divergence(faceFlux,wallFlux)
        out=zeros(n,size(faceFlux,2));m=double(f.cellMinus);p=double(f.cellPlus);in=p>0;
        for a=1:size(faceFlux,2)
            out(:,a)=accumarray(m,f.area./v(m).*faceFlux(:,a),[n 1])- ...
                accumarray(p(in),f.area(in)./v(p(in)).*faceFlux(in,a),[n 1])+ ...
                accumarray(double(w.ownerCell),w.area./v(w.ownerCell).*wallFlux(:,a),[n 1]);
        end
    end
end
function value=total(parts)
names=fieldnames(parts);value=zeros(size(parts.(names{1})));for j=1:numel(names),value=value+parts.(names{j});end
end
function e=relative(a,b)
e=norm(a-b,'fro')/max(1,norm(b,'fro'));
end
