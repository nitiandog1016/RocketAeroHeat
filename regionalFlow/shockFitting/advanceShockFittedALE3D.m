function step=advanceShockFittedALE3D(Q,mesh,nodes0,nodes1,dt,boundary,gas,maximumArrayBytes)
%ADVANCESHOCKFITTEDALE3D 共享面三维五量ALE有限体积一步，含移动壁面和显式激波边界。
% 模块一：内部/储库面用移动系精确Riemann，壁压由镜像Riemann求出；类型4仅允许超声速外流。
% 模块二：激波边界使用上游与实际网格速度的RH通量，另存内部压力相容差，不假装已匹配。
% 模块三：以新旧真实体积更新VQ，核对全域五量账；非正态明确拒绝。
% 模块四：声明缓存和Riemann临时工作数组预算；它不是MATLAB进程峰值测量。
if nargin<8,maximumArrayBytes=128*1024^2;end
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:FluxGas','当前一步内核仅常gamma，不回退NASA7。');
assert(isequal(size(Q),[mesh.cellCount 5]),'ShockFit:State','单元五量数目不符。');
validateattributes(maximumArrayBytes,{'double'},{'real','finite','positive','scalar'});
Nb=numel(mesh.boundaryFaces);
assert(isstruct(boundary)&&all(isfield(boundary,{'kind','exteriorQ','upstream'}))&& ...
    isequal(size(boundary.kind),[Nb 1])&&isequal(size(boundary.exteriorQ),[Nb 5])&& ...
    all(ismember(boundary.kind,[1 2 3 4])),'ShockFit:Boundary','边界须显式为储库1/壁面2/激波3/超声速外流4。');
inputs=whos('Q','mesh','nodes0','nodes1','boundary','gas');inputBytes=sum([inputs.bytes]);
% exactEuler最多预留1000个每面double标量，包含左右态/根/波扇/审计与表达式临时副本。
% 几何预留每单元200/每面120/每节点24个double，另加输出与更新五量副本。
workspace=8*(1000*mesh.faceCount+200*mesh.cellCount+120*mesh.faceCount+24*mesh.nodeCount);
budget=struct('retainedInputBytes',inputBytes,'declaredTemporaryArrayUpperBytes',workspace, ...
    'declaredTotalArrayUpperBytes',inputBytes+workspace,'maximumArrayBytes',maximumArrayBytes, ...
    'fullJacobianBuilt',false,'processPeakMemoryMeasured',false);
assert(budget.declaredTotalArrayUpperBytes<=maximumArrayBytes,'ShockFit:MemoryBudget','声明的全部工作数组预算超限；不启动推进。');
before=shockFittedState3D(Q,gas,'conservative'); %#ok<NASGU>
geo=shockFittedALEGeometry(mesh,nodes0,nodes1,dt);Nf=mesh.faceCount;
flux=nan(Nf,5);interior=find(mesh.neighbour>0);
if ~isempty(interior)
    flux(interior,:)=shockFittedALEFlux3D(Q(mesh.owner(interior),:),Q(mesh.neighbour(interior),:), ...
        geo.normal(interior,:),geo.gridNormalSpeedMS(interior),gas);
end
wallPressure=nan(Nb,1);shockPressure=nan(Nb,1);shockMismatch=nan(Nb,1);wallAudit=struct;shockAudit=struct;outflowAudit=struct;
for kind=1:4
    selected=find(boundary.kind==kind);if isempty(selected),continue;end
    faces=mesh.boundaryFaces(selected);inside=Q(mesh.owner(faces),:);normal=geo.normal(faces,:);D=geo.gridNormalSpeedMS(faces);
    if kind==1
        flux(faces,:)=shockFittedALEFlux3D(inside,boundary.exteriorQ(selected,:),normal,D,gas);
    elseif kind==2
        [flux(faces,:),wallAudit]=shockFittedSlipWallFlux3D(inside,normal,D,gas);wallPressure(selected)=wallAudit.pressurePa;
    elseif kind==3
        shockAudit=movingShockJump3D(boundary.upstream,gas,-normal,-D);
        flux(faces,:)=-shockAudit.downstreamALEFlux;shockPressure(selected)=shockAudit.state.pressure;
        cellState=shockFittedState3D(inside,gas,'conservative');
        shockMismatch(selected)=cellState.pressure./shockAudit.state.pressure-1;
    else
        trace=shockFittedState3D(inside,gas,'conservative');relative=sum(trace.velocity.*normal,2)-D;
        assert(all(relative>trace.soundSpeed),'ShockFit:SubsonicOutlet', ...
            '类型4要求全部法向特征出流；出现亚声速或回流时须另设相容边界，不能给定补压。');
        flux(faces,:)=shockFittedALEFlux3D(inside,inside,normal,D,gas);
        outflowAudit=struct('faceIDs',faces,'relativeNormalMach',relative./trace.soundSpeed, ...
            'pressureGivenAsBoundaryCondition',false,'method','全部法向特征出流/内部迹外推；不指定外压');
    end
end
assert(all(isfinite(flux),'all'),'ShockFit:Flux','存在未计算的共享面通量。');
integrated=geo.averageAreaM2.*flux;divergence=mesh.incidence*integrated;
newQ=(geo.volume0.*Q-dt*divergence)./geo.volume1;
after=shockFittedState3D(newQ,gas,'conservative');
delta=sum(geo.volume1.*newQ-geo.volume0.*Q,1);
boundaryPower=sum(integrated(mesh.boundaryFaces,:),1);
conservation=delta+dt*boundaryPower;
scale=sum(abs(geo.volume0.*Q)+abs(geo.volume1.*newQ),1)+dt*sum(abs(integrated(mesh.boundaryFaces,:)),1);
assert(max(abs(conservation)./max(realmin,scale))<2e-12,'ShockFit:Conservation','全域五量账没有闭合。');
step=struct('Q',newQ,'state',after,'geometry',geo,'faceALEFlux',flux,'boundaryFluxRate',boundaryPower, ...
    'conservedIntegralChange',delta,'conservationResidual',conservation,'conservationScale',scale, ...
    'wallPressurePa',wallPressure,'shockPressurePa',shockPressure,'shockInteriorPressureRelativeMismatch',shockMismatch, ...
    'wallAudit',wallAudit,'shockAudit',shockAudit,'supersonicOutflowAudit',outflowAudit,'budget',budget, ...
    'method','三维四面体共享面ALE五量有限体积/移动系精确Godunov/移动壁面功/同cp激波RH', ...
    'order','空间分片常数、一步显式积分；仅均匀态/GCL和指定基准已可独立核验', ...
    'steadyConverged',false,'shockShapeSolved',false,'bodyMatched',false,'heatFluxSolved',false, ...
    'wallTemperatureSolved',false,'productionQualified',false,'newProductionSolvedAreaM2',0);
end
