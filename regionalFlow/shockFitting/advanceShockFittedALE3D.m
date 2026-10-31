function step=advanceShockFittedALE3D(Q,mesh,nodes0,nodes1,dt,boundary,gas,maximumArrayBytes)
%ADVANCESHOCKFITTEDALE3D 共享面三维五量ALE有限体积一步，含移动壁面和显式激波边界。
% 模块一：内部/储库面用移动系精确Riemann，壁压由镜像Riemann求出；类型4仅允许超声速外流。
% 模块二：激波3使用原RH通量；候选边界5在实际网格速度处取上/下游精确Riemann数值通量。
% 模块五：类型5保留与实际移动RH的五量通量缺陷及原单元压力差，不能据数值通量就授予激波资格。
% 模块三：以新旧真实体积更新VQ，核对全域五量账；非正态明确拒绝。
% 模块四：声明缓存和Riemann临时工作数组预算；它不是MATLAB进程峰值测量。
if nargin<8,maximumArrayBytes=128*1024^2;end
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:FluxGas','当前一步内核仅常gamma，不回退NASA7。');
assert(isequal(size(Q),[mesh.cellCount 5]),'ShockFit:State','单元五量数目不符。');
validateattributes(maximumArrayBytes,{'double'},{'real','finite','positive','scalar'});
Nb=numel(mesh.boundaryFaces);
assert(isstruct(boundary)&&all(isfield(boundary,{'kind','exteriorQ','upstream'}))&& ...
    isequal(size(boundary.kind),[Nb 1])&&isequal(size(boundary.exteriorQ),[Nb 5])&& ...
    all(ismember(boundary.kind,[1 2 3 4 5])),'ShockFit:Boundary','边界须显式为储库1/壁面2/RH激波3/超声速外流4/候选Riemann激波5。');
inputs=whos('Q','mesh','nodes0','nodes1','boundary','gas');inputBytes=sum([inputs.bytes]);
% exactEuler最多预留1000个每面double标量，包含左右态/根/波扇/审计与表达式临时副本。
% 几何预留每单元200/每面120/每节点24个double，另加输出与更新五量副本。
workspace=8*(1000*mesh.faceCount+200*mesh.cellCount+120*mesh.faceCount+24*mesh.nodeCount);
hasRiemannShock=any(boundary.kind==5);
if hasRiemannShock,workspace=workspace+8*80*Nb;end
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
shockFluxDefect=zeros(0,5);shockStarPressure=zeros(0,1);shockNumericalAudit=struct;
if hasRiemannShock,shockFluxDefect=nan(Nb,5);shockStarPressure=nan(Nb,1);end
for kind=1:5
    selected=find(boundary.kind==kind);if isempty(selected),continue;end
    faces=mesh.boundaryFaces(selected);inside=Q(mesh.owner(faces),:);normal=geo.normal(faces,:);D=geo.gridNormalSpeedMS(faces);
    if kind==1
        flux(faces,:)=shockFittedALEFlux3D(inside,boundary.exteriorQ(selected,:),normal,D,gas);
    elseif kind==2
        [flux(faces,:),wallAudit]=shockFittedSlipWallFlux3D(inside,normal,D,gas);wallPressure(selected)=wallAudit.pressurePa;
    elseif kind==3||kind==5
        shockAudit=movingShockJump3D(boundary.upstream,gas,-normal,-D);
        referenceFlux=-shockAudit.downstreamALEFlux;shockPressure(selected)=shockAudit.state.pressure;
        if kind==3
            flux(faces,:)=referenceFlux;
        else
            up=shockFittedState3D([boundary.upstream.density boundary.upstream.velocityVector' ...
                boundary.upstream.temperature],gas,'primitive');
            [flux(faces,:),shockNumericalAudit]=shockFittedALEFlux3D(inside,repmat(up.Q,numel(faces),1),normal,D,gas);
            assert(all(shockNumericalAudit.riemann.pressureStarPa>boundary.upstream.pressure), ...
                'ShockFit:ShockBoundaryWave','候选边界右波不是压缩激波，不能授予拟合激波身份。');
            speedScale=max(norm(boundary.upstream.velocityVector),up.soundSpeed);
            physicalScale=boundary.upstream.density*[speedScale repmat(speedScale^2,1,3) speedScale^3];
            shockFluxDefect(selected,:)=(flux(faces,:)-referenceFlux)./physicalScale;
            shockStarPressure(selected)=shockNumericalAudit.riemann.pressureStarPa;
        end
        cellState=shockFittedState3D(inside,gas,'conservative');
        shockMismatch(selected)=cellState.pressure./shockAudit.state.pressure-1;
    elseif kind==4
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
if hasRiemannShock
    step.shockBoundaryFluxRelativeMismatch=shockFluxDefect;step.shockRiemannStarPressurePa=shockStarPressure;
    step.shockNumericalAudit=shockNumericalAudit;
    step.method='三维共享面ALE/候选激波实际移动面上的上游-内部精确Riemann通量/独立RH通量缺陷';
end
end
