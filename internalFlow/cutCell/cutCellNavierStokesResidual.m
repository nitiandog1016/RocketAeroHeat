function [rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary)
%CUTCELLNAVIERSTOKESRESIDUAL 原三角曲面流体控制体的五守恒量空间残差草稿。
% 模块用途：共享开放面只求一次HLLC/Newton/Fourier通量，按真实fluid
% volume更新两侧；原三角wall使用真实流体外法向和给定无滑移壁温。
% Q为全Cartesian单元数×5：rho,rho*u,rho*v,rho*w,rho*E；固体NaN。
% 气体为定gamma/常mu/常k层流近似，SI；不是NASA7/湍流模型。
% 梯度为面积加权最小二乘，加入邻居、墙值及入射边界；通量方向用
% 完整位移约束修正，不把笛卡尔间距或阶梯法向当真实曲面壁距。
% 这里不推进、不处理小单元稳定化、不证明近壁/网格/边界收敛。
% HLLC合同沿用独立库的Davis波速界及正星压力检查，不裁物理状态。
% wall p/q为当前保守态派生；Tw/储库始终为输入，不计物理解覆盖。

%% 1. 明确物性、几何和流体保守态合同；固体不能用零值冒充流体
required={'gamma','R','mu','Pr','temperatureRange'};
assert(all(isfield(gas,required)),'CutNS:Gas','须明确完整常物性气体。');
validateattributes([gas.gamma gas.R gas.mu gas.Pr],{'numeric'},{'real','finite','positive'});
assert(gas.gamma>1,'CutNS:Gas','gamma须大于1。');
validateattributes(gas.temperatureRange,{'numeric'},{'real','finite','positive','numel',2});
assert(gas.temperatureRange(2)>gas.temperatureRange(1),'CutNS:Gas','物性域须递增。');
cp=gas.gamma*gas.R/(gas.gamma-1); cv=gas.R/(gas.gamma-1); conductivity=gas.mu*cp/gas.Pr;
if isfield(gas,'k'),assert(abs(gas.k/conductivity-1)<1e-12,'CutNS:Gas','k与Pr不一致。');end
volume=cut.cells.volume(:); centroid=cut.cells.centroid; n=numel(volume);
assert(isequal(size(Q),[n 5]) && isreal(Q),'CutNS:Shape','Q须为全单元数×5。');
assert(all(isfinite(volume)&volume>=0) && isequal(size(centroid),[n 3]),'CutNS:Geometry','体积或质心合同错误。');
fluid=volume>0;
assert(any(fluid) && all(isnan(Q(~fluid,:)),'all'),'CutNS:Solid','固体保守量须为NaN。');
assert(all(isfinite(centroid(fluid,:)),'all'),'CutNS:Geometry','流体体积质心未知。');
P=nan(n,8); P(fluid,:)=primitive(Q(fluid,:),gas,cv);
phi=P(:,2:5); faces=cut.openFaces; walls=cut.walls;
nf=numel(faces.area); nw=numel(walls.area);
workPolicy=cutCellViscousWorkPolicy(boundary);arithmeticWork=strcmp(workPolicy,'arithmetic');
interiorPolicy=cutCellInteriorEulerPolicy(boundary);exactEuler=strcmp(interiorPolicy,'exact-godunov');
assert(~any(isfield(boundary,{'velocityGradient','velocityTangentialGradient','velocitySurfaceGradient'})), ...
    'CutNS:WallVelocityBoundary','速度梯度策略只能放在原曲面solidWall内，不接受全局字段。');
wallVelocityPolicy='legacy';
if nw>0&&isfield(boundary,'solidWall'),wallVelocityPolicy=cutCellWallVelocityPolicy(boundary.solidWall);end
if arithmeticWork
    workFaceID=find(faces.area>0&faces.cellMinus>0&faces.cellPlus>0);
    oldWorkVelocity=nan(nf,3);newWorkVelocity=nan(nf,3);workTraction=nan(nf,3);
end
assert(all(isfinite(faces.area)&faces.area>=0) && all(isfinite(walls.area)&walls.area>0), ...
    'CutNS:Geometry','面面积无效；不得补负面积或未知孔隙。');
activeWall=walls.ownerCell>0;
assert(all(walls.ownerCell>=0 & walls.ownerCell<=n) && all(fluid(walls.ownerCell(activeWall))), ...
    'CutNS:Owner','墙必须属于真实流体单元。');

%% 2. 先取得显式边界及最小二乘约束，不隐式补入口剖面
A=zeros(3,3,n); B=zeros(3,4,n);
wallInput=nan(nw,4); wallNormal=zeros(nw,3);
faceBoundary=cell(nf,1); facePrimitive=nan(nf,8); faceIsInflow=false(nf,1);
boundaryNames={'xmin','xmax';'ymin','ymax';'zmin','zmax'};
% 显式新储库按实际Riemann波系批量求迹；没有此type时返回空且保留旧路径。
riemannReservoir=buildCutCellRiemannReservoirBoundaryState(Q,cut,gas,boundary);
areaClosure=zeros(n,3); momentClosure=zeros(3,3,n);
for f=1:nf
    area=faces.area(f); if area==0,continue;end
    m=faces.cellMinus(f); p=faces.cellPlus(f); normal=faces.normal(f,:); point=faces.centroid(f,:);
    assert(isfinite(area) && all(isfinite(point)) && abs(norm(normal)-1)<1e-12,'CutNS:Geometry','开放面几何无效。');
    assert((m==0 || (m>=1&&m<=n&&fluid(m))) && (p==0 || (p>=1&&p<=n&&fluid(p))) && (m~=0||p~=0), ...
        'CutNS:Owner','非零开放面连到固体或非法单元。');
    if m>0
        areaClosure(m,:)=areaClosure(m,:)+area*normal;
        momentClosure(:,:,m)=momentClosure(:,:,m)+(point-centroid(m,:))'*(area*normal);
    end
    if p>0
        areaClosure(p,:)=areaClosure(p,:)-area*normal;
        momentClosure(:,:,p)=momentClosure(:,:,p)-(point-centroid(p,:))'*(area*normal);
    end
    if m>0 && p>0
        delta=centroid(p,:)-centroid(m,:); difference=phi(p,:)-phi(m,:);
        [A(:,:,m),B(:,:,m)]=constraint(A(:,:,m),B(:,:,m),delta,difference,area);
        [A(:,:,p),B(:,:,p)]=constraint(A(:,:,p),B(:,:,p),-delta,-difference,area);
    else
        if m>0,id=m; outward=normal; side=2;else,id=p;outward=-normal;side=1;end
        axis=faces.axis(f); assert(axis>=1 && axis<=3,'CutNS:Geometry','开放面的axis无效。');
        assert(isfield(boundary,boundaryNames{axis,side}),'CutNS:Boundary','缺少显式外域边界。');
        b=boundary.(boundaryNames{axis,side});
        assert(strcmp(cutCellWallVelocityPolicy(b),'legacy'),'CutNS:WallVelocityBoundary', ...
            '当前显式三维墙速度梯度仅支持原曲面solidWall，不接受开放面边界策略。');
        if isfield(b,'reservoirPrimitive') && isa(b.reservoirPrimitive,'function_handle')
            % 明确给定的空间变化储库，用于入口/解析验证，不是从壁面外缘猜全截面。
            b.reservoirPrimitive=b.reservoirPrimitive(point,0);
        end
        if isfield(b,'primitive') && isa(b.primitive,'function_handle'),b.primitive=b.primitive(point,0);end
        if strcmp(b.type,'wall')
            input=wallValues(b,point,0,outward,gas); given=P(id,:); given(2:5)=input;
            incoming=true;
        elseif strcmp(b.type,'characteristic')
            [given5,characteristic]=resolveNSCharacteristicBoundary3D(P(id,1:5),b,outward,gas);
            given=primitive(conservative(given5,cv),gas,cv); incoming=characteristic.isInflow;
        elseif strcmp(b.type,'riemannReservoir')
            slot=riemannReservoir.slot(f);assert(slot>0,'CutNS:RiemannBoundary','真实储库面未取得批量波系状态。');
            given=primitive(conservative(riemannReservoir.primitive(slot,:),cv),gas,cv);incoming=riemannReservoir.incoming(slot);
        elseif strcmp(b.type,'supersonicInflow')
            given=primitive(conservative(b.primitive,cv),gas,cv);
            assert(dot(given(2:4),outward)<-given(7),'CutNS:BoundaryMach','输入并非向内超声速。');incoming=true;
        elseif strcmp(b.type,'supersonicOutflow')
            given=P(id,:);assert(dot(given(2:4),outward)>given(7),'CutNS:BoundaryMach','出口并非全特征向外。');incoming=false;
        else
            error('CutNS:Boundary','不支持或未明确的边界类型；不自动生成周期配对。');
        end
        faceBoundary{f}=b; facePrimitive(f,:)=given; faceIsInflow(f)=incoming;
        if incoming
            [A(:,:,id),B(:,:,id)]=constraint(A(:,:,id),B(:,:,id),point-centroid(id,:),given(2:5)-phi(id,:),area);
        end
    end
end
assert(nw==0 || isfield(boundary,'solidWall'),'CutNS:Boundary','曲面墙须明确solidWall条件。');
for w=1:nw
    if ~activeWall(w),continue;end % owner0为计算域外侧实体边界，没有域内fluid-wall通量。
    id=walls.ownerCell(w); point=walls.centroid(w,:); area=walls.area(w);
    if isfield(walls,'normal'),normal=walls.normal(w,:);else,normal=walls.areaVector(w,:)/area;end
    assert(abs(norm(normal)-1)<1e-10,'CutNS:Geometry','真实墙面积向量与面积不一致。');
    assert(norm(walls.areaVector(w,:)-area*normal)<=128*eps(max(area,realmin)), ...
        'CutNS:Geometry','墙面积向量须与可靠的原三角单位法向一致。');
    wallNormal(w,:)=normal;
    input=wallValues(boundary.solidWall,point,walls.originalFaceID(w),normal,gas);wallInput(w,:)=input;
    [A(:,:,id),B(:,:,id)]=constraint(A(:,:,id),B(:,:,id),point-centroid(id,:),input-phi(id,:),area);
    areaClosure(id,:)=areaClosure(id,:)+walls.areaVector(w,:);
    momentClosure(:,:,id)=momentClosure(:,:,id)+(point-centroid(id,:))'*walls.areaVector(w,:);
end
gradient=nan(4,3,n); geometryError=0; matrixError=0;
for id=find(fluid)'
    areaScale=volume(id)^(2/3);
    geometryError=max(geometryError,norm(areaClosure(id,:))/areaScale);
    matrixError=max(matrixError,norm(momentClosure(:,:,id)-volume(id)*eye(3),'fro')/volume(id));
    assert(rcond(A(:,:,id))>1e-13,'CutNS:GradientRank','流体单元%d梯度约束退化，须改善几何/模板。',id);
    gradient(:,:,id)=(A(:,:,id)\B(:,:,id))';
end
assert(geometryError<1e-7 && matrixError<1e-7,'CutNS:GeometryClosure','逐流体单元几何未闭合，拒绝通量。');

% 显式精确通量按真实内部面批量求解；默认HLLC不经过新数值路径。
exactConvection=[];exactAudit=[];
if exactEuler
    exactIDs=find(faces.area>0&faces.cellMinus>0&faces.cellPlus>0);exactConvection=nan(nf,5);
    if ~isempty(exactIDs)
        [value,exactAudit]=cutCellExactEulerFlux(Q(faces.cellMinus(exactIDs),:),Q(faces.cellPlus(exactIDs),:),faces.normal(exactIDs,:),gas);
        exactConvection(exactIDs,:)=value;exactAudit.faceID=exactIDs;
    else
        exactAudit=struct('kind','explicit-exact-Godunov-no-active-interior-faces','faceID',zeros(0,1));
    end
end

%% 3. 每个开放面共享同一对流/黏性/导热通量；正反号严格成对
integratedRate=zeros(n,5); integratedRate(~fluid,:)=NaN;
faceFlux=nan(nf,5);convectiveFlux=nan(nf,5);viscousFaceFlux=nan(nf,5);
netBoundary=zeros(1,5); openBoundary=zeros(1,5); outerWall=zeros(1,5);
for f=1:nf
    area=faces.area(f);if area==0,continue;end
    m=faces.cellMinus(f); p=faces.cellPlus(f); normal=faces.normal(f,:); point=faces.centroid(f,:);
    if m>0 && p>0
        delta=centroid(p,:)-centroid(m,:);
        G=directionalCorrection(.5*(gradient(:,:,m)+gradient(:,:,p)),delta,phi(p,:)-phi(m,:));
        um=P(m,2:4)+(gradient(1:3,:,m)*(point-centroid(m,:))')';
        up=P(p,2:4)+(gradient(1:3,:,p)*(point-centroid(p,:))')';
        if exactEuler,convection=exactConvection(f,:);
        else,convection=hllc(Q(m,:),Q(p,:),P(m,:),P(p,:),normal,gas,cv);end
        if arithmeticWork
            oldVelocity=.5*(um+up);selectedVelocity=.5*(P(m,2:4)+P(p,2:4));
            diffusion=viscousFlux(G,selectedVelocity,normal,gas.mu,conductivity);
            oldWorkVelocity(f,:)=oldVelocity;newWorkVelocity(f,:)=selectedVelocity;workTraction(f,:)=diffusion(2:4);
        else
            diffusion=viscousFlux(G,.5*(um+up),normal,gas.mu,conductivity);
        end
        F=convection-diffusion;
    else
        if m>0,id=m;outward=normal;sign=1;else,id=p;outward=-normal;sign=-1;end
        b=faceBoundary{f}; given=facePrimitive(f,:); G=gradient(:,:,id);
        if strcmp(b.type,'wall')
            assert(strcmp(cutCellWallThermalPolicy(b),'legacy'),'CutNS:WallThermalBoundary', ...
                '当前显式三维热梯度仅支持原曲面solidWall。');
            G=directionalCorrection(G,point-centroid(id,:),given(2:5)-phi(id,:));
            convection=sign*[0,P(id,6)*outward,0];
            diffusion=sign*viscousFlux(G,given(2:4),outward,gas.mu,conductivity);
            outwardFlux=sign*(convection-diffusion);
            outerWall=outerWall+area*outwardFlux;
        else
            if faceIsInflow(f)
                G=directionalCorrection(G,point-centroid(id,:),given(2:5)-phi(id,:));
            else
                % 显式零法向黏性梯度出口闭合，不冒充完整无反射边界。
                G=G-(G*outward')*outward;
            end
            convection=sign*physicalFlux(conservative(given(1:5),cv),given,outward);
            if strcmp(b.type,'riemannReservoir'),convection=sign*riemannReservoir.outwardConvectiveFlux(riemannReservoir.slot(f),:);end
            diffusion=sign*viscousFlux(G,given(2:4),outward,gas.mu,conductivity);
            outwardFlux=sign*(convection-diffusion);
            openBoundary=openBoundary+area*outwardFlux;
        end
        F=sign*outwardFlux;netBoundary=netBoundary+area*outwardFlux;
    end
    faceFlux(f,:)=F;convectiveFlux(f,:)=convection;viscousFaceFlux(f,:)=diffusion;
    if m>0,integratedRate(m,:)=integratedRate(m,:)-area*F;end
    if p>0,integratedRate(p,:)=integratedRate(p,:)+area*F;end
end

%% 4. 原三角曲面wall通量，法向为fluid向solid；给定Tw与派生p/q分开
wallFlux=nan(nw,5);wallHeat=nan(nw,1);wallPressure=nan(nw,1);wallPower=nan(nw,1);
explicitWallGradient=[];wallThermalAudit=[];
explicitVelocityGradient=[];wallVelocityAudit=[];wallVelocityIndex=zeros(nw,1);
if any(activeWall)&&strcmp(cutCellWallThermalPolicy(boundary.solidWall),'dirichlet-tangent')
    ids=find(activeWall);owner=walls.ownerCell(ids);
    [selectedGradient,wallThermalAudit]=cutCellWallThermalGradient(centroid(owner,:),walls.centroid(ids,:), ...
        wallNormal(ids,:),P(owner,5),wallInput(ids,4),boundary.solidWall,ids);
    explicitWallGradient=nan(nw,3);explicitWallGradient(ids,:)=selectedGradient;
end
if any(activeWall)&&strcmp(wallVelocityPolicy,'dirichlet-tangent')
    ids=find(activeWall);owner=walls.ownerCell(ids);
    [selectedVelocityGradient,wallVelocityAudit]=cutCellWallVelocityGradient(centroid(owner,:),walls.centroid(ids,:), ...
        wallNormal(ids,:),P(owner,2:4),wallInput(ids,1:3),boundary.solidWall,ids);
    explicitVelocityGradient=nan(3,3,nw);explicitVelocityGradient(:,:,ids)=selectedVelocityGradient;
    wallVelocityIndex(ids)=1:numel(ids);wallVelocityAudit.oldGradientXYZ=nan(3,3,numel(ids));
end
wallNet=zeros(1,5);
for w=1:nw
    if ~activeWall(w),continue;end
    id=walls.ownerCell(w);normal=wallNormal(w,:);input=wallInput(w,:);
    G=directionalCorrection(gradient(:,:,id),walls.centroid(w,:)-centroid(id,:),input-phi(id,:));
    if ~isempty(explicitWallGradient),G(4,:)=explicitWallGradient(w,:);end
    if ~isempty(explicitVelocityGradient)
        wallVelocityAudit.oldGradientXYZ(:,:,wallVelocityIndex(w))=G(1:3,:);
        G(1:3,:)=explicitVelocityGradient(:,:,w);
    end
    visc=viscousFlux(G,input(1:3),normal,gas.mu,conductivity);
    F=[0,P(id,6)*normal,0]-visc;
    integratedRate(id,:)=integratedRate(id,:)-walls.area(w)*F;
    wallFlux(w,:)=F;wallPressure(w)=P(id,6);wallHeat(w)=-conductivity*dot(G(4,:),normal);
    wallPower(w)=dot(input(1:3),visc(2:4));
    wallNet=wallNet+walls.area(w)*F;
end
netBoundary=netBoundary+wallNet;
rate=nan(n,5);rate(fluid,:)=integratedRate(fluid,:)./volume(fluid);
assert(all(isfinite(rate(fluid,:)),'all') && all(isfinite(faceFlux(faces.area>0,:)),'all') && ...
    all(isfinite(wallFlux(activeWall,:)),'all'),'CutNS:NonfiniteResidual', ...
    '有限保守态不能保证极小体积上的残差可表示；拒绝溢出，不裁剪残差或删除小单元。');
volumeIntegrated=sum(integratedRate(fluid,:),1);
info=struct('status','experimental-curved-cut-cell-spatial-only','physicalSolution',false, ...
    'eligibleForSurfaceCoverage',false,'steadyConverged',false,'timeIntegratorImplemented',false, ...
    'scope','常gamma/常mu/k层流空间残差；面积加权最小二乘及位移约束梯度；未处理小单元稳定推进/湍流/物性/网格域收敛', ...
    'method','HLLC共享开放面 + Newton应力 + Fourier原三角曲面热流', ...
    'primitive',P(:,1:5),'gradient',gradient,'faceFlux',faceFlux,'convectiveFaceFlux',convectiveFlux, ...
    'viscousFaceFlux',viscousFaceFlux,'wallFlux',wallFlux, ...
    'wallPressure',wallPressure,'wallHeatFluxIntoWall',wallHeat, ...
    'wallMechanicalPowerIntoFluid',wallPower,'wallTemperatureInput',wallInput(:,4), ...
    'wallActiveInFluidDomain',activeWall, ...
    'wallTemperatureValueKind','prescribed','wallPressureValueKind','derived-current-fluid-state', ...
    'wallHeatFluxValueKind','derived-current-fluid-state-not-converged', ...
    'netBoundaryFlux',netBoundary,'openBoundaryFlux',openBoundary,'outerWallBoundaryFlux',outerWall, ...
    'curvedWallBoundaryFlux',wallNet,'volumeIntegratedResidual',volumeIntegrated, ...
    'balanceResidual',volumeIntegrated+netBoundary,'geometryVectorError',geometryError, ...
    'geometryTensorError',matrixError,'positivityClippingCount',0,'stabilityBoundVerified',false);
if ~isempty(riemannReservoir.faceID)
    info.riemannReservoirBoundaryAudit=rmfield(riemannReservoir,'slot');
    info.riemannReservoirBoundaryAudit.method='显式三维储库 / 精确Godunov xi=0；入流按采样质量方向，出流零法向黏性梯度';
    info.riemannReservoirBoundaryAudit.physicalValidation=false;
    info.method=[info.method '；开放储库采用精确Riemann波系，物理方程/边界身份与旧特征组分开'];
end
if ~isempty(wallThermalAudit)
    info.wallThermalGradientAudit=wallThermalAudit;
    info.method=[info.method '；显式三维Dirichlet壁温切向约束'];
end
if ~isempty(wallVelocityAudit)
    info.wallVelocityGradientAudit=wallVelocityAudit;
    info.method=[info.method '；显式三维静止壁速度切向约束'];
end
if arithmeticWork
    ids=workFaceID;m=faces.cellMinus(ids);p=faces.cellPlus(ids);old=oldWorkVelocity(ids,:);selected=newWorkVelocity(ids,:);traction=workTraction(ids,:);
    oldPower=sum(old.*traction,2);selectedPower=sum(selected.*traction,2);pair=sum((P(p,2:4)-P(m,2:4)).*traction,2);
    ownerHalf=sum((selected-P(m,2:4)).*traction,2);neighbourHalf=sum((P(p,2:4)-selected).*traction,2);
    info.viscousWorkAudit=struct('kind','explicit-arithmetic-three-dimensional-interior-viscous-work', ...
        'policy',workPolicy,'faceID',ids,'oldVelocityXYZ',old,'velocityXYZ',selected,'tractionXYZ',traction, ...
        'oldWorkFluxPerArea',oldPower,'workFluxPerArea',selectedPower,'pairMomentumToInternalEnergyPerArea',pair, ...
        'ownerHalfPairDifference',ownerHalf-.5*pair,'neighbourHalfPairDifference',neighbourHalf-.5*pair, ...
        'negativePairCount',nnz(pair<0),'boundaryVelocityChanged',false,'wallHeatChanged',false, ...
        'tractionChanged',false,'cellStateClipped',false,'generalSkewGridPositivityProved',false, ...
        'cellVelocityMeaning','conservative-mean-derived point proxy, not exact FV mean velocity', ...
        'secondOrderLinearReconstructionClaimed',false,'physicalValidation',false);
    info.method=[info.method '；内部黏性功显式一阶算术速度代理（非限幅原迹）'];
end
if exactEuler
    info.exactEulerFluxAudit=exactAudit;
    info.method=strrep(info.method,'HLLC共享开放面','精确理想气体Riemann-Godunov共享内部面（边界迹保持）');
end
end

function [A,B]=constraint(A,B,delta,difference,area)
% 模块用途：以面积/位移平方加权，线性样本可由完整三方向矩阵恢复。
d2=dot(delta,delta);
assert(isfinite(d2)&&d2>0,'CutNS:Distance','约束采样点与体积质心重合，不能构造壁距。');
w=area/d2;A=A+w*(delta'*delta);B=B+w*(delta'*difference);
end

function G=directionalCorrection(G,delta,difference)
% 模块用途：修正完整位移方向导数，保留非正交曲面及cell质心位置。
% 不以可能为负的单一法向投影距离作分母；仍需近壁/网格适用性验证。
d2=dot(delta,delta);assert(isfinite(d2)&&d2>0,'CutNS:Distance','通量采样位移为零。');
G=G+((difference(:)-G*delta')/d2)*delta;
end

function input=wallValues(b,point,faceID,normal,gas)
assert(strcmp(b.type,'wall'),'CutNS:Boundary','solidWall必须明确为wall。');
if isa(b.temperature,'function_handle'),T=b.temperature(point,faceID);else,T=b.temperature;end
if isa(b.velocity,'function_handle'),u=b.velocity(point,faceID);else,u=b.velocity;end
validateattributes(u,{'numeric'},{'real','finite','numel',3});u=double(u(:)');
assert(isscalar(T)&&isfinite(T)&&isreal(T)&&T>=gas.temperatureRange(1)&&T<=gas.temperatureRange(2), ...
    'CutNS:Temperature','给定墙温越过明确物性域。');
assert(abs(dot(u,normal))<1e-11*max(1,norm(u)),'CutNS:WallMotion','固定几何墙法向速度必须零。');
input=[u,T];
end

function P=primitive(Q,gas,cv)
assert(isreal(Q)&&all(isfinite(Q),'all')&&all(Q(:,1)>0),'CutNS:State','流体密度/保守量无效，不裁剪。');
rho=Q(:,1);u=Q(:,2:4)./rho;internal=Q(:,5)-.5*rho.*sum(u.^2,2);
assert(all(internal>0),'CutNS:InternalEnergy','内能非正，不裁剪。');
T=internal./(rho*cv);assert(all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)), ...
    'CutNS:Temperature','温度越过明确常物性模型域，不裁温。');
pressure=rho*gas.R.*T;P=[rho,u,T,pressure,sqrt(gas.gamma*gas.R*T),internal];
assert(all(isfinite(P),'all'),'CutNS:State','气体派生状态非有限，不外推或裁剪。');
end

function Q=conservative(P,cv)
Q=[P(1),P(1)*P(2:4),P(1)*(cv*P(5)+.5*sum(P(2:4).^2))];
end

function F=physicalFlux(Q,P,normal)
un=dot(P(2:4),normal);F=[Q(1)*un,Q(2:4)*un+P(6)*normal,(Q(5)+P(6))*un];
end

function F=viscousFlux(G,u,normal,mu,k)
velocityGradient=G(1:3,:);stress=mu*(velocityGradient+velocityGradient'-(2/3)*trace(velocityGradient)*eye(3));
traction=(stress*normal')';F=[0,traction,dot(u,traction)+k*dot(G(4,:),normal)];
end

function F=hllc(QL,QR,PL,PR,normal,gas,cv)
% 模块用途：三波HLLC，正星压力与实际星内能分别核查，无HLL回退。
uL=dot(PL(2:4),normal);uR=dot(PR(2:4),normal);
SL=min(uL-PL(7),uR-PR(7));SR=max(uL+PL(7),uR+PR(7));
FL=physicalFlux(QL,PL,normal);FR=physicalFlux(QR,PR,normal);
if SL>=0,F=FL;return;elseif SR<=0,F=FR;return;end
denominator=PL(1)*(SL-uL)-PR(1)*(SR-uR);
assert(isfinite(denominator)&&denominator~=0,'CutNS:HLLC','HLLC分母退化。');
SM=(PR(6)-PL(6)+PL(1)*uL*(SL-uL)-PR(1)*uR*(SR-uR))/denominator;
assert(all(isfinite([SL SM SR]))&&SM>SL&&SM<SR,'CutNS:HLLC','HLLC波序无效。');
jumpL=PL(1)*(SL-uL)*(SM-uL);jumpR=PR(1)*(SR-uR)*(SM-uR);
pL=PL(6)+jumpL;pR=PR(6)+jumpR;
assert(all(isfinite([pL pR]))&&pL>0&&pR>0,'CutNS:HLLCStarPressure','星压力无效，不裁剪或自动回退。');
scale=max([1 abs(PL(6)) abs(PR(6)) abs(jumpL) abs(jumpR)]);
assert(abs(pL-pR)<=128*eps(scale),'CutNS:HLLCStarPressureMismatch','左右星压力不满足舍入一致性。');
if SM>=0,QS=star(QL,PL,SL,SM,uL,normal);primitive(QS,gas,cv);F=FL+SL*(QS-QL);
else,QS=star(QR,PR,SR,SM,uR,normal);primitive(QS,gas,cv);F=FR+SR*(QS-QR);end
end

function QS=star(Q,P,S,SM,un,normal)
rho=P(1)*(S-un)/(S-SM);u=P(2:4)+(SM-un)*normal;
energy=Q(5)/P(1)+(SM-un)*(SM+P(6)/(P(1)*(S-un)));
QS=[rho,rho*u,rho*energy];
end
