function out=cutCellBoundaryDirectionalState(Q,directionU,cut,gas,boundary,scale,base)
%CUTCELLBOUNDARYDIRECTIONALSTATE 传播全部真实开放面的局部特征状态导数。
% 模块一：守恒扰动转换为rho/u/v/w/T扰动，不改变给定墙速/壁温/储库。
% 模块二：逐真实面使用原axis及外法向，入/出流分支必须与基态一致。
% 模块三：输出供WLS约束及物理边界通量共同使用，不把边界输入标为区域解。
scale=reshape(scale,1,5);dq=directionU.*scale;rho=Q(:,1);u=Q(:,2:4)./rho;cv=gas.R/(gas.gamma-1);
T=(Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);du=(dq(:,2:4)-u.*dq(:,1))./rho;
dT=(dq(:,5)-sum(u.*dq(:,2:4),2)+(.5*sum(u.^2,2)-cv*T).*dq(:,1))./(rho*cv);
dP=[dq(:,1),du,dT];assert(isequal(size(Q),size(directionU))&&all(isfinite(dP),'all'));
f=cut.openFaces;names={'xmin','xmax';'ymin','ymax';'zmin','zmax'};given=zeros(numel(f.area),5);
counts=struct('supersonicOutflow',0,'supersonicInflow',0,'subsonicOutflow',0,'reservoirBackflow',0,'wall',0);
riemann=buildCutCellRiemannReservoirJacobian(Q,cut,gas,boundary);
if ~isempty(riemann.faceID),counts.riemannReservoir=numel(riemann.faceID);end
for k=find(f.area>0&(f.cellMinus==0|f.cellPlus==0))'
    if f.cellMinus(k)>0,id=f.cellMinus(k);n=f.normal(k,:);side=2;else,id=f.cellPlus(k);n=-f.normal(k,:);side=1;end
    b=boundary.(names{f.axis(k),side});point=f.centroid(k,:);
    for key={'reservoirPrimitive','primitive'}
        if isfield(b,key{1})&&isa(b.(key{1}),'function_handle'),b.(key{1})=b.(key{1})(point,0);end
    end
    if strcmp(b.type,'characteristic')
        [value,tag]=cutCellCharacteristicBoundaryDirection(base.interiorPrimitive(id,:),dP(id,:),b,n,gas);
        assert(isequaln(tag.basePrimitive,base.primitive(k,:))&&tag.isInflow==base.incoming(k), ...
            'CutDirection:BoundaryBase','特征基态或分支没有复现原模板状态。');
        given(k,:)=value;
        switch tag.branch
            case 'supersonic-outflow',counts.supersonicOutflow=counts.supersonicOutflow+1;
            case 'supersonic-inflow',counts.supersonicInflow=counts.supersonicInflow+1;
            case 'subsonic-outflow',counts.subsonicOutflow=counts.subsonicOutflow+1;
            case 'reservoir-backflow',counts.reservoirBackflow=counts.reservoirBackflow+1;
        end
    elseif strcmp(b.type,'riemannReservoir')
        slot=riemann.slot(k);assert(slot>0&&riemann.incoming(slot)==base.incoming(k));
        given(k,:)=reshape(riemann.primitiveJacobian(:,:,slot)*dP(id,:)',1,5);
    elseif strcmp(b.type,'supersonicInflow'),counts.supersonicInflow=counts.supersonicInflow+1;
    elseif strcmp(b.type,'supersonicOutflow'),given(k,:)=dP(id,:);counts.supersonicOutflow=counts.supersonicOutflow+1;
    elseif strcmp(b.type,'wall'),given(k,1)=dP(id,1);counts.wall=counts.wall+1;
    else,error('CutDirection:Boundary','未支持的明确边界类型。');end
end
out=struct('primitiveDerivative',given,'interiorPrimitiveDerivative',dP,'incoming',base.incoming,'wall',base.wall, ...
    'branchCounts',counts,'boundarySensitivityIncluded',true,'selectedBranchOnly',true,'givenInputsChanged',false);
end
