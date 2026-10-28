function s=cutCellBoundaryStencilState(Q,cut,gas,boundary)
%CUTCELLBOUNDARYSTENCILSTATE 取得真实边界分支，供曲面重构模板与步长估计共用。
% 模块用途：显式求特征储库的入/出流状态，不能只按来流轴向猜边界。
% 外边界已知样本在扩散模板中冻结；特征状态的非线性敏感性未由此求得。
cv=gas.R/(gas.gamma-1);rho=Q(:,1);u=Q(:,2:4)./rho;
% 与冻结空间内核primitive逐式同序，避免近零特征速度被等价式舍入翻号。
internal=Q(:,5)-.5*rho.*sum(u.^2,2);T=internal./(rho*cv);P=[rho,u,T];
f=cut.openFaces;nf=numel(f.area);names={'xmin','xmax';'ymin','ymax';'zmin','zmax'};
riemannReservoir=buildCutCellRiemannReservoirBoundaryState(Q,cut,gas,boundary);
s=struct('incoming',false(nf,1),'wall',false(nf,1),'primitive',nan(nf,5),'interiorPrimitive',P);
if any(cut.walls.ownerCell>0)&&strcmp(cutCellWallThermalPolicy(boundary.solidWall),'dirichlet-tangent')
    s.thermalDirichletSolidWall=true;
end
for k=find(f.area>0 & (f.cellMinus==0 | f.cellPlus==0))'
    if f.cellMinus(k)>0,id=f.cellMinus(k);normal=f.normal(k,:);side=2;
    else,id=f.cellPlus(k);normal=-f.normal(k,:);side=1;end
    b=boundary.(names{f.axis(k),side});point=f.centroid(k,:);
    for key={'reservoirPrimitive','primitive'}
        if isfield(b,key{1}) && isa(b.(key{1}),'function_handle'),b.(key{1})=b.(key{1})(point,0);end
    end
    if strcmp(b.type,'wall')
        assert(strcmp(cutCellWallThermalPolicy(b),'legacy'),'CutNS:WallThermalBoundary', ...
            '当前显式三维热梯度仅支持原曲面solidWall；外域wall不得静默忽略新策略。');
        given=P(id,:);s.incoming(k)=true;s.wall(k)=true;
        if isa(b.temperature,'function_handle'),given(5)=b.temperature(point,0);else,given(5)=b.temperature;end
        if isa(b.velocity,'function_handle'),given(2:4)=b.velocity(point,0);else,given(2:4)=b.velocity;end
    elseif strcmp(b.type,'characteristic')
        [given,tag]=resolveNSCharacteristicBoundary3D(P(id,:),b,normal,gas);s.incoming(k)=tag.isInflow;
    elseif strcmp(b.type,'riemannReservoir')
        slot=riemannReservoir.slot(k);assert(slot>0,'CutNS:RiemannBoundary','储库模板面没有匹配真实波系。');
        given=riemannReservoir.primitive(slot,:);s.incoming(k)=riemannReservoir.incoming(slot);
    elseif strcmp(b.type,'supersonicInflow')
        given=b.primitive;s.incoming(k)=true;
    elseif strcmp(b.type,'supersonicOutflow')
        given=P(id,:);
    else
        error('CutAdvance:Boundary','不支持或缺失明确边界类型。');
    end
    s.primitive(k,:)=given;
end
end
