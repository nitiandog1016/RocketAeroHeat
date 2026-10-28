function out=evaluateQuadraticTemperature(meanField,cut,boundary,stencil)
%EVALUATEQUADRATICTEMPERATURE 用明确声明的温度FV平均求局部梯度/Hessian。
% 模块用途：一般可压缩Q的EOS primitive不是本接口可默认接受的温度均值。
% 全部边界点温度标输入；局部质心多项式值与FV平均分开保存。
assert(isstruct(meanField)&&isfield(meanField,'valueKind')&&strcmp(meanField.valueKind,'volume-average-temperature') ...
    &&isfield(meanField,'units')&&strcmp(meanField.units,'K'),'CutQuadratic:MeanIdentity', ...
    '必须明确提供volume-average-temperature [K]，不能默认为EOS primitive点样本。');
T=meanField.values(:);n=numel(cut.cells.volume);fluid=cut.cells.volume>0;
assert(isnumeric(T)&&isreal(T)&&numel(T)==n&&all(isfinite(T(fluid))&T(fluid)>0)&&all(isnan(T(~fluid))), ...
    'CutQuadratic:MeanValues','真实温度均值非法或固体未保持NaN。');
assert(isequaln(stencil.volume,cut.cells.volume(:))&&isequaln(stencil.centroid,cut.cells.centroid), ...
    'CutQuadratic:GeometryIdentity','模板不是当前CV几何。');
assert(isequaln(stencil.geometryOpenFaces,cut.openFaces)&&isequaln(stencil.geometryWalls,cut.walls), ...
    'CutQuadratic:GeometryIdentity','共享面/墙几何或所有权已改变，须重建模板。');
s=stencil.samples;given=zeros(numel(s.owner),1);wallInput=nan(numel(cut.walls.area),1);
wallInputKind=repmat({'unknown-wall-temperature'},numel(wallInput),1);
for k=1:numel(given)
    bc=boundary.(s.boundaryName{k});
    if s.kind(k)==1
        assert(strcmp(bc.type,'temperaturePoint'),'CutQuadratic:BoundaryIdentity','边界类型改变，须重建模板。');
        id=0;if s.isCurvedWall(k),id=cut.walls.originalFaceID(s.sourceID(k));end
        if isa(bc.value,'function_handle'),value=bc.value(s.point(k,:),id);else,value=bc.value;end
        assert(isscalar(value)&&isreal(value)&&isfinite(value)&&value>0,'CutQuadratic:BoundaryValue','给定点温度非法。');given(k)=value;
        if s.isCurvedWall(k)
            wallInput(s.sourceID(k))=value;wallInputKind{s.sourceID(k)}='prescribed-point-boundary-input';
        end
    else,assert(strcmp(bc.type,'zeroNormalGradient'),'CutQuadratic:BoundaryIdentity','边界类型改变，须重建模板。');end
end
gradient=nan(n,3);hessian=nan(n,3,3);pointValue=nan(n,1);fitResidual=nan(n,1);
for id=find(fluid)'
    local=stencil.local{id};bids=local.boundarySampleIDs;rhs=[T(local.neighbors)-T(id);given(bids)];
    point=s.kind(bids)==1;rhs(numel(local.neighbors)+find(point))=given(bids(point))-T(id);
    theta=local.inverse*rhs;L=local.length;gradient(id,:)=theta(1:3)'/L;
    H=[theta(4),theta(7),theta(8);theta(7),theta(5),theta(9);theta(8),theta(9),theta(6)]/L^2;
    hessian(id,:,:)=reshape(H,1,3,3);Cov=reshape(stencil.covariance(id,:,:),3,3);
    pointValue(id)=T(id)-.5*sum(H.*Cov,'all');fitResidual(id)=norm(local.weight.*(local.design*theta-rhs))/max(1,norm(local.weight.*rhs));
end
assert(all(isfinite(gradient(fluid,:)),'all')&&all(isfinite(hessian(fluid,:,:)),'all'),'CutQuadratic:Nonfinite','重构非有限，不裁剪系数。');
out=struct('temperatureMean',T,'cellPolynomialValueAtCentroid',pointValue,'gradient',gradient, ...
    'hessian',hessian,'fitResidual',fitResidual,'wallTemperatureInput',wallInput, ...
    'wallTemperatureValueKind',{wallInputKind}, ...
    'method',stencil.method,'physicalSolution',false,'eligibleForSurfaceCoverage',false);
end
