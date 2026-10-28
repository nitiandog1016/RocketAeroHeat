function value=applyRocketRectangularCoupling(coupling,rhs)
%APPLYROCKETRECTANGULARCOUPLING 使用全部原五量矩形耦合块，不按阈值删除连接。
% 模块一：适用于原F-C、C-F或C-C子块，行列控制体编号分开保留。
% 模块二：批量块乘积后按原所有者求和，不调用原NS。
assert(strcmp(coupling.kind,'all-original-five-field-rectangular-coupling')&& ...
    isequal(size(rhs),[coupling.columnCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'SpatialCoupling:Input','矩形耦合或五量右端不符。');
value=zeros(coupling.rowCount,5);if isempty(coupling.column),return;end
ordered=rhs(coupling.column,:)';product=pagemtimes(coupling.values,reshape(ordered,5,1,[]));
for equation=1:5,value(:,equation)=accumarray(double(coupling.owner),reshape(product(equation,1,:),[],1),[coupling.rowCount 1]);end
assert(all(isfinite(value),'all'),'SpatialCoupling:Nonfinite','完整矩形块作用非有限。');
end
