function coordinates=shockSurfaceMassSource(map,eta)
%SHOCKSURFACEMASSSOURCE 用正向投影质量比寻找同一共同激波线的入口点。
% 模块一：eta=1为当前激波点；eta=0只返回驻点熵来源极限，不是驻点流场解。
% 模块二：对数质量在单一区间内二分反演；不以K矩阵指数替代一般曲面映射。
validateattributes(eta,{'double'},{'real','finite','vector','>=',0,'<=',1});eta=eta(:)';
coordinates=zeros(2,numel(eta));coordinates(:,eta==1)=repmat(map.coordinate,1,nnz(eta==1));
coordinates(:,eta==0)=repmat(map.apexCoordinate,1,nnz(eta==0));
active=eta>0&eta<1;if ~any(active),return;end
targets=log(eta(active));L=map.logMassRatio;
assert(all(targets>=L(end)),'ShockLayer:MassRange','入口质量比超出已分辨回溯范围，不外推。');
indices=zeros(size(targets));
for j=1:numel(targets)
    indices(j)=find(L>=targets(j),1,'last');
    if indices(j)==numel(L),indices(j)=indices(j)-1;end
end
t=map.traceParameter;h=t(indices+1)-t(indices);lo=zeros(size(targets));hi=ones(size(targets));
for iteration=1:48
    mid=(lo+hi)/2;value=hermite(L(indices),L(indices+1), ...
        h.*map.logMassDerivative(indices),h.*map.logMassDerivative(indices+1),mid);
    earlier=value>targets;lo(earlier)=mid(earlier);hi(~earlier)=mid(~earlier);
end
fraction=(lo+hi)/2;coordinates(:,active)=hermite(map.sourceCoordinateGridM(:,indices), ...
    map.sourceCoordinateGridM(:,indices+1),h.*map.coordinateDerivative(:,indices), ...
    h.*map.coordinateDerivative(:,indices+1),fraction);
end
function value=hermite(a,b,da,db,s)
value=(2*s.^3-3*s.^2+1).*a+(-2*s.^3+3*s.^2).*b+(s.^3-2*s.^2+s).*da+(s.^3-s.^2).*db;
end
