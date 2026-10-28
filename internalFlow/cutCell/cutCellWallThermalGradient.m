function [gradient,audit]=cutCellWallThermalGradient(cellCenter,wallCenter,normal,Tcell,Twall,boundary,wallID)
%CUTCELLWALLTHERMALGRADIENT 已知壁面两切向温度导数与真实位移共同确定三维墙梯度。
% 模块一：常数Tw的切向梯度为0；函数Tw须给完整切向向量，不推断缺失数据。
% 模块二：g=g_t+n*(Tw-Tc-g_t·delta)/(n·delta)，要求真实正壁距，不裁梯度/温度。
% 模块三：给出对Towner的精确线性导数，供同一扩散模板/步长核对使用。
assert(strcmp(cutCellWallThermalPolicy(boundary),'dirichlet-tangent'));
N=numel(Tcell);Tcell=Tcell(:);Twall=Twall(:);wallID=wallID(:);
assert(isequal(size(cellCenter),[N 3])&&isequal(size(wallCenter),[N 3])&&isequal(size(normal),[N 3])&& ...
    numel(Twall)==N&&numel(wallID)==N&&all(isfinite([cellCenter,wallCenter,normal,Tcell,Twall]),'all')&& ...
    all(abs(sum(normal.^2,2)-1)<1e-10),'CutNS:WallThermalGeometry','墙温/几何尺寸、有限值或法向不符。');
delta=wallCenter-cellCenter;dn=sum(delta.*normal,2);
assert(all(Tcell>0&Twall>0),'CutNS:WallThermalState','绝对温度必须为正。');
assert(all(isfinite(dn)&dn>0),'CutNS:WallThermalGeometry','真实法向墙距必须严格为正，不能加epsilon。');
tangent=zeros(N,3);hasGradient=isfield(boundary,'temperatureSurfaceGradient');
if hasGradient
    value=boundary.temperatureSurfaceGradient;
    if isa(value,'function_handle')
        for k=1:N
            input=value(wallCenter(k,:),wallID(k));
            validateattributes(input,{'numeric'},{'real','finite','numel',3});tangent(k,:)=double(input(:)');
        end
    else
        validateattributes(value,{'numeric'},{'real','finite','numel',3});tangent=repmat(double(value(:)'),N,1);
    end
end
assert(all(abs(sum(tangent.*normal,2))<=1e-10*max(1,vecnorm(tangent,2,2))), ...
    'CutNS:WallThermalBoundary','所给表面温度梯度必须切向；不自动投影错误输入。');
assert(isa(boundary.temperature,'function_handle')||all(tangent==0,'all'), ...
    'CutNS:WallThermalBoundary','常数墙温不能给非零切向温度导数。');
gn=(Twall-Tcell-sum(tangent.*delta,2))./dn;gradient=tangent+gn.*normal;
assert(all(isfinite(gradient),'all'),'CutNS:WallThermalState','墙梯度非有限，不裁剪。');
audit=struct('kind','explicit-three-dimensional-Dirichlet-tangent-wall-temperature-gradient', ...
    'wallID',wallID,'normal',normal,'normalDistanceM',dn,'delta',delta,'cellTemperatureK',Tcell, ...
    'wallTemperatureInputK',Twall,'surfaceGradientInputKm',tangent,'normalGradientKm',gn, ...
    'gradientXYZ',gradient,'ownerTemperatureDerivative',-normal./dn,'temperatureInputKind','prescribed', ...
    'stateClipped',false,'physicalValidation',false);
end
