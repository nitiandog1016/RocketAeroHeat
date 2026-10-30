function geometry=axisymmetricShockGraph(shape,z)
%AXISYMMETRICSHOCKGRAPH 共同轴对称二次曲线的正则半径平方坐标与曲率。
% 模块一：z=r^2只是局部坐标，绝不是世界坐标Z；K=k*I，B=-k为球面。
% 模块二：解析给出两主曲率、1/G、法向速度导数，轴线直接取有限表达式。
% 模块三：整个计算域必须为单值凸图；不修改形状或裁剪开方自变量。
validateattributes(z,{'double'},{'real','finite','nonnegative'});
validateattributes(shape.K,{'double'},{'real','finite','size',[2 2]});
k=trace(shape.K)/2;
assert(k>0&&norm(shape.K-k*eye(2),'fro')<1e-12*k, ...
    'ShockPatch:Axisymmetry','本模块需要共同激波K=k*I；不能把非轴对称输入取平均后计算。');
B=0;if isfield(shape,'profileB'),B=shape.profileB;end
validateattributes(B,{'double'},{'real','finite','scalar'});
L=1+B*k*z;W=1+k*(B+k)*z;
assert(all(L>0&W>0,'all'),'ShockPatch:Graph','共同激波超出单值凸图范围。');
root=sqrt(L);invG=root./sqrt(W);
geometry=struct('radiusSquaredM2',z,'radiusM',sqrt(z),'axialPositionM',k*z./(1+root), ...
    'inverseMeridionalMetric',invG,'meridionalCurvaturePerM',k./W.^1.5, ...
    'azimuthalCurvaturePerM',k./sqrt(W), ...
    'inverseMetricDerivativeByRadiusSquared',-k^2./(2*root.*W.^1.5), ...
    'apexCurvaturePerM',k,'profileBPerM',B);
end
