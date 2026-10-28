function mapped=mapCutCellWallFields(cut,spatialInfo,steadyConverged)
%MAPCUTCELLWALLFIELDS 原曲面墙片到原三角面的保守面积积分映射（实验草稿）。
% 模块用途：保持原faceID，热功率按Σq*A、压力力按Σp*n_fluid*A计算。
% 每个变量分别统计有效面积；原面有任意未知部分时平均值仍为NaN。
% 给定Tw永远标输入；当前p/q即使有限也只是派生量，不能自动计入覆盖。
% steadyConverged仅描述数值迭代，不代表近壁/网格/物性/实验验证完成。
validateattributes(steadyConverged,{'logical'},{'scalar'});
V=cut.inputGeometry.V;F=cut.inputGeometry.F;n=size(F,1);w=cut.walls;ids=w.originalFaceID(:);A=w.area(:);
assert(all(ids>=1&ids<=n&ids==floor(ids))&&all(isfinite(A)&A>0),'CutMap:Geometry','原面编号或墙片面积无效。');
normal=w.normal;assert(all(abs(vecnorm(normal,2,2)-1)<1e-10),'CutMap:Normal','须使用可靠原曲面单位法向。');
av=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2)/2;
sourceArea=vecnorm(av,2,2);assert(all(isfinite(sourceArea)&sourceArea>0),'CutMap:Geometry','原三角面退化。');
sumArea=accumarray(ids,A,[n 1]);areaRelative=abs(sumArea-sourceArea)./sourceArea;
assert(all(areaRelative<1e-9),'CutMap:MissingGeometry','墙片几何有漏片/重复或原三角身份不一致。');
% 同时校核有向面积；法向相反的错误不能被标量面积和掩盖。
sumAV=zeros(n,3);for d=1:3,sumAV(:,d)=accumarray(ids,A.*normal(:,d),[n 1]);end
vectorRelative=vecnorm(sumAV+av,2,2)./sourceArea;
assert(all(vectorRelative<1e-9),'CutMap:Normal','流体墙法向与原实体外法向不一致。');
required={'wallActiveInFluidDomain','wallPressure','wallHeatFluxIntoWall','wallTemperatureInput'};
assert(all(isfield(spatialInfo,required)),'CutMap:Fields','缺失墙通量或边界身份。');
active=spatialInfo.wallActiveInFluidDomain(:);
assert(isequal(active,w.ownerCell(:)>0),'CutMap:Owner','流体有效墙身份不一致。');
pressure=fieldMap(spatialInfo.wallPressure,'Pa','derived-current-fluid-state','当前真实CV压力的原三角面面积平均');
heat=fieldMap(spatialInfo.wallHeatFluxIntoWall,'W/m^2','derived-current-fluid-state','Fourier真实曲面法向通量，按原墙片热功率守恒映射');
wallInput=fieldMap(spatialInfo.wallTemperatureInput,'K','prescribed','明确给定壁温在原三角面上的面积平均');
forcePartial=zeros(n,3);valid=active&isfinite(spatialInfo.wallPressure(:));p=spatialInfo.wallPressure(:);
for d=1:3,forcePartial(:,d)=accumarray(ids(valid),A(valid).*p(valid).*normal(valid,d),[n 1]);end
force=nan(n,3);force(pressure.complete,:)=forcePartial(pressure.complete,:);
power=nan(n,1);power(heat.complete)=heat.integralKnownPart(heat.complete);
knownPower=sum(heat.integralKnownPart);patchPower=sum(A(active&isfinite(spatialInfo.wallHeatFluxIntoWall(:))).* ...
    spatialInfo.wallHeatFluxIntoWall(active&isfinite(spatialInfo.wallHeatFluxIntoWall(:))));
mapped=struct('originalFaceID',(1:n)','sourceTriangleArea',sourceArea, ...
    'mappedWallPatchArea',sumArea,'geometryAreaMaximumRelativeError',max(areaRelative), ...
    'geometryVectorMaximumRelativeError',max(vectorRelative), ...
    'pressure',pressure,'heatFluxIntoWall',heat,'wallTemperatureInput',wallInput, ...
    'heatPowerIntoWall',power,'pressureForceOnBody',force, ...
    'pressureForceKnownPart',forcePartial,'heatPowerKnownPart',heat.integralKnownPart, ...
    'summedKnownHeatPower',knownPower,'originalPatchKnownHeatPower',patchPower, ...
    'heatPowerMappingDifference',knownPower-patchPower,'steadyConverged',steadyConverged, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'newSolvedSurfaceAreaM2',0, ...
    'scope','实验墙片守恒映射；数值稳态与模型/近壁/网格/实验适用性独立，不能自动覆盖生产工程场');

    function out=fieldMap(values,units,kind,method)
        values=values(:);assert(numel(values)==numel(A)&&isreal(values),'CutMap:Fields','墙场形状或实数合同无效。');
        assert(all(isnan(values(~active))),'CutMap:Inactive','域外无流体墙必须保留NaN，不能给定或补零。');
        valid=active&isfinite(values);knownArea=accumarray(ids(valid),A(valid),[n 1]);
        integral=accumarray(ids(valid),A(valid).*values(valid),[n 1]);
        % 面积容差只处理累加舍入；任意未知正片都使该原面不完整。
        unknownCount=accumarray(ids(~valid),ones(nnz(~valid),1),[n 1]);
        complete=unknownCount==0 & knownArea>0 & abs(knownArea-sourceArea)./sourceArea<1e-9;
        average=nan(n,1);average(complete)=integral(complete)./knownArea(complete);
        status=repmat({'unknown-original-wall-patch'},n,1);
        if strcmp(kind,'prescribed'),label='prescribed-boundary-input';
        elseif steadyConverged,label='derived-numerically-steady-model-unvalidated';
        else,label='derived-not-steady';end
        status(complete)={label};
        out=struct('value',average,'units',units,'valueKind',kind,'method',method, ...
            'status',{status},'complete',complete,'knownArea',knownArea, ...
            'unknownPatchCount',unknownCount,'integralKnownPart',integral, ...
            'eligibleForSurfaceCoverage',false);
    end
end
