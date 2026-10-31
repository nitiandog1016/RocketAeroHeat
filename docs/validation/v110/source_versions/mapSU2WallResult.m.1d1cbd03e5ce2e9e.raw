function result=mapSU2WallResult(data,projection)
%MAPSU2WALLRESULT 将实际VTK节点场映射到原三角形，保留来源与未资格状态。
% 模块一：核对输出节点与已知精确壁网格一一对应，允许VTK文本舍入但不补点。
% 模块二：按面积片常值投影p/q/流体节点T，保留有符号热功率。
% 模块三：无论数值是否有限，本接口不推定收敛/物理资格或给定壁温。
ids=projection.wallVertexIDs;exact=projection.sourceVertices(ids,:);n=size(exact,1);
assert(size(data.vertices,1)==n,'RocketAeroHeat:SU2:MappingNodes','输出壁节点数不匹配。');
% SU2 ASCII坐标为六位小数科学计数，绝对舍入界随坐标量级变化。
tolerance=8e-7*max(1,max(abs(exact(:))));match=zeros(n,1);distance=zeros(n,1);
for k=1:n
    [squared,j]=min(sum((exact-data.vertices(k,:)).^2,2));match(k)=ids(j);distance(k)=sqrt(squared);
end
assert(max(distance)<tolerance && numel(unique(match))==n,'RocketAeroHeat:SU2:MappingNodes', ...
    '壁节点无法唯一匹配，禁止最近面或重复节点补漏。');
actual=cell(numel(data.offsets),1);first=1;
for k=1:numel(actual)
    last=data.offsets(k);actual{k}=sort(match(data.connectivity(first:last)+1))';first=last+1;
end
expected=cellfun(@sort,projection.sourceCells,'UniformOutput',false);
assert(isequal(sortrows(pad(actual)),sortrows(pad(expected))), ...
    'RocketAeroHeat:SU2:MappingConnectivity','输出壁面连接与原体网格不一致。');
nv=size(projection.sourceVertices,1);fields=zeros(nv,3);
fields(match,:)=[data.pressureRaw,data.heatFluxRaw,data.temperatureRaw];mapped=projection.weights*fields;
sourcePower=projection.sourceNodeArea'*fields(:,2);targetPower=projection.targetFaceArea'*mapped(:,2);
scale=max(sum(projection.sourceNodeArea.*abs(fields(:,2))),realmin);
assert(abs(sourcePower-targetPower)/scale<2e-8,'RocketAeroHeat:SU2:MappingPower','有符号热功率不守恒。');
% 引擎使用节点合成面积向量模；保守传递的是其热负荷，不改写原始q或物面平均q。
factor=zeros(nv,1);factor(ids)=projection.engineNodeArea(ids)./projection.sourceNodeArea(ids);
loadDensity=projection.weights*(fields(:,2).*factor);
enginePower=projection.engineNodeArea'*fields(:,2);mappedLoadPower=projection.targetFaceArea'*loadDensity;
engineScale=max(sum(projection.engineNodeArea.*abs(fields(:,2))),realmin);
assert(abs(mappedLoadPower-enginePower)/engineScale<2e-8);
result=struct('pressureRaw',mapped(:,1),'heatFluxRaw',mapped(:,2),'fluidNodeTemperatureMeanRaw',mapped(:,3), ...
    'faceArea',projection.targetFaceArea,'sourceHeatPowerRaw',sourcePower,'targetHeatPowerRaw',targetPower, ...
    'heatPowerRelativeClosure',abs(sourcePower-targetPower)/scale,'maxNodeMatchingDistance',max(distance), ...
    'conservativeHeatLoadDensityRaw',loadDensity,'reconstructedEngineHeatPowerRaw',enginePower, ...
    'conservativeHeatLoadPowerRaw',mappedLoadPower,'heatLoadRelativeClosure',abs(mappedLoadPower-enginePower)/engineScale, ...
    'heatLoadMethod','native SU2 nodal heat load conservatively distributed by intersecting scalar median-dual patches', ...
    'heatFluxAndHeatLoadAreDistinct',true,'engineAreaRule',projection.engineAreaRule, ...
    'projectionMethod',projection.representation,'solverMethod',data.method, ...
    'temperatureIdentity','area-mean fluid-node state; neither prescribed Tw nor wall thermal solution', ...
    'quantityUnits','raw solver units; SI only when case contract verified', ...
    'resultStatus','mapped-research-output; convergence-and-physics-not-inferred', ...
    'eligibleForSurfaceCoverage',false,'productionQualified',false,'newQualifiedRocketArea',0);
end

function matrix=pad(faces)
matrix=zeros(numel(faces),4);
for k=1:numel(faces),matrix(k,1:numel(faces{k}))=faces{k};end
matrix=sort(matrix,2);
end
