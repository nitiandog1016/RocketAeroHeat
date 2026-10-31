function result=attachSU2ThermalBoundary(result,boundary)
%ATTACHSU2THERMALBOUNDARY 把实际使用的热边界带回面元、探针和导出，保持输入/求解身份。
assert(numel(boundary.givenTemperatureByFace)==numel(result.pressureRaw),'边界与原面元数量不匹配。');
if ~boundary.stockExecutableSupported
    assert(isfield(result.receipt,'thermalProfileApplied') && result.receipt.thermalProfileApplied.passed && ...
        strcmp(result.receipt.runtimeFlavor,'rah-axial-wall-v1'),'空间壁温缺少真实逐节点施加回执。');
    applied=result.receipt.thermalProfileApplied;
    assert(isequal(applied.profileX_m(:),boundary.specification.x(:)) && ...
        isequal(applied.profileTemperature_K(:),boundary.specification.temperature(:)), ...
        '空间壁温显示输入与实际施加温度表不一致。');
    result.solverMethod=[result.solverMethod '；RAH原生轴向壁温输入适配v1'];
end
result.thermalBoundary=boundary;result.wallTemperatureInput=boundary.givenTemperatureByFace;
result.heatFluxIsBoundaryInput=boundary.heatFluxIsBoundaryInput;
if boundary.heatFluxIsBoundaryInput
    assert(all(result.heatFluxRaw==0),'绝热原生热流输出与给定零值不一致。');
    result.resultStatus=[result.resultStatus '; adiabatic-q-is-prescribed-not-predicted'];
end
end
