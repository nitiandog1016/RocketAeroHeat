function result=attachSU2ThermalBoundary(result,boundary)
%ATTACHSU2THERMALBOUNDARY 把实际使用的热边界带回面元、探针和导出，保持输入/求解身份。
assert(numel(boundary.givenTemperatureByFace)==numel(result.pressureRaw),'边界与原面元数量不匹配。');
assert(boundary.stockExecutableSupported,'尚未执行的空间壁温不得附给已有结果。');
result.thermalBoundary=boundary;result.wallTemperatureInput=boundary.givenTemperatureByFace;
result.heatFluxIsBoundaryInput=boundary.heatFluxIsBoundaryInput;
if boundary.heatFluxIsBoundaryInput
    assert(all(result.heatFluxRaw==0),'绝热原生热流输出与给定零值不一致。');
    result.resultStatus=[result.resultStatus '; adiabatic-q-is-prescribed-not-predicted'];
end
end
