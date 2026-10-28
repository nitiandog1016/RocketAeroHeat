function value=applyRocketTemperatureProjector(projector,rhs)
%APPLYROCKETTEMPERATUREPROJECTOR 作用固定对称投影Pi，正向与伴随使用同一公式。
assert(strcmp(projector.kind,'fixed-current-temperature-row-nullspace-in-QR-right-coordinates')&& ...
    isequal(size(rhs),[projector.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'), ...
    'TemperatureProjection:Input','固定温度投影或五量输入非法。');
v=rhs(:);
if isfield(projector,'basisCoordinateIndices')
    % 紧凑基只作用于P转置H转置实际涉及的坐标；支撑以外严格保留输入。
    ids=projector.basisCoordinateIndices;
    assert(size(projector.basis,1)==numel(ids)&&numel(unique(ids))==numel(ids)&& ...
        all(ids>=1&ids<=numel(v)),'TemperatureProjection:Input','紧凑投影坐标不符。');
    selected=v(ids);v(ids)=selected-projector.basis*(projector.basis'*selected);
else
    v=v-projector.basis*(projector.basis'*v);
end
value=reshape(v,size(rhs));
assert(all(isfinite(value),'all'),'TemperatureProjection:Nonfinite','固定温度投影非有限。');
end
