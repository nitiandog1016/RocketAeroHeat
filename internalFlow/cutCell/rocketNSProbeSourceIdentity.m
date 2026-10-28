function identity=rocketNSProbeSourceIdentity
%ROCKETNSPROBESOURCEIDENTITY 绑定原NS及显式精确储库依赖，防止跨算子复用方向。
identity=rocketSelectedDirectionSourceIdentity;
root=fileparts(mfilename('fullpath'));
for name={'resolveNSRiemannReservoir3D','buildCutCellRiemannReservoirBoundaryState'}
    identity.(name{1})=rocketStreamingFileSHA256(fullfile(root,[name{1} '.m']));
end
end
