function snapshot=captureAerothermalFlowSnapshot(mesh)
%CAPTUREAEROTHERMALFLOWSNAPSHOT 保存当前工况的四个外流显示字段及真实来源。
% 模块作用：在外流建立时取一次快照，壁温迭代复用它；不保存重复几何、
% 流线或热壁字段，也不从供材料反馈的 referenceLaw 占位值恢复外流。
small=struct('V',mesh.V,'F',mesh.F);
for name={'pressure','boundaryLayerEdge'}
    if isfield(mesh,name{1}),small.(name{1})=mesh.(name{1});end
end
data=buildAerothermalViewerData(small);
fields=data.fields;
for j=1:numel(fields)
    if startsWith(fields(j).key,'edge') && isfield(mesh.boundaryLayerEdge,'scope')
        fields(j).scope=[fields(j).scope '；' mesh.boundaryLayerEdge.scope];
    end
end
fs=mesh.pressure.freestream;
snapshot=struct('fields',fields,'faceCount',size(mesh.F,1), ...
    'condition',[fs.Mach fs.altitude fs.alpha fs.beta], ...
    'flowDirection',fs.unitFlowDirection, ...
    'timeNote',sprintf('准定常外流工况：M=%.4g，H=%.6g m，攻角 %.4g°，侧滑 %.4g°', ...
    fs.Mach,fs.altitude,fs.alpha,fs.beta));
end
