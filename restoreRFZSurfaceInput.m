function [mesh,geometry,input,wall] = restoreRFZSurfaceInput(s,raw,geometry)
% 模块用途：由原始米制几何重建朝向、编号、分区和壁温。通过后才交给工作台。
assert(isstruct(s) && isscalar(s) && strcmp(s.schema,'rfz-rec-surface-input-v1') && strcmp(s.domain,'fullCAD') && ...
    strcmp(s.stlSHA256,raw.rfzBoundary.stlSHA256),'RAH:RFZ:SavedInput','所选文件与 RFZ 完整 CAD 不符。');
pose=str2double(string(s.poseDecimal(:)));
assert(numel(pose)==12 && all(isfinite(pose)),'RAH:RFZ:SavedInput','保存的朝向无效。');
rotation=reshape(pose(1:9),3,3);center=reshape(pose(10:12),1,3);
if isequal(rotation,eye(3)),mesh=raw;else,mesh=alignGeometryAxis(raw,'matrix',rotation,center);end
geometry.rotationMatrix=rotation;geometry.rotationCenter=center;
geometry.isAligned=true;geometry.alignmentMethod='workbenchConfirmed';
geometry.bodyAxis=[1 0 0];geometry.noseDirection=[1 0 0];
[input,wall]=buildRFZSurfaceInput(mesh,geometry,s.cavityMode,s.cavityWallK);
assert(isequaln(input,s),'RAH:RFZ:SavedInput','保存输入的网格、温界或来源与重建内容不一致。');
end
