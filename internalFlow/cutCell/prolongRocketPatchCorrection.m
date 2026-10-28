function fine=prolongRocketPatchCorrection(partition,coarse)
%PROLONGROCKETPATCHCORRECTION 将每组的五量代数校正延拓至所有成员CV。
% 模块一：只在预条件层次中共享组常数，不合并物理单元或平均最终物理场。
% 模块二：保留原物理CV顺序，与固定等权限制配对。
assert(isequal(size(coarse),[partition.patchCount 5])&&isreal(coarse)&&all(isfinite(coarse),'all'),'PatchTransfer:Input','粗层五量尺寸或数值非法。');
fine=coarse(double(partition.groupID),:);
end
