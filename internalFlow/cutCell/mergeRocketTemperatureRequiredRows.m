function rows=mergeRocketTemperatureRequiredRows(old,new)
%MERGEROCKETTEMPERATUREREQUIREDROWS 合并同一当前状态的全部必要温度行，不删旧约束。
% 模块一：共同面必须逐项保持H、符号和原CV/波系元数据；冲突明确拒绝。
% 模块二：新增面按原面号插入，不能按秩选择或近似合并，仍使用原数组预算。
assert(old.cellCount==new.cellCount&&size(old.H,2)==size(new.H,2), ...
    'TemperatureRepreview:Identity','必须是同一当前CV排列的温度行。');
assert(all(isfield(old,{'scale','gas'}))&&all(isfield(new,{'scale','gas'}))&& ...
    isequaln(old.scale,new.scale)&&isequaln(old.gas,new.gas), ...
    'TemperatureRepreview:Identity','温度行气体与五量缩放必须相同，包括空约束集合。');
if isempty(new.faceIDs),rows=old;return;end
if isempty(old.faceIDs)
    rows=new;rows.options=old.options;
    rows.options.temperatureBandK=max(old.options.temperatureBandK,new.options.temperatureBandK);return;
end
keys={'faceIDs','ownerCV','neighbourCV','normal','QL','QR','H','leftRowsK','rightRowsK', ...
    'baseSamplePrimitive','baseTemperatureK','remainingUpperTemperatureK','sampleRegion', ...
    'rowTwoNormK','temperatureBoundSigns'};
[common,a,b]=intersect(old.faceIDs,new.faceIDs); %#ok<ASGLU>
for key=keys
    assert(isequal(old.(key{1})(a,:),new.(key{1})(b,:)), ...
        'TemperatureRepreview:Retention','同一原面出现不同H、侧别或元数据，拒绝覆盖。');
end
fresh=~ismember(new.faceIDs,old.faceIDs);rows=old;
if any(fresh)
    [~,order]=sort([old.faceIDs;new.faceIDs(fresh)]);
    bound=24*(5*old.cellCount+1)+4096*numel(order)+65536;
    assert(bound<=64*1024^2,'TemperatureRepreview:Budget','合并全部行超过原64MiB构造预算，不截断。');
    for key=keys
        joined=[old.(key{1});new.(key{1})(fresh,:)];rows.(key{1})=joined(order,:);
    end
end
rows.options.temperatureBandK=max(old.options.temperatureBandK,new.options.temperatureBandK);
end
