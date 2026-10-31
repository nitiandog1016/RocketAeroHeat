function budget=estimateSU2CheckpointStorage(volume,iterations,every,named,limitMiB,flowModel)
%ESTIMATESU2CHECKPOINTSTORAGE 运行前按节点数和保存频率估计ASCII检查点体量。
% 模块一：计算基础文件与具名文件对的保守数量；末步单独输出也计入。
% 模块二：使用每节点文本上界估算，不把估计称为实测磁盘或进程内存。
% 模块三：默认100MiB保留；完整域可显式配置上限，不能静默放宽旧快照。
if nargin<5,limitMiB=100;end
if nargin<6,flowModel='laminar';end
assert(ismember(flowModel,{'laminar','SST-2003m'}));
validateattributes(limitMiB,{'numeric'},{'scalar','real','finite','>=',1,'<=',1024});
validateattributes(iterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(every,{'numeric'},{'scalar','integer','positive'});
validateattributes(named,{'logical'},{'scalar'});
every=min(iterations,every);wallNodes=numel(unique([volume.wall{:}]));
pairBytes=240*size(volume.vertices,1)+600*wallNodes+1024^2;
% SST新增两列状态及湍流壁面字段，不能沿用层流体量低估保存预算。
if strcmp(flowModel,'SST-2003m'),pairBytes=296*size(volume.vertices,1)+640*wallNodes+1024^2;end
indexedPairs=0;if named,indexedPairs=floor((iterations-1)/every)+2;end
budget=struct('estimatedBytes',(1+indexedPairs)*pairBytes,'estimatedPairBytes',pairBytes, ...
    'maximumIndexedPairs',indexedPairs,'limitBytes',limitMiB*1024^2, ...
    'semantics','ASCII checkpoint estimate; excludes mesh, input snapshot, process memory and filesystem allocation');
budget.withinLimit=budget.estimatedBytes<=budget.limitBytes;
end
