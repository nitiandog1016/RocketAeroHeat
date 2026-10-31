function budget=estimateSU2CheckpointStorage(volume,iterations,every,named)
%ESTIMATESU2CHECKPOINTSTORAGE 运行前按节点数和保存频率估计ASCII检查点体量。
% 模块一：计算基础文件与具名文件对的保守数量；末步单独输出也计入。
% 模块二：使用每节点文本上界估算，不把估计称为实测磁盘或进程内存。
validateattributes(iterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(every,{'numeric'},{'scalar','integer','positive'});
validateattributes(named,{'logical'},{'scalar'});
every=min(iterations,every);wallNodes=numel(unique([volume.wall{:}]));
pairBytes=240*size(volume.vertices,1)+600*wallNodes+1024^2;
indexedPairs=0;if named,indexedPairs=floor((iterations-1)/every)+2;end
budget=struct('estimatedBytes',(1+indexedPairs)*pairBytes,'estimatedPairBytes',pairBytes, ...
    'maximumIndexedPairs',indexedPairs,'limitBytes',100*1024^2, ...
    'semantics','ASCII checkpoint estimate; excludes mesh, input snapshot, process memory and filesystem allocation');
budget.withinLimit=budget.estimatedBytes<=budget.limitBytes;
end
