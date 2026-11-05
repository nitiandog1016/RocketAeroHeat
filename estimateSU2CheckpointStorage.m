function budget=estimateSU2CheckpointStorage(volume,iterations,every,named,limitMiB,flowModel,thermodynamicModel)
%ESTIMATESU2CHECKPOINTSTORAGE 运行前按节点数和保存频率估计ASCII检查点体量。
% 模块一：计算基础文件与具名文件对的保守数量；末步单独输出也计入。
% 模块二：使用每节点文本上界估算，不把估计称为实测磁盘或进程内存。
% 模块三：默认100MiB保留；完整域可显式配置上限，不能静默放宽旧快照。
if nargin<5,limitMiB=100;end
if nargin<6,flowModel='laminar';end
if nargin<7,thermodynamicModel='idealGas';end
assert(ismember(flowModel,{'laminar','SST-2003m'}));
validateattributes(limitMiB,{'numeric'},{'scalar','real','finite','>=',1,'<=',1024});
validateattributes(iterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(every,{'numeric'},{'scalar','integer','positive'});
validateattributes(named,{'logical'},{'scalar'});
every=min(iterations,every);wallNodes=numel(unique([volume.wall{:}]));
pairBytes=240*size(volume.vertices,1)+600*wallNodes+1024^2;
% 新SST非紧凑文件同时保留高精度primitive；旧紧凑估计会低估真实保存成本。
if strcmp(flowModel,'SST-2003m'),pairBytes=640*size(volume.vertices,1)+720*wallNodes+1024^2;end
indexedPairs=0;if named,indexedPairs=floor((iterations-1)/every)+2;end
budget=struct('estimatedBytes',(1+indexedPairs)*pairBytes,'estimatedPairBytes',pairBytes, ...
    'maximumIndexedPairs',indexedPairs,'limitBytes',limitMiB*1024^2, ...
    'semantics','ASCII checkpoint estimate; excludes mesh, input snapshot, process memory and filesystem allocation');
budget.withinLimit=budget.estimatedBytes<=budget.limitBytes;
if strcmp(thermodynamicModel,'idealAir')
    assert(strcmp(flowModel,'SST-2003m') && named,'变比热SST需要具名阶段保存。');
    % 与Python output_plan逐项相同，包含首帧、末次别名和历史/日志预留。
    volumeBytes=size(volume.vertices,1)*32*27+4096;
    surfaceBytes=wallNodes*32*27+numel(volume.wall)*128+8192;
    pairs=min(iterations,2+floor((iterations-1)/every));
    budget.estimatedBytes=(pairs+1)*(volumeBytes+surfaceBytes)+iterations*4096+1024^2;
    budget.estimatedPairBytes=volumeBytes+surfaceBytes;budget.maximumIndexedPairs=pairs;
    budget.withinLimit=budget.estimatedBytes<=budget.limitBytes;
    budget.semantics='stage-aware ASCII file length estimate; excludes initial state, mesh and process memory';
end
end
