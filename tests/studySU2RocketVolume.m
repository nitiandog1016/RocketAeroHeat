function studySU2RocketVolume(out,pythonExecutable)
%STUDYSU2ROCKETVOLUME 真实三维完整箭体几何与SU2预处理检查，不求物理覆盖。
% 模块一：保留全部六类区域及原三角形，检验体域正号/闭合/面积/体积。
% 模块二：0°与20°/5°侧滑进行实际SU2 dryrun，检验来流三分量约定。
if nargin<2,pythonExecutable='python';end
assert(~exist(out,'dir'));mkdir(out);
[surface,truth]=buildRocketTestBody(struct('resolution',24,'protrusionEnabled',false));
volume=buildSU2RocketVolume(surface,truth,struct('radialLayers',8,'maxWallAngleStepDeg',10,'baseRadialSegments',4));
assert(volume.topology.declaredBoundaryExactlyMatched && volume.topology.nonManifoldFaces==0);
assert(sum(strcmp(volume.wallRegion,'base'))>0 && sum(strcmp(volume.wallRegion,'contraction'))>0);
expected={'afterbody';'base';'cone';'contraction';'cylinder';'sphericalNose'};
assert(isequal(sort(unique(volume.wallRegion)),expected));
angles=[0 0;20 5];receipts=cell(2,1);cases=cell(2,1);
for k=1:2
    folder=fullfile(out,sprintf('angle%02d',k));
    cases{k}=writeSU2RocketCase(volume,folder,struct('alpha',angles(k,1),'beta',angles(k,2)));
    receipts{k}=runEmbeddedSU2Case(folder,struct('pythonExecutable',pythonExecutable,'threads',2,'maximumSeconds',120,'dryRun',true));
    log=fileread(fullfile(folder,'engine_stdout.log'));
    assert(contains(log,'Exit Success'),'SU2预处理没有成功终态。');
end
save(fullfile(out,'volume_and_source.mat'),'volume','surface','truth','cases','receipts','-v7');
fprintf('Complete 3D rocket volume: %d cells, %d original faces, all six regions; actual SU2 dryruns completed. No physical qualification.\n',numel(volume.cells),size(surface.F,1));
end
