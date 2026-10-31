function testKinneyProvenanceExport(out,baselineFile,selectedFile)
%TESTKINNEYPROVENANCEEXPORT 修复后核对原方法来源，并从真实保存结果导出。
% 模块一：只重放侧壁公式，不重算外流或网格；逐位核对已保存最终物理场。
% 模块二：逐面核对原算法名称、温度身份和范围，包含鼻区与脱体分支。
% 模块三：拒绝缺失/损坏来源，导出新目录；仅保存小型新增来源，不重复MAT。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(~exist(out,'dir'));mkdir(out);
b=load(baselineFile,'current','cfg','fs');s=load(selectedFile,'current');
cfg=b.cfg;cfg.heating.sidewallModel='kinneyAcreage';base=b.current;
base.heating.sidewallModel='kinneyAcreage';
current=calcKinneyAcreageHeating(base,b.fs,cfg);
for key={'heatFlux','recoveryTemperature','faceMethodID','eligibleFace', ...
        'heatTransferCoefficient','wallTemperatureByFace'}
    assert(isequaln(current.heating.(key{1}),s.current.heating.(key{1})),key{1});
end
assert(isequaln(current.pressure,s.current.pressure));
assert(isequaln(current.boundaryLayerEdge,s.current.boundaryLayerEdge));
data=buildAerothermalViewerData(current);baseline=buildAerothermalViewerData(b.current);
q=data.fields(strcmp({data.fields.key},'baselineHeatFlux'));
original=baseline.fields(strcmp({baseline.fields.key},'heatFlux'));
ids=double(current.heating.kinneyAcreage.baselineFaceMethodID(:));
for i=1:numel(ids)
    assert(q.status(i)==original.status(i) && isequaln(q.values(i),original.values(i)));
    assert(strcmp(q.methodNames{q.methodIndex(i)+1},original.methodNames{original.methodIndex(i)+1}));
    label=q.methodNames{q.methodIndex(i)+1};scope=q.scopeNames{q.scopeIndex(i)+1};
    if ids(i)==3
        assert(contains(label,'余弦延拓') && ~contains(label,'Thwaites'));
        assert(contains(scope,b.current.heating.noseScope));
    elseif ismember(ids(i),[4 5])
        assert(contains(label,'脱体闭合外流') && contains(scope,'受限前体脱体'));
    end
end
assert(nnz(ids==3)==1811 && nnz(ids==4)==335);
bad=current;bad.heating.kinneyAcreage=rmfield(bad.heating.kinneyAcreage,'baselineProvenance');
expect(@()buildAerothermalViewerData(bad));
bad=current;bad.heating.kinneyAcreage.baselineFaceMethodID(find(ids==3,1))=0;
expect(@()buildAerothermalViewerData(bad));
exportAerothermalWeb(current,fullfile(out,'workbench'));
baselineProvenance=current.heating.kinneyAcreage.baselineProvenance;
save(fullfile(out,'provenance.mat'),'baselineProvenance','-v7');
summary=struct('status','passed','faces',numel(ids),'noseSourceFaces',nnz(ids==3), ...
    'detachedLaminarSourceFaces',nnz(ids==4),'rejectionChecks',2, ...
    'physicalFieldsBitwiseUnchanged',true,'newPhysicalArea',0,'fullChainRecomputed',false);
fid=fopen(fullfile(out,'summary.json'),'w','n','UTF-8');fprintf(fid,'%s\n',jsonencode(summary));fclose(fid);
fprintf('Original per-face provenance verified for %d faces; all physical fields unchanged.\n',numel(ids));
end
function expect(action)
try,action();error('Test:UnexpectedSuccess','预期拒绝没有发生');
catch ex,assert(strcmp(ex.identifier,'RocketAeroHeat:Viewer:BaselineProvenance'),ex.identifier);end
end
