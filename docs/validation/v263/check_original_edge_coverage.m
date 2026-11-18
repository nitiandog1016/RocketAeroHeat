function check_original_edge_coverage(folder)
%CHECK_ORIGINAL_EDGE_COVERAGE 从原 MAT 与工程快照核对区域页实际接收的超域和未知面元。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root,fullfile(root,'internalFlow','cutCell'));
native=fullfile(root,'outputs','v104_validation','baseline02','current_case.mat');
archive=fullfile(root,'outputs','v176_validation','main02','engineering_data.json');
assert(strcmp(rocketStreamingFileSHA256(native),'b4889e42a28131c253a56f06472a58c9d2cb290d7c9937935fbb4b5d7349488f'));
assert(strcmp(rocketStreamingFileSHA256(archive),'9518d1d1f44abb83f60664a8c1d072d25797ce16d2182823900cfdf119110b1b'));
s=load(native,'current');saved=jsondecode(fileread(archive));
view=buildAerothermalViewerData(s.current);ledger=buildAerothermalEvidenceLedger(view);
wire=jsondecode(jsonencode(view,'ConvertInfAndNaN',true));
% 模块一：在同一工况、同一原面编号下，八个场量分别与冻结版本的原值和身份对照。
assert(isequal(wire.vertices,saved.vertices) && isequal(wire.faces,saved.faces) && isequal(wire.faceArea,saved.faceArea));
assert(isequaln(wire.resultCondition,saved.resultCondition));
assert(isequal({wire.fields.key},{saved.fields.key}));
for k=1:numel(view.fields)
    assert(isequaln(wire.fields(k).values,saved.fields(k).values) && isequal(wire.fields(k).status,saved.fields(k).status));
end
% 模块二：独立取原生拒绝标志与上次追溯清单。直接失效面、下游原面分别比对。
src=fullfile(root,'docs','validation','v263','evidence','sidewall_temperature');
source=jsondecode(fileread(fullfile(src,'direct_gap_status.json')));direct=[source.faceRecords.faceID]';
tbl=readtable(fullfile(src,'faces.csv'));down=tbl.faceID;
outside=logical(s.current.boundaryLayerEdge.thermalExpansionOutOfRangeFace(:));
assert(numel(direct)==6 && numel(unique(down))==936 && all(outside(direct)) && ~any(outside(down)));
for key={'edgePressure','edgeTemperature'}
    f=ledger.fields(strcmp({ledger.fields.key},key{1}));
    old=saved.evidenceLedger.fields(strcmp({saved.evidenceLedger.fields.key},key{1}));
    assert(isequal(f.availability==4,outside) && all(old.availability(outside)==0));
    assert(all(f.qualification(outside)==1) && all(isnan(f.values(outside))));
    assert(all(f.availability(down)==1) && all(f.qualification(down)==0));
end
for key={'heatFlux','recoveryTemperature'}
    f=ledger.fields(strcmp({ledger.fields.key},key{1}));
    check=[direct;down];assert(all(f.availability(check)==0) && all(f.qualification(check)==0) && all(isnan(f.values(check))));
end
q=ledger.fields(strcmp({ledger.fields.key},'pressure'));
assert(all(q.availability(direct)==3) && all(q.availability(down)==1));
assert(isequal(ledger.joint.pressureHeatFluxEngineeringFaceIDs,saved.evidenceLedger.joint.pressureHeatFluxEngineeringFaceIDs));
assert(ledger.joint.pressureHeatFluxQualifiedArea==0);
% 模块三：导出七个查询示例的实际字段。NaN 在 JSON 中保留为 null。
chosen=[direct;down(1)];probes=struct([]);
for j=1:numel(chosen)
    id=chosen(j);probes(j).faceID=id;probes(j).areaM2=view.faceArea(id);
    for k=1:numel(ledger.fields)
        f=ledger.fields(k);
        probes(j).fields.(f.key)=struct('value',f.values(id),'unit',f.unit, ...
            'method',f.methodNames{f.methodIndex(id)+1},'availability',f.availability(id), ...
            'qualification',f.qualification(id),'missingReason',f.missingReason{id});
    end
end
report=struct('originalFaces',numel(outside),'sourceDomainRefusals',nnz(outside),'sourceDomainAreaM2',sum(view.faceArea(outside)), ...
    'checkedDirectFaceIDs',direct,'checkedDownstreamFaces',numel(down),'matchedValueAndStatusFields',{ {view.fields.key} }, ...
    'jointEngineeringAreaM2',ledger.joint.pressureHeatFluxEngineeringArea,'jointQualifiedAreaM2',ledger.joint.pressureHeatFluxQualifiedArea, ...
    'newJointEngineeringAreaM2',0,'flowSolverRun',false,'browserRun',false,'probeExamples',probes);
fid=fopen(fullfile(folder,'report.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
end
