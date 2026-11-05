function report=testSU2IdealAirLegacyReads(out)
%TESTSU2IDEALAIRLEGACYREADS 新EOS分支加入后，只读复核两份既有小型保存态。
% 模块用途。检查旧HEOS和定比热SST的具名时刻及EOS读取，无新求解或窗口。
assert(~isfolder(out));mkdir(out);
folders={fullfile(pwd,'outputs','v121_validation','capability01'), ...
    fullfile(pwd,'outputs','v124_validation','plate01','plate')};
expected=[0 99];rows=cell(1,2);
for j=1:2
    folder=folders{j};wall=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
    weight=zeros(size(wall.vertices,1),1);start=1;
    for i=1:numel(wall.offsets)
        ids=wall.connectivity(start:wall.offsets(i))+1;assert(numel(ids)==2);
        weight(ids)=weight(ids)+norm(wall.vertices(ids(1),:)-wall.vertices(ids(2),:))/2;
        start=wall.offsets(i)+1;
    end
    receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
    result=struct('receipt',receipt,'reconstructedEngineHeatPowerRaw',weight'*wall.heatFluxRaw);
    audit=auditSU2ResultConvergence(folder,result);
    assert(audit.checkpointIdentityConfirmed,audit.reason);
    assert(audit.checkpointIteration==expected(j) && ~audit.productionQualified);
    rows{j}=audit;
end
report=struct('passed',true,'cases',{rows},'nativeRuns',0,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
