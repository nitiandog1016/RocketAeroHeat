function report=reviewSchaeferCollarStudy(out)
%REVIEWSCHAEFERCOLLARSTUDY 复用两份终态场，检查原面映射和生产结果组件；不再运行求解器。
% 模块一：原图11读数和缺测/重叠/测站状态，输出科学表而非自动构造实验边界。
% 模块二：两网格都映射同一3648原面，保留绝热给定零与未收敛身份。
% 模块三：复用工作台结果组件的真实WebGL和CSV出口，无新增永久平行界面。
assert(~isfolder(out));mkdir(out);data=readSchaeferHeatTransferData();
assert(height(data.points)==40 && nnz(data.points.missing)==4 && nnz(data.points.broken)==4);
assert(all(isnan(data.points.St(data.points.missing))) && ~data.physicsQualified);
assert(~data.heatFluxConversionPerformed && ~data.wallTemperatureExtended);
assert(any(~data.points.nominalStationMatched & ~data.points.missing));
writetable(data.points,fullfile(out,'paired_heat_data.csv'));save(fullfile(out,'paired_heat_data.mat'),'data','-v7');
[surface,~]=buildSchaeferFergusonBody();names={'control','candidate'};records=cell(2,1);
fig=uifigure('Name','RocketAeroHeat · 近壁网格对照复核','Position',[20 30 1460 850]);guard=onCleanup(@()close(fig));
host=uipanel(fig,'Units','normalized','Position',[0 0 1 1]);
for j=1:2
    folder=fullfile('outputs','v116_validation','study01',names{j});saved=load(fullfile(folder,'projection.mat'));
    receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
    assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut);
    r=mapSU2WallResult(readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk')),saved.projection);
    r.V=surface.V;r.F=surface.F;r.receipt=receipt;r.flight=saved.manifest.flight;
    % 独立研究读取用本次唯一任务标识，不冒用主工作流的输入修订号。
    r.sourceRevision='v116-study01-research-replay';r.solveId=receipt.request.jobToken;
    boundary=resolveSU2ThermalBoundary(surface,struct('mode','adiabatic'));r=attachSU2ThermalBoundary(r,boundary);
    r.convergenceAudit=auditSU2ResultConvergence(folder,r);
    assert(strcmp(r.convergenceAudit.status,'history-audited-checkpoint-unconfirmed'));
    assert(r.convergenceAudit.terminalIteration==599 && ~r.convergenceAudit.terminalResidualPassed);
    assert(isnan(r.convergenceAudit.checkpointIteration) && ~r.productionQualified);
    % 本轮关闭了具名输出，绝热HF全零不能唯一绑定迭代；不能伪造599步检查点身份。
    delete(host.Children);ui=createSU2ResearchView(host,r);
    ui.select.Value='pressureRaw';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(.158496);waitRender(ui);
    assert(isequal(ui.geometry.html.UserData.payload.faceField.values,r.pressureRaw));
    ui.inspectFace(1454,[]);assert(contains(strjoin(ui.probe.Value,newline),'未指定（绝热壁）'));
    ui.exportTo(fullfile(out,[names{j} '_fields.csv']));ui.exportSectionTo(fullfile(out,[names{j} '_section.csv']));
    save(fullfile(out,[names{j} '_result.mat']),'r','-v7');
    records{j}=struct('case',names{j},'originalFaces',numel(r.pressureRaw), ...
        'rendered',ui.geometry.html.UserData.rendered,'convergence',r.convergenceAudit, ...
        'newQualifiedArea',0,'physicsQualified',false);
end
report=struct('cases',{records},'pairedRecords',height(data.points),'missingRecords',nnz(data.points.missing), ...
    'brokenRecords',nnz(data.points.broken),'stationMatchedRecords',nnz(data.points.nominalStationMatched), ...
    'productionResultComponentReused',true,'mainWorkflowAcceptanceClaimed',false, ...
    'solverOrMesherRerun',false,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Paired data and original-face production component exports checked; numerical cases remain unqualified.\n');
end
function waitRender(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision&&s.rendered.webglError==0&&~isempty(ui.getSection()));
end
