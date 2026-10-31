function report=reviewActualRocketNearWallStudy(out)
%REVIEWACTUALROCKETNEARWALLSTUDY 复用科研组件核验同源两网格结果，不重复运行流场。
% 模块一：原19752面完整映射、给定Tw/未知求解Tw及具名检查点身份。
% 模块二：同网格累计与改网格局部身份分别进入云图、探针、截面和CSV。
assert(~isfolder(out));mkdir(out);old=load('outputs/v119_validation/review01/result.mat','r');
report=struct;
for name={'control','candidate'}
    label=name{1};folder=fullfile('outputs/v120_validation/study01',label);
    receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
    assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut);
    p=load(fullfile(folder,'projection.mat'),'projection');
    r=mapSU2WallResult(readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk')),p.projection);
    manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
    r.V=old.r.V;r.F=old.r.F;r.receipt=receipt;r.flight=manifest.flight;
    r.sourceRevision=['v120-fixed-exterior-nearwall-' label];r.solveId=receipt.request.jobToken;
    r=attachSU2ThermalBoundary(r,old.r.thermalBoundary);r.convergenceAudit=auditSU2ResultConvergence(folder,r);
    a=r.convergenceAudit;assert(a.checkpointIdentityConfirmed,a.reason);assert(a.checkpointIteration==199);
    if strcmp(label,'control'),assert(a.cumulativeCheckpointIteration==799);
    else,assert(~isfield(a,'cumulativeCheckpointIteration') && a.sourceCumulativeIteration==599);end
    assert(~r.productionQualified && numel(r.heatFluxRaw)==19752 && all(r.wallTemperatureInput==300));
    fig=uifigure('Name',['RocketAeroHeat · 近壁对照 ' label],'Position',[30 30 1460 850]);g=onCleanup(@()close(fig));
    host=uipanel(fig,'Units','normalized','Position',[0 0 1 1]);ui=createSU2ResearchView(host,r);
    if strcmp(label,'candidate'),assert(contains(ui.auditText.Text,'新网格局部'));end
    ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(-1.5);
    t=tic;while toc(t)<40
        drawnow;s=ui.geometry.html.UserData;
        if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
        pause(.1);
    end
    assert(s.rendered.webglError==0 && s.rendered.revision==s.revision && ~isempty(ui.getSection()));
    ui.inspectFace(12000,[]);ui.exportTo(fullfile(out,[label '_fields.csv']));ui.exportSectionTo(fullfile(out,[label '_section.csv']));
    for suffix={'_fields.csv.json','_section.csv.json'}
        m=jsondecode(fileread(fullfile(out,[label suffix{1}])));assert(isequaln(m.convergenceAudit,jsondecode(jsonencode(a))));
    end
    save(fullfile(out,[label '_result.mat']),'r','-v7');section=ui.getSection();
    report.(label)=struct('audit',a,'originalFaces',19752,'fullArea',sum(r.faceArea), ...
        'componentWebGLPassed',true,'sectionRows',size(section.segments,1),'physicsQualified',false,'newQualifiedArea',0);
    clear g;
end
lines=describeSU2Convergence(old.r);assert(contains(lines{2},'累计 599'));
report.legacySnapshotPreserved=true;report.mainEndToEndLaunchClaimed=false;report.solverRerun=false;
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
