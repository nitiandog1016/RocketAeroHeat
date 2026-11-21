function check_domain_ui(out)
% 模块一：按带线模型和中段加密网走主侧栏，先读取给定壁温，再切换到整件净功率。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,1);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
for k=1:2
    mode='givenWall';
    if k==2
        mode='givenPower';press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;choose(c.beadCondition,'power');
        assert(c.beadPower.Value==.004 && isempty(fig.UserData.workflowState.result));
    end
    press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
    s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
    r=s.result;p=r.thermistorOuterDomain;assert(p.available && strcmp(p.mode,mode));
    if k==1,assert(p.wall30K==305.15 && p.wall20K==305.15 && isfield(r,'thermistorLeadGrid'));
    else,assert(abs(p.outward30W-.004)<1e-13 && abs(p.outward20W-.004)<1e-13 && isfield(r,'hybridSolvedWallTemperatureK'));end
    assert(abs(sum([p.regions.outward20W])-p.outward20W)<1e-13 && abs(sum([p.regions.outward30W])-p.outward30W)<1e-13);
    web=fig.UserData.webUI;t=tic;while isempty(web.html.UserData.rendered) && toc(t)<35,drawnow;pause(.1);end
    gl=web.html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==10272);
    v=web.evidenceView;d=v.getData();assert(isequal(d.comparisonColorRanges.heatFlux,d.comparisonColorRanges.outerHeat30));
    names={'outerHeat30','outerDeltaHeat','outerDeltaPressure','outerDeltaFirst'};
    values={p.heat30Wm2/1000,p.deltaHeatWm2,p.deltaPressurePa,p.deltaFirstK};
    for j=1:numel(names)
        f=d.fields(strcmp({d.fields.key},names{j}));assert(isequal(f.values(:),values{j}) && all(f.status(:)==1) && contains(f.stateNote,p.scope));
    end
    % 模块二：用实际按钮定位最大热流差，核对探针、整表、同坐标截面和离线页。
    press(fig.UserData.controls.outerHeat);[bound,peak]=max(abs(p.deltaHeatWm2));state=v.getState();
    assert(strcmp(state.fieldKey,'outerDeltaHeat') && state.selectedFaceID==peak && isequal(d.comparisonColorRanges.outerDeltaHeat,[-bound bound]));
    folder=fullfile(out,mode);mkdir(folder);v.exportTo(fullfile(folder,'faces.csv'));
    x=.000172;v.setSectionX(x);t=tic;
    while toc(t)<35
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'outerDeltaHeat') && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'outerDeltaHeat') && abs(cut.x-x)<1e-12);
    path=fullfile(folder,'section.csv');v.exportSectionTo(path);csv=readtable(path,'Encoding','UTF-8');
    assert(height(csv)>0 && max(abs(csv.value-p.deltaHeatWm2(csv.face_id)))<1e-9);
    page=exportAerothermalWeb(r,fullfile(folder,'workbench'));proof.(mode)=struct('page',page,'faces',gl.faces,'webglError',gl.webglError, ...
        'visibleSamples',gl.visiblePixelSamples,'probeFace',peak,'sectionX',x,'cutRows',height(csv),'summary',d.outerDomainComparison);
end
% 模块三：零功率在两域都是物理零；20D 允许而 30D 超域时只撤销该对照。
cfg=s.cfg;cfg.thermistorResultSource='archive';cfg.thermistorRunDirectory='';cfg.thermistorBoundary.netHeatingPowerW=0;
zero=runAeroHeatSingleCase(s.mesh,cfg);o=zero.thermistorOuterDomain;
assert(o.available && all(o.deltaHeatWm2==0) && o.wall20K==293.15 && o.wall30K==293.15);
if p.conductance20WK>p.conductance30WK
    cfg.thermistorBoundary.netHeatingPowerW=6*(p.conductance20WK+p.conductance30WK);
    edge=runAeroHeatSingleCase(s.mesh,cfg);assert(~edge.thermistorOuterDomain.available && all(isfinite(edge.hybridWall.heatFluxWm2)));
    proof.onlyReferenceOutOfRange=true;
end
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,val)
c.Value=val;c.ValueChangedFcn(c,[]);drawnow;
end
