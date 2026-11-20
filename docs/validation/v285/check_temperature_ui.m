function check_temperature_ui(out,cells)
% 两档网格各自从工作台输入三个壁温。核对主场、截面、实际色标、PNG 与成对 CSV。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.beadModel);
if cells==4,choose(fig.UserData.controls.beadCells,4);end
press(fig.UserData.controls.generate);
assert(~isempty(fig.UserData.workflowState.mesh),fig.UserData.status.Text);
saved=struct;results=struct;
% 使用原配置回放并留存独立对照。
for temp=[305.15 299.15 293.15]
    press(fig.UserData.workflowButtons(3));
    c=fig.UserData.controls;
    if strcmp(c.beadCondition.Value,'none'),choose(c.beadCondition,'published');end
    if c.beadWall.Value~=temp
        choose(c.beadWall,temp);s=fig.UserData.workflowState;
        assert(isempty(s.result) && isempty(fig.UserData.webUI) && ~s.steps(5).available && ~s.boundaryConfirmed);
    end
    press(c.confirm);s=fig.UserData.workflowState;
    assert(s.currentStep==4 && s.cfg.thermistorBoundary.wallTemperatureK==temp);
    press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
    assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
    r=s.result;rp=fig.UserData.caseReport;view=fig.UserData.webUI;
    w=r.hybridWall;n=height(w);assert(n==5216 && rp.fullSurfaceComputed);
    deltaT=temp-293.15;scale=deltaT/12;
    if scale==1
        saved=r;[referenceFolder,cleanReference]=openThermistorArchive(cells);
        reference=readtable(fullfile(referenceFolder,'wall.csv'));clear cleanReference;
        assert(isequal(w.faceID,reference.faceID) && isequal(w.pressurePa,reference.pressurePa) && ...
            isequal(w.heatFluxWm2,reference.heatFluxWm2) && isequal(w.nearWallCellTemperatureK,reference.nearWallCellTemperatureK));
    end
    assert(isequal(w.pressurePa,saved.hybridWall.pressurePa) && isequal(w.heatFluxWm2,scale*saved.hybridWall.heatFluxWm2));
    assert(max(abs(w.nearWallCellTemperatureK-(293.15+scale*(saved.hybridWall.nearWallCellTemperatureK-293.15))))<1e-12);
    assert(all(r.hybridWallTemperatureK==temp) && r.hybridCondition.wallTemperatureK==temp);
    assert(isfield(r.thermistorBead,'experiment')==(scale==1) && isfield(r.thermistorBead,'relativeDifference')==(scale==1));
    if cells==4
        assert(isequal(r.beadAngularDelta.pressurePa,saved.beadAngularDelta.pressurePa));
        assert(max(abs(r.beadAngularDelta.heatFluxWm2-scale*saved.beadAngularDelta.heatFluxWm2))<1e-9);
    end
    if scale~=1
        assert(r.evidenceIdentity.sourceFrame.thermalConversion.wallTemperatureK==temp && ~isfield(rp,'native') && isfield(rp,'referenceNative'));
        if cells==2
            if scale==0,name='zero';else,name='half';end
            other=readtable(fullfile(root,'outputs','v285_validation','direct01',name,'direct_wall.csv'));
            assert(isequal(w.faceID,other.faceID) && max(abs(w.heatFluxWm2-other.heatFluxWm2))<.012 && ...
                max(abs(w.nearWallCellTemperatureK-other.nearWallCellTemperatureK))<1e-5 && max(abs(w.pressurePa-other.pressurePa))<.001);
        end
    end
    started=tic;while isempty(view.html.UserData.rendered) && toc(started)<35,drawnow;pause(.1);end
    gl=view.html.UserData.rendered;assert(~isempty(gl) && gl.faces==n && gl.webglError==0 && gl.visiblePixelSamples>0);
    v=view.evidenceView;d=v.getData();
    q=d.fields(strcmp({d.fields.key},'heatFlux'));ts=d.fields(strcmp({d.fields.key},'wallSolved'));
    assert(all(q.status==1) && max(abs(q.values*1000-w.heatFluxWm2))<1e-9 && all(ts.status==0 & isnan(ts.values)));
    assert(contains(q.stateNote,'图 7')==(scale==1));
    v.selectField('heatFlux');v.inspectFace(5204,[]);
    assert(v.getState().selectedFaceID==5204 && q.values(5204)==w.heatFluxWm2(5204)/1000);
    target=fullfile(out,sprintf('wall_%.2f.csv',temp));receipt=v.exportWallPairTo(target);
    wt=readtable(target);assert(isequal(wt.face_id,w.faceID+1) && max(abs(wt.heat_flux*1000-w.heatFluxWm2))<1e-9);
    assert(receipt.joint.pressureHeatFluxQualifiedArea==0 && receipt.identity.sourceTraceable);
    row=struct('wallK',temp,'powerW',rp.wholeSurfaceHeatPowerW,'areaM2',rp.jointComputedAreaM2, ...
        'faces',n,'computedZeroFaces',nnz(q.values==0 & q.status==1),'legend',gl.legend,'conditionSHA256',receipt.identity.conditionSHA256);
    if scale~=1
        imageBytes=matlab.net.base64decode(extractAfter(string(gl.png),'base64,'));
        assert(isequal(reshape(imageBytes(1:4),1,[]),uint8([137 80 78 71])));
        if scale==0,pngName='zero_heat.png';else,pngName='half_heat.png';end
        fid=fopen(fullfile(out,pngName),'wb');assert(fid>=0);fwrite(fid,imageBytes,'uint8');fclose(fid);
        row.pngBytes=numel(imageBytes);
    end
    if scale==0
        fid=fopen(fullfile(out,'zero_legend.json'),'w','n','UTF-8');assert(fid>=0);
        fprintf(fid,'%s\n',jsonencode(struct('legend',gl.legend,'exportRange',gl.exportRange,'colors',gl.exportColors)));fclose(fid);
        assert(all(w.nearWallCellTemperatureK==temp) && all(w.heatFluxWm2==0) && rp.wholeSurfaceHeatPowerW==0 && isnan(r.thermistorBead.meanH));
        assert(strcmp(gl.legend.low,'0') && isempty(gl.legend.high) && contains(gl.legend.note,'均为 0') && ~contains(gl.legend.swatch,'gradient'));
        assert(isequal(gl.exportRange(:),[0;0]) && isequal(gl.exportColors(1,:),gl.exportColors(2,:)));
        view.evidenceTab.Parent.SelectedTab=view.evidenceTab;v.setSectionX(.000172);start=tic;
        while toc(start)<30
            drawnow;sec=v.getSection();
            if isstruct(sec) && isfield(sec,'csv') && isfield(sec,'x') && abs(sec.x-.000172)<1e-11 && strcmp(sec.fieldKey,'heatFlux'),break;end
            pause(.1);
        end
        sectionPath=fullfile(out,'zero_section.csv');v.exportSectionTo(sectionPath);
        cut=readtable(sectionPath);assert(height(cut)>0 && all(cut.value==0) && all(string(cut.availability)=="研究计算"));
        row.sectionRows=height(cut);
        exportAerothermalWeb(r,fullfile(out,'zero_workbench'));
    else
        assert(isfinite(r.thermistorBead.meanH) && ~isempty(gl.legend.high));
        if scale~=1
            assert(gl.exportRange(1)<gl.exportRange(2));
            exportAerothermalWeb(r,fullfile(out,'half_workbench'));
        end
    end
    results.(sprintf('dT%d',deltaT))=row;
end
assert(numel(unique(string({results.dT0.conditionSHA256,results.dT6.conditionSHA256,results.dT12.conditionSHA256})))==3);
for invalid={293.14,305.16,NaN,Inf,[299 300],complex(300,1)}
    b=thermistorWorkbenchBoundary();b.wallTemperatureK=invalid{1};
    try,thermistorWorkbenchBoundary(b);error('RAH:Test:Boundary','接受了无效壁温。');
    catch ex,assert(strcmp(ex.identifier,'RAH:Bead:Boundary'));end
end
b=thermistorWorkbenchBoundary();b.wallTemperatureK=single(300);b=thermistorWorkbenchBoundary(b);assert(isequal(b.wallTemperatureK,300));
results.boundaryRangeChecked=true;
fid=fopen(fullfile(out,'view.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(results,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,v)
c.Value=v;c.ValueChangedFcn(c,[]);drawnow;
end
