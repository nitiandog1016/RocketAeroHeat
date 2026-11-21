function check_wired_ui(out)
% 模块一：从带线模型的网格页选择中段加密，按五步侧栏读取当前冻结结果。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,1);
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(s.mesh.thermistorPrepared.cells==23224,fig.UserData.status.Text);
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
r=s.result;p=r.thermistorPaired;assert(strcmp(p.kind,'wiredRefined') && size(r.F,1)==10272);
web=fig.UserData.webUI;t=tic;while isempty(web.html.UserData.rendered) && toc(t)<35,drawnow;pause(.1);end
gl=web.html.UserData.rendered;assert(gl.webglError==0 && gl.faces==10272 && gl.visiblePixelSamples>0);
v=web.evidenceView;d=v.getData();g=r.thermistorLeadGrid;mask=isfinite(p.heatFluxWm2);assert(nnz(mask)==6144);
fields={'pairedDeltaHeat',p.deltaHeatWm2;'pairedDeltaPressure',p.deltaPressurePa;'pairedReverseChange',p.reverseChange; ...
    'leadGridDeltaHeat',g.deltaHeatWm2;'leadGridDeltaPressure',g.deltaPressurePa;'leadGridDeltaTemperature',g.deltaNearWallK};
for j=1:size(fields,1)
    f=d.fields(strcmp({d.fields.key},fields{j,1}));
    assert(isequaln(f.values(:),fields{j,2}) && isequal(f.status(:)>0,isfinite(fields{j,2})));
end
c=d.comparisonColorRanges;assert(isequal(c.heatFlux,c.leadGridBaseHeat,c.pairedSmoothHeat));
% 模块二：分别用侧栏按钮拾取最大的网距差和几何差，验证两种曲线及全表面导出。
data={fig.UserData.controls.gridHeat,'leadGridDeltaHeat',g.deltaHeatWm2,'grid'; ...
      fig.UserData.controls.pairedHeat,'pairedDeltaHeat',p.deltaHeatWm2,'geometry'};
peaks=zeros(1,2);counts=zeros(1,2);x=.000172;
for j=1:2
    press(data{j,1});key=data{j,2};values=data{j,3};name=data{j,4};
    ids=find(isfinite(values));[bound,k]=max(abs(values(ids)));peaks(j)=ids(k);
    assert(v.getState().selectedFaceID==peaks(j) && strcmp(v.getState().fieldKey,key) && isequal(c.(key),[-bound bound]));
    v.exportTo(fullfile(out,[name '_faces.csv']));
    v.setSectionX(x);t=tic;
    while toc(t)<35
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12);
    file=fullfile(out,[name '_cut.csv']);v.exportSectionTo(file);csv=readtable(file,'Encoding','UTF-8');
    counts(j)=height(csv);assert(counts(j)>0 && max(abs(csv.value-values(csv.face_id)))<1e-9);
end
absent=find(~mask,1);v.inspectFace(absent,[]);assert(v.getState().selectedFaceID==absent);
page=exportAerothermalWeb(r,fullfile(out,'workbench'));
% 模块三：切回旧网，检查已有求解快照和边界确认的失效，并回读原带线结果。
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.leadGrid,0);s=fig.UserData.workflowState;
assert(isempty(s.mesh) && isempty(s.result) && ~s.boundaryConfirmed);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && s.mesh.thermistorPrepared.cells==16264 && strcmp(s.result.thermistorPaired.kind,'wired') && ...
    ~isfield(s.result,'thermistorLeadGrid'));
proof=struct('faces',gl.faces,'webglError',gl.webglError,'visiblePixelSamples',gl.visiblePixelSamples,'commonFaces',nnz(mask), ...
    'sectionX',x,'gridPeak',peaks(1),'geometryPeak',peaks(2),'gridCutRows',counts(1),'geometryCutRows',counts(2), ...
    'sourceFrame',r.evidenceIdentity.sourceFrame.sha256,'coarseFrame',s.result.evidenceIdentity.sourceFrame.sha256, ...
    'smoothFrame',p.summary.sourceFrames.smooth,'meshSHA256',r.evidenceIdentity.meshSHA256,'heatRangeKWm2',c.heatFlux,'page',page);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
function choose(control,value)
control.Value=value;control.ValueChangedFcn(control,[]);drawnow;
end
