function check_ui(out,mode,otherUI)
% 模块一：经真实侧栏选用带线试件，检查模型、体网、边界与计算来源的顺序。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);s=fig.UserData.workflowState;
assert(s.currentStep==2 && strcmp(s.cfg.thermistorVariant,'leads'));
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh),fig.UserData.status.Text);p=s.mesh.thermistorPrepared;
assert(p.cells==16264 && p.nodes==12231 && p.wallFaces==432 && isempty(s.result));
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);
s=fig.UserData.workflowState;assert(s.currentStep==4 && s.freestream.speed==3 && s.cfg.thermistorBoundary.wallTemperatureK==305.15);
if strcmp(mode,'solve'),press(fig.UserData.controls.run);else,press(fig.UserData.controls.readBead);end
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
bead=s.result;w=bead.hybridWall;report=fig.UserData.caseReport;
assert(height(w)==10272 && report.fullSurfaceComputed && bead.evidenceIdentity.reusedArchive==strcmp(mode,'archive'));
if strcmp(mode,'solve')
    assert(strcmp(report.outputDirectory,bead.evidenceIdentity.runDirectory) && ...
        strcmp(bead.evidenceIdentity.sourceFrame.sha256,cm011FileSHA256(fullfile(report.outputDirectory,'restart.csv'))));
    wall=readtable(fullfile(report.outputDirectory,'surface','wall.csv'));
else
    [folder,guard]=openThermistorArchive('leads3'); %#ok<NASGU>
    wall=readtable(fullfile(folder,'wall.csv'));
end
assert(isequal(w.faceID,wall.faceID) && max(abs(w.heatFluxWm2-wall.heatFluxWm2))<1e-9 && ...
    max(abs(w.pressurePa-wall.pressurePa))<1e-9 && max(abs(w.nearWallCellTemperatureK-wall.nearWallCellTemperatureK))<1e-9);

% 模块二：核对七个区域的三维探针、逐面下载和独立离线页面的实际数据。
ui=fig.UserData.webUI;t=tic;
while isempty(ui.html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(~isempty(gl) && gl.faces==height(w) && gl.webglError==0 && gl.visiblePixelSamples>0);
view=ui.evidenceView;view.selectField('heatFlux');data=view.getData();
regions=bead.thermistorLead.regions;peaks=[regions.peakFaceID];
for j=1:numel(regions)
    view.inspectFace(peaks(j),[]);assert(view.getState().selectedFaceID==peaks(j));
    assert(strcmp(string(data.region(peaks(j))),string(regions(j).name)));
end
view.exportWallPairTo(fullfile(out,'pairs.csv'));
flow=bead.thermistorLeadFlow;view.selectField('leadReverse');
reverse=data.fields(strcmp({data.fields.key},'leadReverse'));
assert(isequal(reverse.values(:),flow.reverse) && nnz(reverse.values)==flow.summary.reverseFaces);
view.inspectFace(find(reverse.values==1,1),[]);view.exportTo(fullfile(out,'reverse.csv'));
if nargin==3
    % 单独以实际侧栏新算产生的终场回归新的诊断入口，复用已通过的原生推进。
    other=jsondecode(fileread(otherUI));assert(strcmp(other.mode,'solve') && strcmp(other.sourceFrame.sha256,bead.evidenceIdentity.sourceFrame.sha256));
    checked=attachThermistorLeadNearWall(bead,other.runDirectory,out);
    assert(isequal(checked.thermistorLeadFlow.reverse,flow.reverse) && ...
        isequal(checked.thermistorLeadFlow.downstreamCf,flow.downstreamCf) && ...
        isequal(checked.thermistorLeadFlow.estimatedYPlus,flow.estimatedYPlus));
end
view.selectField('heatFlux');view.inspectFace(bead.thermistorLead.peakFaceID,[]);
page=exportAerothermalWeb(bead,fullfile(out,'workbench'));assert(isfile(page));
proof=struct('mode',mode,'faces',height(w),'areaM2',sum(wall.areaM2),'regionPeaks',peaks,'sourceFrame',bead.evidenceIdentity.sourceFrame, ...
    'runDirectory',report.outputDirectory,'wholePowerIntoW',report.wholeSurfaceHeatPowerW,'page',page,'webglError',gl.webglError, ...
    'reverseFaces',flow.summary.reverseFaces,'freshDiagnosticsCompared',nargin==3);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end
