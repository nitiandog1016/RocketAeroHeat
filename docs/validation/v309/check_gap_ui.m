function check_gap_ui(out,level)
% 模块一：走实际侧栏输入 299.15 K，核查结果页、给定壁温和热流说明。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,level);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
c=fig.UserData.controls;choose(c.beadWall,299.15);press(c.confirm);press(fig.UserData.controls.readBead);
state=fig.UserData.workflowState;assert(state.currentStep==5 && ~isempty(state.result),fig.UserData.status.Text);
r=state.result;view=fig.UserData.webUI.evidenceView;data=view.getData();
heat=data.fields(strcmp({data.fields.key},'heatFlux'));obs=r.thermistorLead.beadMeanComparison;
assert(strcmp(obs.status,'noMatchingPoint') && obs.deltaTK==6 && isnan(obs.relativeDifference) && isnan(obs.measuredHWm2K));
assert(contains(heat.stateNote,'温差 6 K 处无匹配测点') && contains(heat.stateNote,sprintf('%.2f W/(m²·K)',obs.predictedHWm2K)));
assert(abs(obs.predictedHWm2K-r.thermistorLead.beadAndJointsW/(6*3.094e-6))<1e-10 && r.thermistorLead.jointExperimentVerifiedM2==0);
view.selectField('heatFlux');view.inspectFace(207,[]);assert(view.getState().selectedFaceID==207);
receipt=view.exportTo(fullfile(out,'heat.csv'));assert(receipt.rows==10272);
t=tic;html=fig.UserData.webUI.html;while isempty(html.UserData.rendered) && toc(t)<30,drawnow;pause(.1);end
gl=html.UserData.rendered;assert(gl.faces==10272 && gl.webglError==0 && gl.visiblePixelSamples>0);
page=exportAerothermalWeb(r,fullfile(out,'workbench'));

% 模块二：原端点、空档外的给定温差和同温差的功率模式各保留自己的证据。
cfg=state.cfg;cfg.thermistorRunDirectory='';cfg.thermistorResultSource='archive';
for wall=[296.15 301.15 305.15]
    cfg.thermistorBoundary.wallTemperatureK=wall;
    [trial,~,report]=runAeroHeatSingleCase(state.mesh,cfg);
    hasGap=isfield(trial.thermistorLead,'beadMeanComparison');
    assert(hasGap==(wall==301.15) && isfield(report,'beadMeanComparison')==hasGap);
    if hasGap,assert(report.beadMeanComparison.deltaTK==8);end
    if wall==305.15,assert(abs(trial.thermistorLead.experimentMeanHWm2K-218.35690235690237)<1e-9);end
end
cfg.thermistorBoundary=rmfield(cfg.thermistorBoundary,'wallTemperatureK');
cfg.thermistorBoundary.netHeatingPowerW=r.thermistorLead.referenceNative.totalOutwardW/2;
trial=runAeroHeatSingleCase(state.mesh,cfg);
assert(abs(trial.hybridSolvedWallTemperatureK(1)-299.15)<1e-10 && ~isfield(trial.thermistorLead,'beadMeanComparison'));
proof=struct('gridLevel',level,'page',page,'faces',gl.faces,'visibleSamples',gl.visiblePixelSamples, ...
    'webglError',gl.webglError,'heatNote',heat.stateNote,'chosenFace',207,'comparison',obs, ...
    'wholeOutwardW',r.thermistorLead.totalOutwardW,'sourceFrame',r.evidenceIdentity.sourceFrame.sha256, ...
    'otherGivenWallsK',[296.15 301.15 305.15],'solvedWallAtHalfPowerK',trial.hybridSolvedWallTemperatureK(1));
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
