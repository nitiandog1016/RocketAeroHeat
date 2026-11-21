function check_mean_ui(out,level)
% 模块一：从实际侧栏输入 299.15 K，在结果页读出方法表、热流探针与导出说明。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,level);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
c=fig.UserData.controls;choose(c.beadWall,299.15);press(c.confirm);press(fig.UserData.controls.readBead);
state=fig.UserData.workflowState;assert(state.currentStep==5 && ~isempty(state.result),fig.UserData.status.Text);
r=state.result;m=r.thermistorLead.meanConvection;tab=fig.UserData.controls.beadMeanMethods;
assert(m.available && m.filmTemperatureK==296.15 && strcmp(tab.Visible,'on'));
shownRows=tab.Data;
assert(isequal(tab.Data(:,2:3),{sprintf('%.2f',m.meanHWm2K),sprintf('%.4f',1000*m.outwardW); ...
    sprintf('%.2f',m.surfaceMeanHWm2K),sprintf('%.4f',1000*r.thermistorLead.beadAndJointsW)}));
view=fig.UserData.webUI.evidenceView;d=view.getData();heat=d.fields(strcmp({d.fields.key},'heatFlux'));
assert(contains(heat.stateNote,'温差 6 K 处无匹配测点') && contains(heat.stateNote,'表 1 极限 Nu'));
view.selectField('heatFlux');view.inspectFace(207,[]);assert(view.getState().selectedFaceID==207);
receipt=view.exportTo(fullfile(out,'heat.csv'));assert(receipt.rows==10272);
html=fig.UserData.webUI.html;timer=tic;
while isempty(html.UserData.rendered) && toc(timer)<30,drawnow;pause(.1);end
gl=html.UserData.rendered;assert(gl.faces==10272 && gl.webglError==0 && gl.visiblePixelSamples>0);
page=exportAerothermalWeb(r,fullfile(out,'workbench'));

% 模块二：用同一生产入口检查零温差与 12 K 端点，并确认旧网格比较仍然可用。
cfg=state.cfg;cfg.thermistorRunDirectory='';cfg.thermistorResultSource='archive';
cfg.thermistorBoundary.wallTemperatureK=293.15;
[zero,~,report]=runAeroHeatSingleCase(state.mesh,cfg);zm=report.meanConvection;
assert(~zm.available && contains(zm.reason,'正的') && isnan(zm.meanHWm2K) && isnan(zm.surfaceMeanHWm2K));
assert(zm.surfaceOutwardW==0 && all(zero.hybridWall.heatFluxWm2==0));
cfg.thermistorBoundary.wallTemperatureK=305.15;
[endpoint,~,report]=runAeroHeatSingleCase(state.mesh,cfg);ep=report.meanConvection;
assert(ep.available && ep.filmTemperatureK==299.15 && abs(endpoint.thermistorLead.experimentMeanHWm2K-218.35690235690237)<1e-9);
assert(isfield(endpoint,'thermistorLayerBudget')==(level==2));

% 模块三：切换净功率清除旧表。测试工程方法对低 Re、超温和错误试件的原因。
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;choose(c.beadCondition,'power');
assert(isempty(fig.UserData.workflowState.result));press(c.confirm);press(fig.UserData.controls.readBead);
state=fig.UserData.workflowState;assert(~isempty(state.result) && ~isfield(state.result.thermistorLead,'meanConvection'));
assert(~isfield(fig.UserData.controls,'beadMeanMethods') && isfield(state.result,'hybridSolvedWallTemperatureK'));
[folder,archiveGuard]=openThermistorArchive(thermistorLeadGridChoice(level)); %#ok<NASGU>
spec=jsondecode(fileread(fullfile(folder,'case.json')));changed=spec;changed.speedMps=.1;
low=thermistorBeadMeanConvection(changed,299.15);assert(~low.available && contains(low.reason,'Re_L'));
outside=thermistorBeadMeanConvection(spec,5000);assert(~outside.available && contains(outside.reason,'空气物性'));
changed=spec;changed.caseKind='other';wrong=thermistorBeadMeanConvection(changed,299.15);assert(~wrong.available && contains(wrong.reason,'表 1'));
proof=struct('level',level,'page',page,'comparison',m,'shownRows',{shownRows}, ...
    'webglError',gl.webglError,'visibleSamples',gl.visiblePixelSamples,'probeFace',207,'faces',gl.faces, ...
    'sourceFrame',r.evidenceIdentity.sourceFrame.sha256,'endpoint',ep,'rejections',{{zm.reason,low.reason,outside.reason,wrong.reason}});
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
