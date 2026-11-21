function check_wall_ui(out,level)
% 模块一：沿生产侧栏准备网格、输入壁温，确认结果绑定本次冻结的边界。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);edit(fig.UserData.controls.leadGrid,level);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
c=fig.UserData.controls;assert(strcmp(c.beadWall.Enable,'on'));
edit(c.beadWall,299.15);press(c.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
r=s.result;web=fig.UserData.webUI;v=web.evidenceView;d=v.getData();
[folder,archiveGuard]=openThermistorArchive(thermistorLeadGridChoice(level)); %#ok<NASGU>
w=readtable(fullfile(folder,'wall.csv'));ref=jsondecode(fileread(fullfile(folder,'summary.json')));
checkWall(r,w,ref,299.15,level);
assert(d.resultCondition.wallUniformValueK==299.15 && contains(d.caseNote,'给定壁温 299.15 K'));
clock=tic;while isempty(web.html.UserData.rendered) && toc(clock)<40,drawnow;pause(.1);end
render=web.html.UserData.rendered;
assert(render.faces==10272 && render.webglError==0 && render.visiblePixelSamples>0);

% 模块二：导出同一快照的压力/热流联合表、首层流体温度、截面和离线三维页。
pair=v.exportWallPairTo(fullfile(out,'wall.csv'));
assert(pair.joint.pressureHeatFluxQualifiedArea==0);
v.selectField('nearWallCellTemperature');v.inspectFace(r.thermistorLead.peakFaceID,[]);
assert(v.getState().selectedFaceID==r.thermistorLead.peakFaceID);
v.exportTo(fullfile(out,'first.csv'));v.selectField('heatFlux');x=.000172;v.setSectionX(x);clock=tic;
while toc(clock)<40
    drawnow;cut=v.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-12);
v.exportSectionTo(fullfile(out,'section.csv'));page=exportAerothermalWeb(r,fullfile(out,'workbench'));
proof=struct('level',level,'page',page,'regionW',[r.thermistorLead.regions.outwardW], ...
    'totalOutwardW',r.thermistorLead.totalOutwardW,'peakFaceID',r.thermistorLead.peakFaceID, ...
    'thermal',r.thermistorLead.thermal,'source',r.evidenceIdentity.sourceFrame,'faces',render.faces, ...
    'visibleSamples',render.visiblePixelSamples,'webglError',render.webglError);

% 模块三：验证区间、零值的单值色标和参照端点，检查既有比较各自的适用温度。
cfg=s.cfg;cfg.thermistorResultSource='archive';cfg.thermistorRunDirectory='';
for wall=[293.15 301.75 305.15]
    cfg.thermistorBoundary.wallTemperatureK=wall;
    trial=runAeroHeatSingleCase(s.mesh,cfg);checkWall(trial,w,ref,wall,level);
    if wall==293.15
        zero=buildAerothermalViewerData(trial);
        assert(zero.uniformFieldValues.heatFlux.value==0 && zero.uniformFieldValues.wallInput.value==wall);
        assert(isempty(trial.thermistorLead.peakFaceID) && isnan(trial.thermistorLead.beadMeanHWm2K) && ...
            all(arrayfun(@(g)isempty(g.peakFaceID),trial.thermistorLead.regions)));
    end
end
for invalid={293.15-eps(293.15),305.15+eps(305.15),NaN,Inf,1i,[299 300]}
    cfg.thermistorBoundary.wallTemperatureK=invalid{1};
    try
        thermistorWorkbenchBoundary(cfg.thermistorBoundary,4,cfg);error('RAH:Test:Accepted','非法壁温通过');
    catch ex,assert(strcmp(ex.identifier,'RAH:Bead:Boundary'),ex.message);end
end
% 模块四：切换净功率清除旧快照，再计算集总壁温并核对两种边界的身份。
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;edit(c.beadCondition,'power');
state=fig.UserData.workflowState;
assert(isempty(state.result) && ~state.boundaryConfirmed && isempty(fig.UserData.webUI) && strcmp(c.beadWall.Enable,'off'));
press(c.confirm);press(fig.UserData.controls.readBead);state=fig.UserData.workflowState;
assert(state.currentStep==5 && ~isempty(state.result),fig.UserData.status.Text);
p=state.result;assert(~isfield(p,'hybridWallTemperatureK') && isfield(p,'hybridSolvedWallTemperatureK') && ...
    abs(p.thermistorLead.totalOutwardW-.004)<1e-13 && ~isfield(p,'thermistorLayerBudget'));
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;edit(c.beadCondition,'published');
assert(isempty(fig.UserData.workflowState.result) && strcmp(c.beadWall.Enable,'on'));
proof.otherWallsK=[293.15 301.75 305.15];proof.invalidRejected=6;proof.powerSolvedWallK=p.hybridSolvedWallTemperatureK(1);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function checkWall(r,ref,s,wall,level)
% 每个区域通过独立面积分验证；分带热量账和光滑配对只在原生端点展示。
scale=(wall-293.15)/12;w=r.hybridWall;g=r.thermistorLead;
assert(isequal(ref.faceID,w.faceID) && isequal(ref.pressurePa,w.pressurePa) && ...
    all(r.hybridWallTemperatureK==wall) && r.hybridCondition.wallTemperatureK==wall && ~isfield(r,'hybridSolvedWallTemperatureK'));
assert(isequal(w.heatFluxWm2,scale*ref.heatFluxWm2) && ...
    max(abs(w.nearWallCellTemperatureK-(293.15+scale*(ref.nearWallCellTemperatureK-293.15))))<1e-12);
assert(abs(sum(ref.areaM2)-3.19665061256752e-6)<1e-17 && numel(g.regions)==7 && ...
    abs(g.totalOutwardW-s.totalOutwardW*scale)<1e-13);
for j=1:7
    use=ref.regionID==j;assert(abs(g.regions(j).outwardW+ref.areaM2(use)'*w.heatFluxWm2(use))<1e-13);
end
use=logical(r.thermistorLeadFlow.reverse);
assert(abs(r.thermistorLeadFlow.summary.reverseReleaseW+ref.areaM2(use)'*w.heatFluxWm2(use))<1e-13);
if scale~=1
    assert(isequaln(s,g.referenceNative) && strcmp(g.thermal.mode,'givenWall') && ...
        ~isfield(g.thermal,'netHeatingPowerW') && isequal(r.evidenceIdentity.sourceFrame.thermalConversion,g.thermal));
end
if level==1
    pair=r.thermistorOuterDomain;assert(pair.wall30K==wall && pair.wall20K==wall && ...
        abs(pair.outward30W-pair.conductance30WK*(wall-293.15))<1e-12);
    grid=r.thermistorLeadGrid;assert(grid.wallK==wall && grid.scale==scale && abs(grid.outwardW-grid.referenceOutwardW*scale)<1e-12);
    assert(isfield(r,'thermistorPaired')==(scale==1));
else
    grid=r.thermistorThreeGrids.summary;assert(all([grid.grids.wallK]==wall) && grid.boundaryValue==wall);
    assert(max(abs([grid.grids.totalOutwardW]-[grid.grids.conductanceWK]*(wall-293.15)))<1e-12);
    assert(isfield(r,'thermistorLayerBudget')==(scale==1));
end
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function edit(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
