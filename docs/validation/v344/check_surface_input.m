function check_surface_input(out)
% 模块一：从算例按钮确认原面与内外壁温，用真实工作台输出工况与全表。
dialog=fullfile(out,'dialog');mkdir(dialog);
writeDialog(fullfile(dialog,'uiputfile.m'),"function [f,p]=uiputfile(~,~,f),p=getappdata(0,'rfzTestOut');end");
writeDialog(fullfile(dialog,'uigetfile.m'),"function [f,p]=uigetfile(varargin),p=getappdata(0,'rfzTestOut');f=getappdata(0,'rfzTestFile');end");
setappdata(0,'rfzTestOut',out);setappdata(0,'rfzTestFile','RFZ-REC-input.json');
addpath(dialog,'-begin');clean=onCleanup(@()restoreDialogs(dialog)); %#ok<NASGU>
fig=rocketAeroHeatWorkbench();first=onCleanup(@()closeFig(fig)); %#ok<NASGU>
press(fig.UserData.controls.rfzModel);s=fig.UserData.workflowState;raw=s.rawMesh;
assert(isfield(raw,'rfzBoundary'),'%s',fig.UserData.status.Text);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.generate);
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);
canonical=fig.UserData.workflowState;
assert(isequal(canonical.mesh.V,raw.V) && isequal(canonical.mesh.F,raw.F));
constant=waitFrame(fig.UserData.geometryView.html);
assert(isequal(constant.fieldRange(:),[300;300]) && contains(constant.givenWallLegend,'全面 300 K'));
press(fig.UserData.controls.rfzSave);
copyfile(fullfile(out,'RFZ-REC-input.json'),fullfile(out,'full_cad_input.json'));
press(fig.UserData.workflowButtons(1));
c=fig.UserData.controls;c.rx.Value=7;press(c.rotate);press(fig.UserData.controls.confirm);
c=fig.UserData.controls;assert(strcmp(c.meshBackend.Value,'rfzSurface'));press(c.generate);
s=fig.UserData.workflowState;m=s.mesh;
assert(isequal(m.F,raw.F) && isequal(m.V,s.orientedMesh.V) && ~isequal(m.V,raw.V) && ~s.boundaryConfirmed);
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);
s=fig.UserData.workflowState;assert(s.currentStep==4 && s.boundaryConfirmed,'%s',fig.UserData.status.Text);
assert(all(s.cfg.heating.wallTemperatureByFace==300) && strcmp(s.cfg.rfzSurfaceInput.cavityMode,'matchExterior'));
defaultSHA=s.cfg.rfzSurfaceInput.wallSHA256;
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;
c.rfzMode.Value='specified';c.rfzMode.ValueChangedFcn(c.rfzMode,[]);drawnow;
assert(~fig.UserData.workflowState.boundaryConfirmed && ~fig.UserData.workflowState.steps(4).available);
assert(~isfield(fig.UserData.geometryView.html.UserData.payload,'faceField'));
c.innerK.Value=0;c.innerK.ValueChangedFcn(c.innerK,[]);drawnow;
invalid=fig.UserData.workflowState;assert(~invalid.boundaryConfirmed && ~isfield(invalid.cfg,'rfzSurfaceInput') && ...
    ~isfield(invalid.cfg.heating,'wallTemperatureByFace'));
c.innerK.Value=315;c.innerK.ValueChangedFcn(c.innerK,[]);press(c.confirm);
s=fig.UserData.workflowState;assert(s.currentStep==4 && s.boundaryConfirmed,'%s',fig.UserData.status.Text);
c=fig.UserData.controls;press(c.rfzSave);press(c.rfzWallCSV);
input=jsondecode(fileread(fullfile(out,'RFZ-REC-input.json')));
tableA=readtable(fullfile(out,'RFZ-REC-wall-input.csv'),'TextType','string');
g=raw.rfzBoundary.regions;ids=repelem(g.cadFace,g.lastFace-g.firstFace+1);inside=g.referenceRole(ids)=="封盖内侧 CAD 壁";
assert(isequal(tableA.faceID,(1:130044)') && isequal(tableA.cadFace,ids) && ...
    all(tableA.givenWallK(~inside)==300) && all(tableA.givenWallK(inside)==315) && nnz(inside)==14982);
assert(all(tableA.temperatureSource(inside)=="用户单独给定内壁温度") && all(tableA.temperatureSource(~inside)=="论文给定外壁温度"));
assert(max(abs(tableA.areaM2-m.area(:)))<1e-12 && abs(sum(tableA.areaM2)-740.7290034031088)<1e-8);
assert(isequal(s.cfg.heating.wallTemperatureByFace,tableA.givenWallK) && ~strcmp(input.wallSHA256,defaultSHA));
assert(isempty(s.result) && ~s.steps(5).available);
% 同一原面上的给定温度与稿件记录进入冻结输入，应用层未执行流场求解。
frozen=aeroHeatWorkflowState('beginSolve',s);
assert(isequal(frozen.frozenInputs.cfg.rfzSurfaceInput,input) && isequal(frozen.frozenInputs.mesh.F,raw.F));
assert(~isfield(frozen.frozenInputs.cfg,'flight') && isequal(frozen.frozenInputs.freestream,raw.rfzBoundary.condition.publishedCondition));
bad=s;bad.cfg.heating.wallTemperatureByFace(1)=999;
rejected=false;try,aeroHeatWorkflowState('beginSolve',bad);catch ex,rejected=strcmp(ex.identifier,'RAH:RFZ:Boundary');end
assert(rejected);
closeFig(fig);

% 模块二：全新窗口从模型浏览入口重开，检查旋转、原面、探针和损坏输入的隔离。
fig=rocketAeroHeatWorkbench();second=onCleanup(@()closeFig(fig)); %#ok<NASGU>
press(fig.UserData.controls.browse);restored=fig.UserData.workflowState;
assert(restored.currentStep==4 && restored.boundaryConfirmed,'%s',fig.UserData.status.Text);
assert(isequal(restored.mesh.V,m.V) && isequal(restored.mesh.F,raw.F) && ...
    isequal(restored.cfg.rfzSurfaceInput,input) && isequal(restored.cfg.heating.wallTemperatureByFace,tableA.givenWallK));
assert(isequal(restored.freestream,raw.rfzBoundary.condition.publishedCondition) && isempty(restored.result));
html=fig.UserData.geometryView.html;gl=waitFrame(html);assert(gl.faces==130044);
assert(isequal(gl.fieldRange(:),[300;315]) && contains(gl.givenWallLegend,'300—315'));
assert(isequal(html.UserData.payload.faceField.values(:),tableA.givenWallK));
id=find(inside,1);html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',id-1);
html.DataChangedFcn(html,[]);assert(contains(fig.UserData.status.Text,'315 K（用户给定内壁）'));
bad=input;bad.cavityWallK=316;writeDialog(fullfile(out,'bad.json'),jsonencode(bad));setappdata(0,'rfzTestFile','bad.json');
press(fig.UserData.controls.rfzReopen);after=fig.UserData.workflowState;
assert(isequaln(after,restored) && contains(fig.UserData.status.Text,'不一致'));
record=struct('domain',input.domain,'faces',height(tableA),'totalAreaM2',sum(tableA.areaM2), ...
    'innerFaces',nnz(inside),'innerAreaM2',sum(tableA.areaM2(inside)),'defaultWallK',300,'testInnerK',315, ...
    'meshSHA256',input.meshSHA256,'defaultWallSHA256',defaultSHA,'testWallSHA256',input.wallSHA256, ...
    'csvSHA256',cm011FileSHA256(fullfile(out,'RFZ-REC-wall-input.csv')),'probeFace',id,'pixels',gl.visiblePixelSamples, ...
    'uniformLegend',constant.givenWallLegend,'splitLegend',gl.givenWallLegend, ...
    'fullReopen',true,'badInputIsolated',true,'confirmedInputFaces',height(tableA),'computedFaces',0);
writeDialog(fullfile(out,'workbench.json'),jsonencode(record,'PrettyPrint',true));
end
function writeDialog(file,content)
f=fopen(file,'w','n','UTF-8');assert(f>=0);c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',content);
end
function restoreDialogs(folder)
rmpath(folder);rmappdata(0,'rfzTestOut');rmappdata(0,'rfzTestFile');
end
function closeFig(fig)
if isgraphics(fig),delete(fig);end
end
function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end
function gl=waitFrame(html)
t=tic;
while toc(t)<55
    drawnow;s=html.UserData;assert(isempty(s.error),'%s',jsonencode(s.error));gl=s.rendered;
    if ~isempty(gl) && gl.revision==s.revision && gl.webglError==0 && gl.visiblePixelSamples>0,return;end
    pause(.1);
end
error('未收到当前原面几何的绘制回执。');
end
