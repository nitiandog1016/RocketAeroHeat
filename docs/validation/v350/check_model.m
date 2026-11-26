function check_model(out,mode)
if nargin<2,mode='archive';end
% 模块一：从默认工作台导入实际 STL，逐个角点核对几何、尺度、面序和分区。
evalin('base','main');app=evalin('base','appFigure');cleanup=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';c=app.UserData.controls;cfg=app.UserData.workflowState.cfg;
definition=jsondecode(fileread(fullfile(fileparts(c.file.Value),'geometry.json')));
assert(strcmp(cm011FileSHA256(c.file.Value),definition.fileSHA256) && c.scale.Value==1);
assert(isempty(app.UserData.workflowState.rawMesh));
press(c.load);s=app.UserData.workflowState;m=s.rawMesh;assert(size(m.F,1)==definition.faces);
[original,truth]=buildRocketTestBody(definition.options);
for j=1:3
    corners=original.V(original.F(:,j),:);
    assert(max(abs(m.V(m.F(:,j),:)-double(single(corners))),[],'all')<1e-10);
end
first=1;
for group=definition.regions'
    assert(group.firstFace==first && all(strcmp(truth.faceRegion(first:group.lastFace),group.name)));
    assert(abs(sum(m.area(first:group.lastFace))-group.areaM2)<1e-10);
    first=group.lastFace+1;
end
assert(first==size(m.F,1)+1);
[rawArea,rawVolume]=checkClosed(m);
assert(abs(rawArea-definition.areaM2)<1e-10 && abs(rawVolume-definition.volumeM3)<1e-10);
% 模块二：保留默认朝向和网格方法，验证关闭的壳体经过普通导入流程仍然闭合。
press(app.UserData.controls.confirm);press(app.UserData.controls.generate);
s=app.UserData.workflowState;assert(~isempty(s.mesh) && isempty(s.result),app.UserData.status.Text);
[remeshArea,remeshVolume]=checkClosed(s.mesh);
assert(abs(remeshArea/rawArea-1)<.02 && abs(remeshVolume/rawVolume-1)<.02);
press(app.UserData.workflowButtons(3));s=app.UserData.workflowState;assert(s.currentStep==3 && ~s.busy);
proof=struct('surfaceSHA256',s.modelSource.sha256,'inputFaces',size(m.F,1),'inputVertices',size(m.V,1), ...
    'inputAreaM2',rawArea,'inputVolumeM3',rawVolume,'remeshFaces',size(s.mesh.F,1),'remeshAreaM2',remeshArea, ...
    'remeshVolumeM3',remeshVolume,'closedConnectedOutward',true,'verifiedRegions',numel(definition.regions), ...
    'stepBeforeSolve',s.currentStep);
% 模块三：固定低速工况一次确认。对全部原面核对显示值、符号和总热功率。
press(app.UserData.controls.confirm);start=tic;
if strcmp(mode,'archive'),press(app.UserData.controls.readLauncher);else,press(app.UserData.controls.run);end
s=app.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),app.UserData.status.Text);
wall=s.result.hybridWall;report=app.UserData.caseReport;
view=app.UserData.webUI;timer=tic;
while isempty(view.html.UserData.rendered) && toc(timer)<45,drawnow;pause(.1);end
gl=view.html.UserData.rendered;assert(gl.faces==definition.faces && gl.visiblePixelSamples>0 && gl.webglError==0);
proof.resultWaitSeconds=toc(start);
data=buildAerothermalViewerData(s.result);flux=data.fields(strcmp({data.fields.key},'heatFlux'));
assert(strcmp(flux.unit,'W/m²') && isequal(flux.values(:),wall.heatFluxWm2(:)) && all(flux.status==1));
assert(data.wallCoverage.fullSurfaceComputed && abs(data.wallCoverage.wholeSurfaceHeatPowerW-report.wholeSurfaceHeatPowerW)<1e-12);
view.evidenceTab.Parent.SelectedTab=view.evidenceTab;ev=view.evidenceView;ev.selectField('heatFlux');
for g=definition.regions'
    face=g.firstFace;ev.inspectFace(face,[]);probe=ev.getState();assert(probe.selectedFaceID==face);
end
ev.setSectionX(.7);timer=tic;
while toc(timer)<40
    drawnow;cut=ev.getSection();if isstruct(cut) && isfield(cut,'csv') && abs(cut.x-.7)<1e-10,break;end
    pause(.1);
end
assert(isstruct(cut) && isfield(cut,'csv') && abs(cut.x-.7)<1e-10);
rec=ev.exportWallPairTo(fullfile(fileparts(out),'wall_pair.csv'));assert(rec.rows==definition.faces);
exportAerothermalWeb(s.result,fullfile(fileparts(out),'workbench'));
proof.mode=mode;proof.finalResultPublished=true;proof.finalFaces=report.jointComputedFaceCount;
proof.resultAreaM2=report.totalArea;proof.outwardW=-report.wholeSurfaceHeatPowerW;
proof.probedRegions=numel(definition.regions);proof.sectionXM=cut.x;proof.exportedRows=rec.rows;
proof.visiblePixels=gl.visiblePixelSamples;proof.runDirectory=report.runDirectory;
writelines(string(jsonencode(proof,'PrettyPrint',true)),out,'Encoding','UTF-8');
end

function [area,volume]=checkClosed(m)
% 模块四：由三角形独立核对边的使用次数、方向、连通性和有向体积。
f=m.F;v=m.V;edges=[f(:,[1 2]);f(:,[2 3]);f(:,[3 1])];
[pairs,~,ids]=unique(sort(edges,2),'rows');
assert(all(accumarray(ids,1)==2) && all(accumarray(ids,sign(edges(:,2)-edges(:,1)))==0));
assert(max(conncomp(graph(pairs(:,1),pairs(:,2),[],size(v,1))))==1);
crosses=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
areas=vecnorm(crosses,2,2)/2;assert(all(areas>0));area=sum(areas);
volume=sum(dot(v(f(:,1),:),crosses,2))/6;assert(volume>0);
end

function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
