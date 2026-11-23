function check_import(out,mode)
% 模块一：从普通文件输入，确认坐标并选择体网。测试真值只在验收端读取。
root=fileparts(which('rocketAeroHeatWorkbench'));fixture=fullfile(root,'outputs','v319_validation','fixture');
path=fullfile(fixture,'probe-mm.stl');order=readmatrix(fullfile(fixture,'reference_ids.txt'));
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=path;c.scale.Value=.001;press(c.load);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'thermistorB');
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh),'%s',fig.UserData.status.Text);b=s.mesh.thermistorImport;p=s.mesh.thermistorPrepared;
assert(isequal(b.referenceFaceID,order) && ~p.reusedGrid && p.cells==558048 && p.nodes==316821);
assert(isequal(s.mesh.F,s.rawMesh.F) && isequal(s.mesh.V,s.rawMesh.V));
press(fig.UserData.controls.nextBoundary);
if strcmp(mode,'native')
    choose(fig.UserData.controls.wireBoundary,'longWire');actions={'native'};
else
    choose(fig.UserData.controls.wireBoundary,'resistance');actions={'solve','archive','uniform'};
end
press(fig.UserData.controls.confirm);
proof=struct('meshSHA256',p.meshSHA256,'prepared',p.directory,'source',b.source,'maxVertexErrorM',b.maxVertexErrorM, ...
    'maxRelativeAreaError',b.maxRelativeAreaError,'minNormalCosine',b.minNormalCosine);
for j=1:numel(actions)
    action=actions{j};dest=fullfile(out,action);mkdir(dest);
    if strcmp(action,'native')
        fig.UserData.controls.wireSeed.Value=fullfile(root,'outputs','v316_validation','run01');
    elseif strcmp(action,'uniform')
        press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.wireBoundary,'uniform');
        assert(isempty(fig.UserData.workflowState.result));press(fig.UserData.controls.confirm);
    end
    if strcmp(action,'archive'),press(fig.UserData.controls.readBead);else,press(fig.UserData.controls.run);end
    s=fig.UserData.workflowState;
    assert(s.currentStep==5 && ~isempty(s.result),'%s',fig.UserData.status.Text);
    % 模块二：从实际渲染器收到像素和截线回执，在七个区域及热峰拾取 STL 面号。
    res=s.result;w=res.hybridWall;assert(isequal(res.F,s.rawMesh.F) && isequal(res.V,s.rawMesh.V));
    report=fig.UserData.caseReport;ui=fig.UserData.webUI;v=ui.evidenceView;
    timer=tic;while isempty(ui.html.UserData.rendered) && toc(timer)<40,drawnow;pause(.1);end
    gl=ui.html.UserData.rendered;assert(gl.webglError==0 && gl.faces==10272 && gl.visiblePixelSamples>0);
    ids=zeros(7,1);v.selectField('heatFlux');
    for k=1:7
        ids(k)=find(strcmp(res.localGeometry.region,report.regions(k).name),1);
        v.inspectFace(ids(k),[]);u=v.getState();assert(u.selectedFaceID==ids(k));
    end
    [~,peak]=min(w.heatFluxWm2);assert(peak==res.thermistorSidewind.peakFaceID);
    v.inspectFace(peak,[]);u=v.getState();assert(u.selectedFaceID==peak);
    wall=v.exportWallPairTo(fullfile(dest,'wall.csv'));assert(wall.rows==numel(order) && wall.joint.pressureHeatFluxQualifiedArea==0);
    x=.000172;v.setSectionX(x);timer=tic;
    while toc(timer)<40
        drawnow;cut=v.getSection();
        if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-12);
    v.exportSectionTo(fullfile(dest,'section.csv'));
    page=exportAerothermalWeb(res,fullfile(dest,'workbench'));
    identity=res.evidenceIdentity;
    proof.(action)=struct('page',page,'frame',identity.sourceFrame,'directory',identity.runDirectory,'archiveSource',identity.reusedArchive, ...
        'regions',report.regions,'peakFaceID',peak,'referencePeak',order(peak),'peakFluxWm2',w.heatFluxWm2(peak), ...
        'areaM2',report.totalArea,'inwardW',report.wholeSurfaceHeatPowerW,'faces',gl.faces,'webGLError',gl.webglError, ...
        'pixels',gl.visiblePixelSamples,'probeIDs',ids,'probeReferenceIDs',order(ids));
    if j<numel(actions),press(fig.UserData.workflowButtons(4));end
end
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.meshBackend,'engineering');
s=fig.UserData.workflowState;assert(isempty(s.mesh) && isempty(s.result) && ~s.boundaryConfirmed);
press(fig.UserData.workflowButtons(1));assert(isfield(fig.UserData.controls,'file'));
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
