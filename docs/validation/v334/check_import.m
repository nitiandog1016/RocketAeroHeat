function check_import(out,branch)
% 模块一：通过真实五步侧栏导入标准或换单位、转向、重排的 STL。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(out,'specimen.stl');scale=1;
if strcmp(branch,'half'),scale=.001;end
c.scale.Value=scale;press(c.load);
if strcmp(branch,'half')
    c=fig.UserData.controls;c.rz.Value=-90;press(c.rotate);
end
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'shortSupport65');
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh),'%s',fig.UserData.status.Text);b=s.mesh.shortSupportImport;
order=readmatrix(fullfile(out,'expected_order.txt'));assert(isequal(order,b.referenceFaceID));
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.shortBranch,branch);
press(fig.UserData.controls.confirm);ready=fig.UserData.workflowState;
assert(ready.currentStep==4 && ready.boundaryConfirmed && ~isempty(fig.UserData.controls.run.ButtonPushedFcn));
new=aeroHeatWorkflowState('beginSolve',ready,struct('shortSupportAction','solve'));
assert(strcmp(new.frozenInputs.cfg.shortSupportAction,'solve') && ~isfolder(new.frozenInputs.cfg.shortSupportJob));
if strcmp(branch,'full')
    probe=tempname;mkdir(probe);copyfile(fullfile(fileparts(mfilename('fullpath')),'dispatch_probe.m.in'),fullfile(probe,'runCM011Python.m'));
    addpath(probe,'-begin');cleanPath=onCleanup(@()removeProbe(probe));
    press(fig.UserData.controls.run);stopped=fig.UserData.workflowState;
    assert(strcmp(stopped.error.identifier,'RAH:ShortTest:DispatchObserved') && isempty(stopped.result),'%s',fig.UserData.status.Text);
    clear cleanPath;
end
press(fig.UserData.controls.readShort);s=fig.UserData.workflowState;
assert(s.currentStep==5 && strcmp(s.status,'solved'),'%s',fig.UserData.status.Text);
result=s.result;report=fig.UserData.caseReport;
assert(result.evidenceIdentity.reusedArchive && isempty(result.shortSupportJob) && isequal(result.F,s.orientedMesh.F));
% 模块二：独立取渲染回执，覆盖五区与热峰探针、七站、截线和逐面 CSV。
ui=fig.UserData.webUI;timer=tic;
while isempty(ui.html.UserData.rendered) && toc(timer)<50,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.faces==32480 && gl.visiblePixelSamples>0 && gl.webglError==0);
ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;v=ui.evidenceView;d=v.getData();
assert(isequal(d.shortSupportImport.referenceFaceID,order) && height(report.stations)==7);
assert(isequal(fig.UserData.controls.stations.Data(:,2),report.stations.fileFaceID));
c=fig.UserData.controls;press(c.nativePeak);st=v.getState();assert(st.selectedFaceID==report.nativePeakFileFaceID);
press(c.filePeak);st=v.getState();assert(st.selectedFaceID==report.peakOutFace);
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    field=d.fields(strcmp({d.fields.key},key{1}));assert(all(isnan(field.values)));
end
ids=arrayfun(@(k)find(b.regionID==k,1),1:5);v.selectField('heatFlux');
for i=[ids report.peakOutFace report.stations.fileFaceID']
    v.inspectFace(i,[]);st=v.getState();assert(st.selectedFaceID==i);
end
x=.00137+b.translationM(1);v.setSectionX(x);timer=tic;
while toc(timer)<50
    drawnow;cut=v.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-13,break;end
    pause(.1);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-13);
v.exportSectionTo(fullfile(out,'section.csv'));wall=v.exportWallPairTo(fullfile(out,'wall.csv'));
assert(wall.rows==32480 && wall.joint.pressureHeatFluxQualifiedArea==0);
writematrix(report.reverseFaces,fullfile(out,'reverse.txt'));
proof=struct('branch',branch,'stlSHA256',b.source.sha256,'maxVertexErrorM',b.maxVertexErrorM, ...
    'maxRelativeAreaError',b.maxRelativeAreaError,'translationM',b.translationM,'rotation',b.rotationMatrix, ...
    'peak',report.peakOutFace,'peakReference',report.peakReferenceFaceID,'nativePeakFile',report.nativePeakFileFaceID, ...
    'nativePeakReference',report.nativePeakReferenceFaceID,'stations',report.stations.fileFaceID, ...
    'pixels',gl.visiblePixelSamples,'glError',gl.webglError,'sourceFrame',result.evidenceIdentity.sourceFrame.sha256, ...
    'areaM2',report.jointComputedAreaM2,'reverseFaces',numel(report.reverseFaces),'reverseAreaM2',report.reverseAreaM2, ...
    'elapsedS',s.freestream.elapsedS,'stepS',s.freestream.dtS,'probedFaces',ids);
% 模块三：拒收与冻结请求不符的源，改步长后结果清空；非法输入单独留证。
payload=struct('result',result,'solveId',new.solveId,'sourceRevision',new.solveInputRevision,'freestream',s.freestream);
proof.wrongSource=reject(@()aeroHeatWorkflowState('completeSolve',new,payload),'RAH:Short:Snapshot');
press(fig.UserData.workflowButtons(3));other='half';if strcmp(branch,'half'),other='full';end
choose(fig.UserData.controls.shortBranch,other);changed=fig.UserData.workflowState;
assert(isempty(changed.result) && ~changed.boundaryConfirmed && ~changed.steps(5).available);
if strcmp(branch,'full')
    source=b.source;geo=s.cfg.geometry;mesh=s.orientedMesh;
    geo.scale=.001;changedSource=source;changedSource.scale=.001;bad=mesh;bad.V=bad.V*.001;
    proof.wrongUnits=reject(@()bindShortSupportSTL(bad,changedSource,geo),'RAH:Short:Vertex');
    geo=s.cfg.geometry;geo.rotationMatrix=[0 -1 0;1 0 0;0 0 1];
    bad.V=(mesh.V-geo.rotationCenter)*geo.rotationMatrix'+geo.rotationCenter;
    proof.wrongOrientation=reject(@()bindShortSupportSTL(bad,source,geo),'RAH:Short:Vertex');
    geo=s.cfg.geometry;bad=mesh;bad.F(1,[2 3])=bad.F(1,[3 2]);
    proof.faceChanged=reject(@()bindShortSupportSTL(bad,source,geo),'RAH:Short:FileOrder');
    bad=mesh;bad.F(1,:)=[];
    proof.missingFace=reject(@()bindShortSupportSTL(bad,source,geo),'RAH:Short:Count');
    for name={'reversed','duplicate'}
        kind=name{1};cfg=s.cfg;cfg.geometry.stlFile=fullfile(out,[kind '.stl']);
        damaged=loadSTLMesh(cfg,true);src=source;src.file=cfg.geometry.stlFile;src.sha256=cm011FileSHA256(src.file);
        direct=mesh;
        if strcmp(kind,'duplicate')
            direct.F(1,:)=direct.F(2,:);expected='RAH:Short:Topology';
        else
            direct.F(1,[2 3])=direct.F(1,[3 2]);expected='RAH:Short:Winding';
        end
        % 读取器可能先修正不一致的法向或省去重复面；按文件角点核对仍必须拒绝。
        proof.([kind 'Importer'])=reject(@()bindShortSupportSTL(damaged,src,cfg.geometry),{expected,'RAH:Short:FileOrder','RAH:Short:Count'});
        proof.(kind)=reject(@()bindShortSupportSTL(direct,src,cfg.geometry),expected);
    end
    badCfg=ready.cfg;badCfg.shortSupportBoundary.airK=300;
    invalid=aeroHeatWorkflowState('setBoundary',ready,struct('cfg',badCfg));
    proof.alteredBoundary=reject(@()aeroHeatWorkflowState('confirmBoundary',invalid),'RAH:Short:Boundary');
    proof.runDispatch=stopped.error.identifier;
else
    % 归档导出在本机保留一份，复核所有面上的热流、温度缺项与近壁诊断。
    exportAerothermalWeb(result,fullfile(out,'workbench'));
end
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);
close=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true)); %#ok<NASGU>
end
function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
function choose(control,value)
control.Value=value;control.ValueChangedFcn(control,[]);drawnow;
end
function id=reject(action,expected)
try,action();catch ex,id=ex.identifier;assert(any(strcmp(id,expected)),'%s %s',id,ex.message);return;end
error('负例意外通过。');
end
function removeProbe(folder)
% 测试替身在系统临时目录，用后清除，日常项目路径始终使用真实原生入口。
rmpath(folder);delete(fullfile(folder,'runCM011Python.m'));rmdir(folder);
end
