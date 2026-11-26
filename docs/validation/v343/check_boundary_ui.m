function check_boundary_ui(out)
% 模块一：只将保存对话框替换为本次验收目录，实际按钮、导出与 WebGL 均使用工作台实现。
dialog=fullfile(out,'dialog');mkdir(dialog);
f=fopen(fullfile(dialog,'uiputfile.m'),'w');assert(f>=0);
fprintf(f,'function [file,folder]=uiputfile(~,~,file)\nfolder=getappdata(0,''rfzCSVFolder'');\nend\n');fclose(f);
setappdata(0,'rfzCSVFolder',out);addpath(dialog,'-begin');
restore=onCleanup(@()restoreDialog(dialog)); %#ok<NASGU>
fig=rocketAeroHeatWorkbench();closeWindow=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.rfzModel);s=fig.UserData.workflowState;
assert(isfield(s.rawMesh,'rfzBoundary'),'%s',fig.UserData.status.Text);
html=fig.UserData.geometryView.html;b=s.rawMesh.rfzBoundary;regions=b.regions;
assert(isequal(s.rawMesh.rfzRegions,regions(:,s.rawMesh.rfzRegions.Properties.VariableNames)));
assert(nnz(regions.referenceRole=="封盖内侧 CAD 壁")==19 && nnz(regions.referenceRole=="实体唇环")==9);
assert(b.condition.publishedCondition.prescribedIsothermalWallK==300 && ~s.boundaryConfirmed);
press(fig.UserData.controls.rfzList);c=fig.UserData.controls;
assert(height(c.rfzTable.Data)==204);
r=regions(105,:);select(c.rfzTable,105);gl=waitFrame(html);
assert(gl.selectedFaceCount==r.lastFace-r.firstFace+1 && gl.faces==130044);
expected=r.firstFace:r.lastFace;
assert(isequal(cell2mat(html.UserData.payload.selectedFaceIDs),expected-1));
probe(html,r.firstFace,html.UserData.revision);
assert(contains(fig.UserData.status.Text,'封盖内侧 CAD 壁'));
press(c.rfzExport);cad=readtable(fullfile(out,'RFZ-REC-deg190-cad.csv'),'TextType','string');
assert(isequal(cad.referenceRole,regions.referenceRole) && isequal(cad.cadFace,regions.cadFace) && ...
    isequal(cad.firstFace,regions.firstFace) && isequal(cad.lastFace,regions.lastFace) && ...
    isequaln(cad.rimCADFace,regions.rimCADFace) && max(abs(cad.stlAreaM2-regions.stlAreaM2))<1e-10);
cadPixels=gl.visiblePixelSamples;

% 模块二：人工封盖行只定位相邻唇环，不向当前网格增加面，也不写入任何物理结果。
c.rfzTableKind.Value='cap';c.rfzTableKind.ValueChangedFcn(c.rfzTableKind,[]);drawnow;
assert(height(c.rfzTable.Data)==9);
j=find(b.caps.rimCADFace==104);select(c.rfzTable,j);gl=waitFrame(html);rim=regions(104,:);
assert(gl.selectedFaceCount==rim.lastFace-rim.firstFace+1 && gl.faces==130044);
probe(html,rim.firstFace,html.UserData.revision);assert(contains(fig.UserData.status.Text,'实体唇环'));
last=fig.UserData.status.Text;probe(html,1,html.UserData.revision-1);assert(strcmp(last,fig.UserData.status.Text));
press(c.rfzExport);caps=readtable(fullfile(out,'RFZ-REC-deg190-cap.csv'));
assert(isequal(caps.rimCADFace,b.caps.rimCADFace) && isequal(caps.referenceInteriorPoints,b.caps.referenceInteriorPoints) && ...
    isequal(caps.farPoints,b.caps.farPoints) && max(abs(caps.polygonAreaM2-b.caps.polygonAreaM2))<1e-12);
capPixels=gl.visiblePixelSamples;
c.rfzTableKind.Value='condition';c.rfzTableKind.ValueChangedFcn(c.rfzTableKind,[]);drawnow;
assert(isempty(c.rfzTable.CellSelectionCallback) && height(c.rfzTable.Data)==10);
press(c.rfzExport);condition=readtable(fullfile(out,'RFZ-REC-deg190-condition.csv'),'TextType','string');
assert(isequal(condition,c.rfzTable.Data) && contains(condition.value(8),'零点待核实') && contains(condition.value(9),'正方向待核实'));
press(c.rfzBack);gl=waitFrame(html);assert(gl.selectedFaceCount==0);
current=fig.UserData.workflowState;
assert(current.revision==s.revision && isempty(current.mesh) && isempty(current.result) && isempty(current.freestream));
c=fig.UserData.controls;c.scale.Value=.001;press(c.load);
assert(~isfield(fig.UserData.controls,'rfzList') && ~isfield(fig.UserData.workflowState.rawMesh,'rfzBoundary'));
record=struct('cadRows',height(cad),'capRows',height(caps),'cadAreaM2',sum(cad.stlAreaM2),'fullReconstructedCapAreaM2',sum(caps.polygonAreaM2), ...
    'farSamples',sum(caps.farPoints),'cadSelectionPixels',cadPixels,'capSelectionPixels',capPixels, ...
    'cavityProbe',105,'rimProbe',104,'conditionRows',height(condition),'expiredProbeIgnored',true,'wrongScaleUnbound',true,'computedFaces',0);
f=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(f>=0);cleanup=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(record,'PrettyPrint',true));
end

function restoreDialog(folder)
rmpath(folder);rmappdata(0,'rfzCSVFolder');
end

function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end

function select(tableControl,index)
tableControl.CellSelectionCallback(tableControl,struct('Indices',[index 1]));drawnow;
end

function probe(html,id,rev)
html.Data=struct('kind','geometry-face-probe','revision',rev,'faceIndex',id-1);
html.DataChangedFcn(html,[]);
end

function gl=waitFrame(html)
t=tic;
while toc(t)<50
    drawnow;s=html.UserData;
    assert(isempty(s.error),'%s',jsonencode(s.error));gl=s.rendered;
    if ~isempty(gl) && gl.revision==s.revision && gl.webglError==0 && gl.visiblePixelSamples>0,return;end
    pause(.1);
end
error('未收到当前选区的绘制回执。');
end
