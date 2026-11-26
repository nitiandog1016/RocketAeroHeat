function check_reference_ui(out)
% 模块一：打开完整工况、绘制作者点云，检查原值、色标和输入状态隔离。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dialog=fullfile(out,'dialog');mkdir(dialog);
writeText(fullfile(dialog,'uigetfile.m'),"function [f,p]=uigetfile(varargin),[p,n,e]=fileparts(getappdata(0,'sampleCase'));f=[n e];end");
writeText(fullfile(dialog,'uiputfile.m'),"function [f,p]=uiputfile(varargin),[p,n,e]=fileparts(getappdata(0,'sampleExport'));f=[n e];end");
setappdata(0,'sampleCase',fullfile(root,'cases','rfz_rec','full_cad_input.json'));
export=fullfile(out,'samples.csv');setappdata(0,'sampleExport',export);
addpath(dialog,'-begin');finish=onCleanup(@()restore(dialog)); %#ok<NASGU>
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.browse);state=fig.UserData.workflowState;
assert(state.currentStep==4 && state.boundaryConfirmed && isempty(state.result));
html=fig.UserData.geometryView.html;wall=frame(html);assert(wall.faces==130044);
press(fig.UserData.controls.rfzSamples);c=fig.UserData.controls;
assert(isfield(c,'sampleProbe'),'%s',fig.UserData.status.Text);
qframe=frame(html);q=html.UserData.payload.pointField.values;
assert(qframe.sampleCount==72706 && qframe.faces==0 && contains(qframe.sampleLegend,'W/m²'));
assert(qframe.fieldRange(1)>min(q) && qframe.fieldRange(2)<max(q) && nnz(q<0)==2944);
pngFile=fullfile(out,'reference.png');setappdata(0,'sampleExport',pngFile);press(c.samplePNG);
t=tic;while ~isfile(pngFile) && toc(t)<25,drawnow;pause(.05);end
assert(isfile(pngFile),'%s',fig.UserData.status.Text);png=imfinfo(pngFile);
assert(strcmpi(png.Format,'png') && png.Width>=760 && png.Height>150);setappdata(0,'sampleExport',export);
c.sampleScale.Value='full';c.sampleScale.ValueChangedFcn(c.sampleScale,[]);qframe=frame(html);
assert(isequal(qframe.fieldRange(:),[min(q);max(q)]));
xyz=html.UserData.payload.vertices;[~,id]=min(q);
c.sampleID.Value=id;c.sampleID.ValueChangedFcn(c.sampleID,[]);highlight=waitSelection(html,id);
assert(contains(c.sampleProbe.Text,sprintf('q = %.10g W/m²',q(id))));
pick(html,72706);assert(c.sampleID.Value==72706);
stale=html.UserData.revision;
c.sampleScale.Value='central';c.sampleScale.ValueChangedFcn(c.sampleScale,[]);trim=frame(html);
assert(trim.fieldRange(1)>qframe.fieldRange(1) && trim.fieldRange(2)<qframe.fieldRange(2));
assert(isequal(html.UserData.payload.pointField.values,q));
c.sampleScale.Value='fixed';c.sampleScale.ValueChangedFcn(c.sampleScale,[]);frame(html);
c.sampleLo.Value=-1000;c.sampleHi.Value=25000;press(c.sampleApply);fixed=frame(html);
assert(isequal(fixed.fieldRange(:),[-1000;25000]));
c.sampleLo.Value=26000;before=html.UserData.revision;press(c.sampleApply);
assert(html.UserData.revision==before && contains(fig.UserData.status.Text,'上限须大于'));
c.sampleField.Value='p';c.sampleField.ValueChangedFcn(c.sampleField,[]);pframe=frame(html);
p=html.UserData.payload.pointField.values;assert(isequal(size(q),size(p)) && ~isequal(q,p));
assert(isequal(pframe.fieldRange(:),[min(p);max(p)]) && contains(pframe.sampleLegend,'Pa'));
html.Data=struct('kind','geometry-point-probe','revision',stale,'sampleIndex',0);
html.DataChangedFcn(html,[]);assert(c.sampleID.Value==72706);
c.sampleField.Value='q';c.sampleField.ValueChangedFcn(c.sampleField,[]);again=frame(html);
assert(isequal(again.fieldRange(:),[-1000;25000]) && c.sampleID.Value==72706);
press(c.sampleExport);assert(isfile(export));
press(c.samplesBack);back=frame(html);assert(back.faces==130044 && ~isfield(back,'sampleCount'));
assert(isequaln(state,fig.UserData.workflowState));

% 模块二：从另一给温和姿态重开，作者的参考坐标、范围和温度仍属固定发布物。
setappdata(0,'sampleCase',fullfile(root,'docs','validation','v344','evidence','rotated_input.json'));
press(fig.UserData.controls.rfzReopen);changed=fig.UserData.workflowState;
assert(changed.cfg.rfzSurfaceInput.cavityWallK==315);
press(fig.UserData.controls.rfzSamples);frame(html);
assert(isequal(xyz,html.UserData.payload.vertices) && isequal(q,html.UserData.payload.pointField.values));
c=fig.UserData.controls;assert(contains(c.sampleProbe.Text,'300 K'));
press(c.samplesBack);split=frame(html);assert(isequal(split.fieldRange(:),[300;315]));
assert(isequaln(changed,fig.UserData.workflowState));
receipt=struct('samples',numel(q),'renderedTriangles',qframe.faces,'pixels',qframe.visiblePixelSamples, ...
    'negativeHeatFluxSamples',nnz(q<0),'qRange',qframe.fieldRange,'pRange',pframe.fieldRange, ...
    'displayCentralRange',trim.fieldRange,'fixedRange',fixed.fieldRange,'minimumSampleID',id, ...
    'highlightedSampleID',highlight.selectedSampleIndex+1,'exportSHA256',cm011FileSHA256(export), ...
    'projectComputedFaces',0,'projectComputedAreaM2',0,'sourceCoordinatesSurvivePoseChange',true);
writeText(fullfile(out,'ui.json'),jsonencode(receipt,'PrettyPrint',true));
end

function v=waitSelection(html,id)
t=tic;
while toc(t)<30
    v=frame(html);if v.selectedSampleIndex==id-1,return;end
    pause(.05);
end
error('样点标记未更新。');
end
function pick(html,id)
html.Data=struct('kind','geometry-point-probe','revision',html.UserData.revision,'sampleIndex',id-1);
html.DataChangedFcn(html,[]);drawnow;
end
function v=frame(html)
t=tic;
while toc(t)<45
    drawnow;s=html.UserData;assert(isempty(s.error),'%s',jsonencode(s.error));v=s.rendered;
    if ~isempty(v) && v.revision==s.revision && v.webglError==0 && v.visiblePixelSamples>0,return;end
    pause(.05)
end
error('未收到当前参考画面的回执。');
end
function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end
function writeText(path,value)
f=fopen(path,'w','n','UTF-8');assert(f>=0);c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',value);
end
function restore(dialog)
rmpath(dialog);rmappdata(0,'sampleCase');rmappdata(0,'sampleExport');
end
