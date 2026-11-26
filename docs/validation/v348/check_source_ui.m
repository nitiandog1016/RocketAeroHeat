function check_source_ui(out)
% 模块一：通过工作台打开参考数据，核对显示、两种导出格式以及返回后的工况。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dialog=fullfile(out,'dialog');mkdir(dialog);
writeText(fullfile(dialog,'uigetfile.m'),"function [f,p]=uigetfile(varargin),[p,n,e]=fileparts(getappdata(0,'refSourceInput'));f=[n e];end");
writeText(fullfile(dialog,'uiputfile.m'),"function [f,p]=uiputfile(varargin),[p,n,e]=fileparts(getappdata(0,'refSourceExport'));f=[n e];end");
setappdata(0,'refSourceInput',fullfile(root,'cases','rfz_rec','full_cad_input.json'));
setappdata(0,'refSourceExport','');addpath(dialog,'-begin');
cleanup=onCleanup(@()restore(dialog)); %#ok<NASGU>
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.browse);s=fig.UserData.workflowState;
assert(s.boundaryConfirmed && s.currentStep==4 && isempty(s.result));
press(fig.UserData.controls.rfzSamples);c=fig.UserData.controls;h=fig.UserData.geometryView.html;
frame(h);caption=h.UserData.payload.caption;note=c.sampleSource.Text;
assert(contains(caption,'5.35（发布说明）') && contains(caption,'5.38（工况表）'));
assert(contains(note,'表 4 编号 3') && contains(note,'Case 2（第 4 行）'));
assert(strcmp(h.UserData.payload.pointField.convention,'作者热流正方向待核实'));
c.sampleField.Value='p';c.sampleField.ValueChangedFcn(c.sampleField,[]);frame(h);
assert(strcmp(h.UserData.payload.pointField.convention,'作者压力零点待核实'));
png=fullfile(out,'pressure.png');setappdata(0,'refSourceExport',png);press(c.samplePNG);
t=tic;while ~isfile(png) && toc(t)<30,drawnow;pause(.05);end
assert(isfile(png),'%s',fig.UserData.status.Text);imfinfo(png);
bundle=fullfile(out,'reference.zip');setappdata(0,'refSourceExport',bundle);press(c.sampleExport);
assert(isfile(bundle),'%s',fig.UserData.status.Text);
csv=fullfile(out,'direct.csv');setappdata(0,'refSourceExport',csv);press(c.sampleExport);assert(isfile(csv));
press(c.samplesBack);frame(h);assert(isequaln(s,fig.UserData.workflowState));

% 模块二：给错误字段附来源时应拒绝，已有压缩包仍可使用。
before=cm011FileSHA256(bundle);
bad=readRFZReferenceSamples();bad.authorPressure_Pa(1)=bad.authorPressure_Pa(1)+1;
rejected=false;try,saveRFZReferenceBundle(bundle,bad);catch ex,rejected=strcmp(ex.identifier,'RAH:RFZ:ReferenceSource');end
assert(rejected && strcmp(before,cm011FileSHA256(bundle)));
writeText(fullfile(out,'ui.json'),jsonencode(struct('caption',caption,'sourceNote',note, ...
    'csvSHA256',cm011FileSHA256(csv),'bundleSHA256',before,'pngSHA256',cm011FileSHA256(png), ...
    'alteredDataRejected',rejected,'projectComputedFaces',0,'projectComputedAreaM2',0),'PrettyPrint',true));
end

function frame(h)
t=tic;
while toc(t)<45
    drawnow;s=h.UserData;assert(isempty(s.error),'%s',jsonencode(s.error));
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.webglError==0 && s.rendered.visiblePixelSamples>0,return;end
    pause(.05)
end
error('等待当前绘图回执超时。');
end
function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end
function writeText(path,value)
f=fopen(path,'w','n','UTF-8');assert(f>=0);c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',value);
end
function restore(dialog)
rmpath(dialog);rmappdata(0,'refSourceInput');rmappdata(0,'refSourceExport');
end
