function check_air_workbench(out,ticketFile)
% 模块一：从工作台打开标准工况及已有的旋转、分温工况，验证显式绑定和撤销。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dialog=fullfile(out,'dialog');mkdir(dialog);
f=fopen(fullfile(dialog,'uigetfile.m'),'w');
fprintf(f,'function [f,p]=uigetfile(varargin)\n[p,n,e]=fileparts(getappdata(0,''airTestInput''));f=[n e];\nend\n');fclose(f);
addpath(dialog,'-begin');restore=onCleanup(@()removeDialog(dialog)); %#ok<NASGU>
files={fullfile(root,'cases','rfz_rec','full_cad_input.json'), ...
    fullfile(root,'docs','validation','v344','evidence','rotated_input.json')};
ticket=jsondecode(fileread(ticketFile));
d=load(fullfile(ticket.directory,'air.mat'),'wallOwner','cellSICN');
probeIDs=[1 71861 130044];evidence=struct([]);
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
for k=1:2
    setappdata(0,'airTestInput',files{k});
    if k==1,press(fig.UserData.controls.browse);else,press(fig.UserData.controls.rfzReopen);end
    state=fig.UserData.workflowState;assert(state.currentStep==4 && state.boundaryConfirmed);
    html=fig.UserData.geometryView.html;
    clock=tic;
    while toc(clock)<55
        drawnow;gl=html.UserData.rendered;
        if ~isempty(gl) && gl.revision==html.UserData.revision && gl.visiblePixelSamples>0,break;end
        pause(.1);
    end
    assert(~isempty(gl) && gl.revision==html.UserData.revision && gl.webglError==0 && gl.faces==130044);
    before=probe(fig,html,probeIDs(2));
    assert(~contains(before,'空气单元') && ~contains(fig.UserData.controls.rfzAirInfo.Text,'个四面体'));
    press(fig.UserData.controls.rfzAir);
    assert(contains(fig.UserData.controls.rfzAirInfo.Text,'185078 个四面体'),'%s',fig.UserData.status.Text);
    labels=cell(size(probeIDs));
    for j=1:numel(probeIDs)
        id=probeIDs(j);owner=d.wallOwner(id);labels{j}=probe(fig,html,id);
        assert(contains(labels{j},sprintf('空气单元 %d · SICN %.5g',owner,d.cellSICN(owner))));
        assert(contains(labels{j},sprintf('给定 %.6g K',state.cfg.heating.wallTemperatureByFace(id))));
    end
    after=fig.UserData.workflowState;
    assert(isequaln(state,after) && isempty(after.result) && ~after.steps(5).available);
    evidence(k).meshSHA256=state.cfg.rfzSurfaceInput.meshSHA256;
    evidence(k).wallSHA256=state.cfg.rfzSurfaceInput.wallSHA256;
    evidence(k).cellIDs=d.wallOwner(probeIDs);
    evidence(k).probeLabels=labels;
    evidence(k).pixels=gl.visiblePixelSamples;
end

% 模块二：修改已确认的边界后，原有探针失效；重确认后仍需要显式绑定。
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;c.innerK.Value=316;
c.innerK.ValueChangedFcn(c.innerK,[]);drawnow;
assert(~fig.UserData.workflowState.boundaryConfirmed);
press(c.confirm);assert(fig.UserData.workflowState.currentStep==4);
label=probe(fig,html,probeIDs(2));assert(contains(label,'316 K') && ~contains(label,'空气单元'));
rec=struct('cases',evidence,'faceIDs',probeIDs,'afterBoundaryChange',label, ...
    'computedFaces',0,'computedAreaM2',0);
f=fopen(fullfile(out,'workbench.json'),'w','n','UTF-8');
fprintf(f,'%s\n',jsonencode(rec,'PrettyPrint',true));fclose(f);
end

function label=probe(fig,html,id)
html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',id-1);
html.DataChangedFcn(html,[]);label=fig.UserData.status.Text;
end
function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end
function removeDialog(folder)
rmpath(folder);rmappdata(0,'airTestInput');
end
