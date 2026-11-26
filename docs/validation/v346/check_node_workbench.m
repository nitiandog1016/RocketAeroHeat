function check_node_workbench(out)
% 模块一：模拟文件选择，检查真实工作台的原面温界、选区与导出，不启动流场计算。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
dialog=fullfile(out,'dialog');mkdir(dialog);
writeText(fullfile(dialog,'uigetfile.m'),"function [f,p]=uigetfile(varargin),[p,n,e]=fileparts(getappdata(0,'nodeInput'));f=[n e];end");
writeText(fullfile(dialog,'uiputfile.m'),"function [f,p]=uiputfile(varargin),[p,n,e]=fileparts(getappdata(0,'nodeOutput'));f=[n e];end");
addpath(dialog,'-begin');clean=onCleanup(@()removeDialogs(dialog)); %#ok<NASGU>
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
inputs={fullfile(root,'cases','rfz_rec','full_cad_input.json'), ...
    fullfile(root,'docs','validation','v344','evidence','rotated_input.json')};
names={'uniform','split'};receipt=struct([]);
for k=1:2
    setappdata(0,'nodeInput',inputs{k});c=fig.UserData.controls;
    if k==1,press(c.browse);else,press(c.rfzReopen);end
    c=fig.UserData.controls;assert(~isfield(c,'rfzNodes'));press(c.rfzAir);
    c=fig.UserData.controls;assert(isfield(c,'rfzNodes'),'%s',fig.UserData.status.Text);
    press(c.rfzNodes);c=fig.UserData.controls;
    t=c.nodeTable;assert(c.nodeScope.Value==13 && sum(t.Data.faces)==22852);
    same=find(t.Data.constraint=="四壁点同温");assert(t.Data.faces(same)==2409);
    pick(t,same);html=fig.UserData.geometryView.html;gl=waitFrame(html);
    ids=cell2mat(html.UserData.payload.selectedFaceIDs)+1;
    assert(numel(ids)==2409 && gl.selectedFaceCount==2409 && all(ids>=4727 & ids<=27578));
    probe(html,ids(1));assert(contains(c.nodeProbe.Text,'四壁点同温') && contains(c.nodeProbe.Text,'实体壁'));
    receipt(k).sample=c.nodeProbe.Text;receipt(k).cad13Pixels=gl.visiblePixelSamples;
    target=fullfile(out,[names{k} '-cad13.csv']);setappdata(0,'nodeOutput',target);press(c.nodeFacesCSV);
    region=readtable(target,'TextType','string');
    assert(isequal(region.faceID,(4727:27578)') && all(region.cadFace==13));
    assert(isequal(ids(:),region.faceID(region.nodalConstraint=="四壁点同温")));
    scope=c.nodeScope;scope.Value=0;scope.ValueChangedFcn(scope,[]);gl=waitFrame(html);
    assert(gl.selectedFaceCount==0 && sum(t.Data.faces)==130044);
    expired=html.UserData.revision;pick(t,find(t.Data.constraint=="节点温界冲突"));gl=waitFrame(html);
    assert(gl.selectedFaceCount==1170*(k==2));
    if k==2
        mixed=cell2mat(html.UserData.payload.selectedFaceIDs)+1;probe(html,mixed(1));
        assert(contains(c.nodeProbe.Text,'节点温界冲突'));receipt(k).mixedSample=c.nodeProbe.Text;
        old=c.nodeProbe.Text;
        html.Data=struct('kind','geometry-face-probe','revision',expired,'faceIndex',0);
        html.DataChangedFcn(html,[]);assert(strcmp(c.nodeProbe.Text,old));
    end
    setappdata(0,'nodeOutput',fullfile(out,[names{k} '-faces.csv']));press(c.nodeFacesCSV);
    setappdata(0,'nodeOutput',fullfile(out,[names{k} '-nodes.csv']));press(c.wallNodesCSV);
    receipt(k).summary=c.nodeSummary.Text;receipt(k).computedFaces=0;
    s=fig.UserData.workflowState;assert(s.boundaryConfirmed && isempty(s.result) && ~s.steps(5).available);
    press(c.nodesBack);gl=waitFrame(html);assert(gl.selectedFaceCount==0 && ~isempty(gl.fieldRange));
end
% 模块二：修改边界后旧绑定、选区及明细不再可用。
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;c.innerK.Value=316;
c.innerK.ValueChangedFcn(c.innerK,[]);press(c.confirm);
assert(~isfield(fig.UserData.controls,'rfzNodes'));
writeText(fullfile(out,'ui.json'),jsonencode(receipt,'PrettyPrint',true));
end

function pick(t,i)
t.CellSelectionCallback(t,struct('Indices',[i 1]));drawnow;
end
function probe(html,id)
html.Data=struct('kind','geometry-face-probe','revision',html.UserData.revision,'faceIndex',id-1);
html.DataChangedFcn(html,[]);drawnow;
end
function gl=waitFrame(html)
clock=tic;
while toc(clock)<55
    drawnow;s=html.UserData;assert(isempty(s.error),'%s',jsonencode(s.error));gl=s.rendered;
    if ~isempty(gl) && gl.revision==s.revision && gl.webglError==0 && gl.visiblePixelSamples>0,return;end
    pause(.1);
end
error('未收到当前选区的绘制回执。');
end
function press(b)
b.ButtonPushedFcn(b,[]);drawnow;
end
function writeText(path,value)
f=fopen(path,'w','n','UTF-8');assert(f>=0);clean=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',value);
end
function removeDialogs(path)
rmpath(path);rmappdata(0,'nodeInput');if isappdata(0,'nodeOutput'),rmappdata(0,'nodeOutput');end
end
