function testStaticAirProbeWorkbench(out)
%TESTSTATICAIRPROBEWORKBENCH 从主窗口按钮进入真实网页，并核对切换回执及 CSV。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(~isfolder(out));mkdir(out);
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig));
before=fig.UserData.workflowState;
button=fig.UserData.controls.staticAir;button.ButtonPushedFcn(button,[]);
html=findall(fig.UserData.resultHost,'Tag','static-air-probe');
assert(isscalar(html),'StaticAir:Host','节点页未创建：%s',fig.UserData.status.Text);
assert(strcmp(fig.UserData.resultHost.Visible,'on') && isempty(fig.UserData.engineJob));
writeText(fullfile(out,'index.html'),html.HTMLSource);
% 模块一：真实 uihtml 的选择确认，保留内部节点、边界输入及对应源帧。
requests={struct,struct('node',62),struct('sample','stable'),struct('frame','terminal'),struct('node',124)};
cases={'default','default','stable','stable','stable'};frames={'initial','initial','initial','terminal','terminal'};
ids=[0 62 62 62 124];observed=cell(1,numel(ids));
for k=1:numel(ids)
    if k>1,command=requests{k};command.kind='static-air-select';html.Data=command;end
    t=tic;found=false;
    while toc(t)<20
        drawnow;event=html.Data;
        found=isstruct(event) && isfield(event,'kind') && strcmp(event.kind,'static-air-selection') && ...
            strcmp(event.sample,cases{k}) && strcmp(event.frame,frames{k}) && event.node==ids(k);
        if found,break;end
        pause(.05);
    end
    assert(found,'StaticAir:Selection','节点页没有确认 %s/%s 的原节点 %d。',cases{k},frames{k},ids(k));
    file=fullfile(out,event.filename);writeText(file,event.csv);
    observed{k}=event.filename;
end
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.engineJob));
fig.UserData.actions.selectStep(1);drawnow;
assert(isempty(findall(fig.UserData.resultHost,'Tag','static-air-probe-host')));
assert(strcmp(fig.UserData.resultHost.Visible,'off') && isequaln(before,fig.UserData.workflowState));
% 模块二：复用实际清单，配置字节更改时生成器须拒绝，原归档保持只读。
evidence=fullfile(root,'docs','validation','v261','evidence');
bad=fullfile(out,'changed-source');mkdir(bad);mkdir(fullfile(bad,'default'));
copyfile(fullfile(evidence,'manifest.json'),fullfile(bad,'manifest.json'));
writeText(fullfile(bad,'default','case.cfg'),'changed input');
reason='';
try,buildStaticAirProbePage(bad);catch ex,reason=ex.identifier;end
assert(strcmp(reason,'StaticAir:SourceChanged'),'篡改配置未按来源校验拒绝。');
report=struct('mainButton',true,'embeddedSelections',{observed},'returnedToWorkflow',true, ...
    'sourceRejection',reason,'solverRuns',0,'addedQualifiedAreaM2',0);
writeText(fullfile(out,'report.json'),jsonencode(report));
end

function writeText(file,text)
fid=fopen(file,'wb');assert(fid>=0);guard=onCleanup(@()fclose(fid));
bytes=unicode2native(text,'UTF-8');assert(fwrite(fid,bytes,'uint8')==numel(bytes));
end
