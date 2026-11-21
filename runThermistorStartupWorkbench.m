function [mesh,report]=runThermistorStartupWorkbench(input,cfg,progress)
% 模块一：核对确认的原面、全周体网和边界，每次显式新算均使用独立目录。
root=fileparts(mfilename('fullpath'));base=fullfile(root,'outputs','v294_validation');
verifyThermistorWorkbenchMesh(input,48);
thermistorWorkbenchBoundary(cfg.thermistorBoundary,48);
action='run';if isfield(cfg,'thermistorResultSource'),action=cfg.thermistorResultSource;end
switch action
    case {'solve','run'}
        action='run';
        if isfield(cfg,'thermistorRunDirectory'),folder=cfg.thermistorRunDirectory;
        else,folder=fullfile(base,'workbench',char(java.util.UUID.randomUUID()));end
    case {'read','resume'}
        assert(isfield(cfg,'thermistorRunDirectory') && isfolder(cfg.thermistorRunDirectory), ...
            'RAH:Bead:Directory','请选择已有 6 m/s 计算的根目录。');
        folder=char(java.io.File(cfg.thermistorRunDirectory).getCanonicalPath());
    otherwise
        error('RAH:Bead:Mode','6 m/s 工况须选择新算、读取或接续。');
end
work=fullfile(root,'outputs',['thermistor_view_' char(java.util.UUID.randomUUID())]);
mkdir(work);receipt=fullfile(work,'result.json');logFile=fullfile(work,'process.log');
script=fullfile(root,'scripts','runThermistorStartupJob.py');
stateFile=fullfile(folder,'progress.json');last='';started=tic;tick=tic;
notify(['6 m/s 热敏珠计算目录 ' folder],'started');
runCM011Python({script,action,folder,'--receipt',receipt},logFile,@heartbeat,true);

% 模块二：只装入经两半步实际终场与守恒审核确认的逐面数据，绑定到本次工作流。
job=jsondecode(fileread(receipt));
assert(strcmpi(job.directory,folder),'RAH:Bead:Source','计算回执与所选目录不符。');
[mesh,report]=readThermistorStartup('nearWall6',job);
assert(isequal(input.V,mesh.V) && isequal(input.F,mesh.F),'RAH:Bead:Geometry','本次结果的原面不符。');
p=input.thermistorPrepared;
mesh.evidenceIdentity.meshSHA256=p.meshSHA256;mesh.evidenceIdentity.preparedSHA256=p.manifestSHA256;
mesh.evidenceIdentity.runDirectory=folder;mesh.evidenceIdentity.reusedArchive=false;
mesh.evidenceIdentity.startupAction=action;report.outputDirectory=folder;
report.startupAction=action;report.resultReceipt=receipt;
report.workbenchPage=exportAerothermalWeb(mesh,fullfile(work,'workbench'));

    function heartbeat()
        % 仅在阶段改变时更新日志，持续处理窗口关闭和界面事件。
        drawnow limitrate;
        if toc(tick)<1,return;end
        tick=tic;stage='';
        if isfile(stateFile)
            record=jsondecode(fileread(stateFile));stage=char(record.stage);
        end
        if ~isempty(stage) && ~strcmp(stage,last)
            names={'prepare_inputs','first01','first_state','first_observe','first_release','first_energy','second01', ...
                'second_state','second_observe','second_release','second_energy','compare'};
            labels={'核对或重建启动初场与观察器','第一半步推进','第一半步全域状态','第一半步固定端点回放','第一半步独立壁反力','第一半步热量收支', ...
                '第二半步推进','第二半步全域状态','第二半步固定端点回放','第二半步独立壁反力','第二半步热量收支','同面近壁对照'};
            k=find(strcmp(names,stage),1);label=stage;if ~isempty(k),label=labels{k};end
            notify(label,'started');last=stage;
        else
            notify('','running');
        end
    end
    function notify(message,status)
        if isempty(progress),return;end
        progress(struct('phase','thermistorStartup','status',status,'fraction',NaN, ...
            'message',message,'elapsedSeconds',toc(started),'errorIdentifier','','errorMessage',''));
    end
end
