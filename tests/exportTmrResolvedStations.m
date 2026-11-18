function report=exportTmrResolvedStations(folder)
%EXPORTTMRRESOLVEDSTATIONS 将已审核的 NASA 风洞测站网格终态接入离线工作台。
% 模块一：从通过独立审核的文件身份中读回原生初场和终场。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'), ...
    fullfile(root,'internalFlow','cutCell'));
runFolder=fullfile(root,'outputs','v263_validation','tmr_resolved01');
archive=fullfile(root,'docs','validation','v263','evidence','tmr_resolved');
receipt=jsondecode(fileread(fullfile(archive,'audit_process.json')));
summary=jsondecode(fileread(fullfile(archive,'summary.json')));
assert(receipt.exitCode==0 && receipt.inputsStable && ~summary.physicalQualification);
digests=struct2cell(receipt.inputSHA256);
assert(~isfolder(folder),'输出目录已经存在。');
start=readCheckpoint('initial.mat');
terminal=readCheckpoint('terminal.mat');
r=terminal.implicit;mesh=terminal.mesh;volume=mesh.volumePerRadian;
n=terminal.report.acceptedSteps;
assert(isequaln(start.mesh,mesh) && isequaln(start.gas,terminal.gas) && isequaln(start.boundary,terminal.boundary));
assert(isequal(r.residualScale,start.scale) && n==summary.chain.acceptedEndpoints && ...
    numel(dir(fullfile(runFolder,'step_*.mat')))==n && terminal.report.physicalTimeAdvanced==0);

% 模块二：按原始保存步拼接历史，累计守恒用整段初场作为基准。
% 各步末行的 Newton 缺陷、目标和单步收支来自该步原记录。
residuals=r.history.residual(:);
h=struct('residual',residuals,'stepDefect',zeros(n+1,5),'balance',zeros(n+1,1), ...
    'newtonDefect',zeros(n+1,1),'newtonTarget',zeros(n+1,1));
assert(numel(residuals)==n+1 && residuals(1)==max(abs(start.rate)./start.scale,[],'all'));
dt=zeros(n,1);
for j=1:n
    step=readCheckpoint(sprintf('step_%03d.mat',j));
    s=step.history;
    assert(step.steps==j && isequal(s.residual(:),residuals(j:j+1)) && ...
        s.dt(end)==step.dt && step.dt>0 && isequal(size(s.stepDefect),[2 5]));
    if j==1
        h.stepDefect(1,:)=s.stepDefect(1,:);
        h.newtonDefect(1)=s.newtonDefect(1);h.newtonTarget(1)=s.newtonTarget(1);
        h.balance(1)=s.balance(1);
    end
    h.stepDefect(j+1,:)=s.stepDefect(end,:);
    h.newtonDefect(j+1)=s.newtonDefect(end);h.newtonTarget(j+1)=s.newtonTarget(end);
    h.balance(j+1)=step.balance;
    dt(j)=step.dt;
end
assert(isequal(r.Q,step.Q) && step.balance==r.globalBalanceRelative);
r.history=h;r.initialQ=start.initialQ;r.initialConservedIntegral=start.initial;
r.finalConservedIntegral=sum(r.Q.*volume,1);
r.accumulatedBoundaryIntegral=step.net;r.accumulatedSourceIntegral=step.src;
r.acceptedSteps=n;r.stopReason=terminal.report.stopReason;
% 原驱动每次只执行一个接受步，初始 CFL 也会更新。这里只传递整段共用的验收参数。
criteria={'steadyResidualTolerance','steadyRequiredSteps','minimumSteadySteps', ...
    'conservationTolerance','nonlinearAbsoluteTolerance','nonlinearRelativeTolerance'};
r.options=rmfield(start.opts,setdiff(fieldnames(start.opts),criteria));
r.geometryBinding=struct('X',start.mesh.X,'R',start.mesh.R,'cells',start.mesh.cells,'volumePerRadian',volume);
r.pseudoTime=sum(dt);r.physicalTimeAdvanced=terminal.report.physicalTimeAdvanced;r.pseudoTimeIsPhysical=false;
r.gas=terminal.gas;r.wallTemperatureValueKind='prescribed';

% 模块三：重用已有工作台导出入口，逐值复核网格、完整空间算子、五量收支和壁面量。
source=['NASA ACC-SBLI, tmr_resolved01/terminal.mat, SHA256 ' summary.terminalMAT_SHA256];
opts=struct('name',sprintf('NASA 风洞 · %d 个壁面环带 · 双量均在实验估计内 %d/%d 站', ...
    summary.wallRings,numel(summary.insideAtJointStationsCm),summary.jointStations), ...
    'title','RocketAeroHeat · NASA 风洞壁面测站','evidenceSource',source);
report=exportTmrAxisymmetricSAStudy(mesh,r,terminal.boundary,folder,opts);

    function frame=readCheckpoint(name)
        filename=fullfile(runFolder,name);
        digest=rocketStreamingFileSHA256(filename);
        assert(any(strcmp(digests,digest)),'文件 SHA 与通过审核的保存链不符：%s',name);
        if strcmp(name,'terminal.mat')
            assert(strcmp(digest,summary.terminalMAT_SHA256));
        end
        frame=load(filename);
    end
end
