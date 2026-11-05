function report=testSU2IdealAirWorkbench(out)
%TESTSU2IDEALAIRWORKBENCH 从主五步骤入口读取两种真实变比热小例并核对导出。
% 模块一。验证运行包、逐变量温度身份、局部步号和界面返回后的输入状态。
% 模块二。使用MATLAB独立物性复核全壁节点数值；错误模型必须拒绝。
assert(~isfolder(out));mkdir(out);
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()close(fig));before=fig.UserData.workflowState;
root=fileparts(fileparts(mfilename('fullpath')));modes={'laminar','sst'};summaries=cell(1,2);
for j=1:2
    folder=fullfile(root,'outputs','v154_validation','cases01',modes{j});
    fig.UserData.actions.readIdealAirReview(folder);drawnow;
    view=fig.UserData.idealAirReview;
    assert(view.summary.interfacePassed && view.summary.nodes==5265 && view.summary.wallNodes==65);
    assert(view.summary.checkpointLocalIteration==9 && ~view.summary.physicsQualified);
    assert(view.convergence.checkpointIdentityConfirmed && ~view.convergence.checkpointResidualPassed);
    assert(strcmp(view.runtime.mode,'ideal-air-variable-cp') && contains(view.method.Text,'OTH'));
    assert(contains(view.method.Text,'Roe') && contains(view.heading.Text,'60–2000'));
    if j==2,assert(contains(view.method.Text,'SST-2003m'));end
    assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.engineJob));
    view.tabs.SelectedTab=view.tableTab;drawnow;
    file=fullfile(out,[modes{j} '.csv']);view.exportTo(file);data=readtable(file,'TextType','string');
    assert(height(data)==65 && all(data.given_Tw_K==300) && all(isnan(data.solvedWallTemperature_K)));
    assert(all(~data.PhysicsQualified) && all(data.SavedLocalIteration==9));
    assert(all(data.Method==string(view.runtime.method)));
    assert(isequal(data.nodeID,view.wall.nodeID));
    assert(max(abs(data.q_W_m2-view.wall.q_W_m2)./abs(view.wall.q_W_m2))<1e-12);
    % 独立MATLAB实现从p/T恢复壁面密度与Cp，避免再次使用定比热公式。
    states=idealAirThermodynamicState('PT',view.wall.p_Pa,view.wall.fluid_T_K);
    assert(all(states.cp>1000));
    receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
    bad=receipt;bad.thermodynamicModel.energyEntropyReference='unknown';rejected=false;
    try,readSU2IdealAirStateAudit(folder,bad);catch,rejected=true;end
    assert(rejected);
    summaries{j}=struct('mode',modes{j},'rows',height(data),'checkpoint',9, ...
        'modelRejected',rejected,'unknownSolvedWallTemperature',true,'summary',view.summary);
end
fid=fopen(fullfile(out,'ready_for_capture.json'),'w');assert(fid>=0);fprintf(fid,'{"ready":true}');fclose(fid);
t=tic;
while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<100,drawnow;pause(.2);end
captured=isfile(fullfile(out,'capture_done.json'));
% 原定比热SST和HEOS回执的描述仍按其各自方法显示。
legacy=jsondecode(fileread(fullfile(root,'outputs','v121_validation','capability01','engine_process.json')));
info=describeSU2Runtime(legacy);assert(strcmp(info.mode,'compressible-heos-air'));
fig.UserData.actions.selectStep(1);drawnow;
assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.idealAirReview));
report=struct('mainWorkbenchEntry',true,'cases',{summaries},'inputStatePreserved',true, ...
    'legacyHEOSDescription',true,'solverRuns',0,'nativeMouseObserved',captured,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);guardFile=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
