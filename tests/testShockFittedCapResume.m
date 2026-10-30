function report=testShockFittedCapResume(folder,sourceFile)
%TESTSHOCKFITTEDCAPRESUME 检查检查点续算逐位相同、输入不符拒绝和区域数组预算。
% 模块一：同一小球鼻六步整算，与三步保存后三步续算直接比较全部五量/节点/历史。
% 模块二：改上游工况须拒绝；数组预算不足不能启动实际推进。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupShockFittedFlow;
assert(~exist(folder,'dir'),'ShockFit:Overwrite','续算测试使用新目录。');mkdir(folder);
source=load(sourceFile,'fs');fs=source.fs;
up=struct('pressure',double(fs.pressure),'temperature',double(fs.temperature),'density',double(fs.density), ...
    'velocityVector',[-norm(double(fs.flowVelocityVector));0;0],'method','零攻角球鼻续算一致性检查');
gas=createSurfaceFlowGas('caloricallyPerfect',double(fs.R),1.4);
cap=buildShockFittedSphericalCap(.05,60,2,6,.008,2);
allOptions=struct('maximumSteps',6,'endTimeS',1e-4);complete=solveShockFittedCap3D(up,gas,cap,allOptions);
partOptions=allOptions;partOptions.maximumSteps=3;partOptions.checkpointFile=fullfile(folder,'first_checkpoint.mat');
part=solveShockFittedCap3D(up,gas,cap,partOptions);
assert(exist(partOptions.checkpointFile,'file')==2,'ShockFit:Checkpoint','检查点没有实际落盘。');
continueOptions=partOptions;continueOptions.resumeFile=partOptions.checkpointFile;
continueOptions.checkpointFile=fullfile(folder,'resumed_checkpoint.mat');resumed=solveShockFittedCap3D(up,gas,cap,continueOptions);
assert(isequal(complete.Q,resumed.Q)&&isequal(complete.nodes,resumed.nodes)&& ...
    isequal(complete.history,resumed.history)&&isequal(complete.conservedIntegrals,resumed.conservedIntegrals)&& ...
    isequal(complete.boundaryFluxRates,resumed.boundaryFluxRates),'ShockFit:ResumeReplay','分段续算没有逐位复现同一推进序列。');
assert(~exist([partOptions.checkpointFile '.writing.mat'],'file')&& ...
    ~exist([continueOptions.checkpointFile '.writing.mat'],'file'),'ShockFit:Checkpoint','完整提交后仍残留临时检查点。');
changed=up;changed.velocityVector(1)=changed.velocityVector(1)-1;
refusals=struct;
refusals.inputMismatch=reject(@()solveShockFittedCap3D(changed,gas,cap,continueOptions),'ShockFit:Resume');
budget=allOptions;budget.maximumArrayBytes=1000;
refusals.regionBudget=reject(@()solveShockFittedCap3D(up,gas,cap,budget),'ShockFit:RegionMemoryBudget');
report=struct('status','tested-bitwise-curved-region-checkpoint-resume-and-refusals','steps',6, ...
    'bitwiseQ',true,'bitwiseNodes',true,'bitwiseHistory',true,'bitwiseConservationLogs',true, ...
    'checkpointActuallyWritten',true,'productionQualified',false,'engineProcessID',feature('getpid'));
save(fullfile(folder,'shock_resume_evidence.mat'),'report','complete','part','resumed','refusals','-v7');
fid=fopen(fullfile(folder,'shock_resume_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));disp(report);
end
function id=reject(action,expected)
id='';try,action();catch ex,id=ex.identifier;end
assert(strcmp(id,expected),'ShockFit:RefusalTest','拒绝身份不符：实际%s，需要%s。',id,expected);
end
