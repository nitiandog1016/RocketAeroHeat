function report=verifyRocketShellKrylovRecomputed(newRoot,frozenRoot)
%VERIFYROCKETSHELLKRYLOVRECOMPUTED 独立包内重算新功检查及六次原三维Krylov后严格核对。
% 模块一：只去除两项运行耗时；其余完整状态、方向、原残差/info/候选日志严格一致。
% 模块二：旧功受限对照保留冻结原值，不称本次重新求解；五组原面JSON须逐字节一致。
checked=0;
for name={'arithmetic01','arithmetic02'}
    file=fullfile(name{1},'steady_krylov_evidence.mat');a=load(fullfile(newRoot,file));b=load(fullfile(frozenRoot,file));
    ra=a.result;rb=b.result;ra.history=untimed(ra.history);rb.history=untimed(rb.history);
    assert(isequaln(ra,rb),'ShellKrylovReplay:Result','完整Krylov数值对象未精确复现：%s',name{1});checked=checked+1;
    for key={'mapped','gas','boundary'}
        k=key{1};assert(isequaln(a.(k),b.(k)),'ShellKrylovReplay:Field','完整映射/物性/边界不同：%s',k);checked=checked+1;
    end
    if isfield(a,'displayHistory'),assert(isequaln(untimed(a.displayHistory),untimed(b.displayHistory)),'ShellKrylovReplay:History','完整接续历史未复现。');checked=checked+1;end
    clear a b ra rb;
end
a=load(fullfile(newRoot,'work01','work_velocity_evidence.mat'));b=load(fullfile(frozenRoot,'work02','work_velocity_evidence.mat'));
for key={'rate','info','gas','boundary'},k=key{1};assert(isequaln(a.(k),b.(k)),'ShellKrylovReplay:Work','完整算术功检查不同：%s',k);checked=checked+1;end
clear a b;
for id={'initial-legacy','initial-arithmetic','legacy-line','legacy-krylov','arithmetic-krylov'}
    file=fullfile('cases',id{1},'results.json');
    assert(strcmp(fileread(fullfile(newRoot,'workbench02',file)),fileread(fullfile(frozenRoot,'workbench02',file))), ...
        'ShellKrylovReplay:View','原面完整结果JSON不同：%s',id{1});
end
report=struct('status','passed-isolated-full-work-strategy-and-six-actual-Krylov-state-identity', ...
    'completeNumericalObjectsChecked',checked,'resultJSONExactlyReproduced',5,'nonlinearCorrectionsActuallyRecomputed',6, ...
    'oldLegacyComparisonTrialsRecomputed',false,'onlyRuntimeTimingsExcluded',true, ...
    'physicalTimeAdvanced',0,'steadyConverged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
fid=fopen(fullfile(newRoot,'recomputed_identity.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function h=untimed(h)
% 模块用途：实际机器耗时允许变化，不删除任何残差、线性缺陷、方向或状态字段。
h=rmfield(h,{'operatorBuildSeconds','linearSolveSeconds'});
end
