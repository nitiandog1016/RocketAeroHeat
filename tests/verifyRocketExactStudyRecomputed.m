function report=verifyRocketExactStudyRecomputed(newRoot,frozenRoot)
%VERIFYROCKETEXACTSTUDYRECOMPUTED 独立工程包内重算后逐字段核对完整数值对象和原面JSON。
% 模块一：精确通量、三维集成及两条各三次Krylov完整结果均须严格一致。
% 模块二：仅排除真实机器耗时和引擎PID；保留Q、所有通量/Jvp/失败候选及计数。
% 模块三：固定方向合法/越域试算与缺省完整源回放分别核对，不授予物理面积资格。
checked=0;
for spec={{'flux02','exact_flux_evidence.mat'}, {'integrate01','exact_integration_evidence.mat'}, ...
        {'exact_krylov01','steady_krylov_evidence.mat'}, {'primitive_krylov01','steady_krylov_evidence.mat'}, ...
        {'primitive02','primitive_path_evidence.mat'}, {'trial_compatibility01','trial_compatibility_evidence.mat'}}
    item=spec{1};file=fullfile(item{1},item{2});a=load(fullfile(newRoot,file));b=load(fullfile(frozenRoot,file));
    assert(isequal(sort(fieldnames(a)),sort(fieldnames(b))),'ExactReplay:Fields','保存变量集合不同：%s',file);
    for name=fieldnames(a)'
        key=name{1};va=a.(key);vb=b.(key);
        if isstruct(va)&&isfield(va,'history')&&ismember(key,{'result','reference','actual','selected'})
            va=untimed(va);vb=untimed(vb);
        end
        if strcmp(key,'report'),va=reportUntimed(va);vb=reportUntimed(vb);end
        assert(isequaln(va,vb),'ExactReplay:Identity','完整数值对象未精确复现：%s / %s',file,key);checked=checked+1;
    end
    clear a b va vb;
end
jsonCount=0;
for folder={'workbench02','workbench03'}
    ids={'hllc-source','exact-baseline','exact-krylov'};
    if strcmp(folder{1},'workbench02'),ids{3}='exact-primitive-krylov';end
    for id=ids
        file=fullfile(folder{1},'cases',id{1},'results.json');
        assert(strcmp(fileread(fullfile(newRoot,file)),fileread(fullfile(frozenRoot,file))), ...
            'ExactReplay:View','完整原面结果JSON未逐字节复现：%s',file);jsonCount=jsonCount+1;
    end
end
report=struct('status','passed-isolated-complete-exact-flux-and-two-real-three-correction-chain-identity', ...
    'completeNumericalObjectsChecked',checked,'resultJSONExactlyReproduced',jsonCount, ...
    'actualNonlinearCorrectionsRecomputed',6,'oldHLLCSixCorrectionsRecomputed',false, ...
    'runtimeAndEnginePIDOnlyExcluded',true,'physicalTimeAdvanced',0,'steadyConverged',false, ...
    'physicalSolution',false,'wholeGoalComplete',false,'newProductionSolvedAreaM2',0);
fid=fopen(fullfile(newRoot,'recomputed_identity.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function r=untimed(r)
r.history=rmfield(r.history,{'operatorBuildSeconds','linearSolveSeconds'});
end
function r=reportUntimed(r)
for name={'elapsedSeconds','engineProcessID'}
    if isfield(r,name{1}),r=rmfield(r,name{1});end
end
if isfield(r,'rows')&&isstruct(r.rows)&&isfield(r.rows,'elapsedSeconds'),r.rows=rmfield(r.rows,'elapsedSeconds');end
end
