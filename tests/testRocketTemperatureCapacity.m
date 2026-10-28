function report=testRocketTemperatureCapacity(folder)
%TESTROCKETTEMPERATURECAPACITY 验收超16行紧凑投影、全约束和完整原目标。
% 模块一：旧七行/三条件及完整秩逐值回放；扩展到17行时保留所有标签。
% 模块二：两种原CV排列下比较完整P、紧凑投影、伴随和独立稠密约束解。
% 模块三：真实Q34全部17行与旧分批证据逐字段核对；不执行新NS或授予面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
legacy=testRocketTemperatureDuplicateRows(fullfile(folder,'legacy'));
s=load(fullfile(folder,'legacy/temperature_duplicate_unit.mat'),'cases');cases=cell(2,1);
for j=1:2
    old=s.cases{j};n=old.factor.cellCount;N=5*n;map=mod((0:16)',7)+1;
    rows=old.rows;rows.H=rows.H(map,:);rows.faceIDs=(1:17)';
    projector=buildRocketTemperatureProjector(old.factor,rows);
    assert(projector.constraintCount==17&&projector.independentConstraintCount==3&& ...
        isfield(projector,'basisCoordinateIndices')&&projector.additionalConstructionArrayBoundBytes<=448*1024^2);
    Pi=zeros(N);
    for k=1:N
        e=zeros(n,5);e(k)=1;v=applyRocketTemperatureProjector(projector,e);Pi(:,k)=v(:);
    end
    assert(norm(Pi-old.Pi,'fro')<1e-12&&norm(Pi-Pi','fro')<1e-12&&norm(Pi*Pi-Pi,'fro')<1e-12);
    rhs=old.source.base.rhs;
    [value,detail,coordinates]=solveRocketTemperatureProjectedDirection(old.packed,old.factor,rhs,old.source.metric,projector,rows,900);
    difference=norm(value(:)-old.dense)/max(1,norm(old.dense));assert(difference<1e-8);
    dependent=rows;dependent.H=[rows.H;rows.H(1,:)+rows.H(2,:)];dependent.faceIDs=(1:18)';
    rejected=false;try,buildRocketTemperatureProjector(old.factor,dependent);catch e,rejected=strcmp(e.identifier,'TemperatureProjection:Rank');end
    assert(rejected,'TemperatureCapacity:Rank','非重复相关约束不能被紧凑表示截断。');
    cases{j}=struct('source',old.source,'factor',old.factor,'rows',rows,'projector',projector,'P',old.P,'Pi',Pi, ...
        'value',value,'detail',detail,'coordinates',coordinates,'dense',old.dense,'denseDifference',difference, ...
        'nonduplicateDependenceRejected',rejected);
end
% 使用真实当前Q34检查选择器和分批行，而非沿用Q33缓存。
stateFile=fullfile(root,'outputs/v83_validation/krylov34/wall_velocity_krylov_evidence.mat');
geometry=fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat');
sourceStateSHA256=rocketStreamingFileSHA256(stateFile);geometrySHA256=rocketStreamingFileSHA256(geometry);
q=load(stateFile,'result','gas');g=load(geometry,'domain');
realRows=selectRocketTemperatureBandRows(q.result.Q,q.result.info,g.domain,q.gas,q.result.residualScale);
capacity=load(fullfile(root,'outputs/v89_validation/capacity01/temperature_capacity.mat'),'rows');
assert(numel(realRows.faceIDs)==17);
for key=fieldnames(capacity.rows)'
    assert(isequaln(realRows.(key{1}),capacity.rows.(key{1})),'TemperatureCapacity:Rows','真实行字段%s与旧分批结果不同。',key{1});
end
rejected=false;try
    selectRocketTemperatureBandRows(q.result.Q,q.result.info,g.domain,q.gas,q.result.residualScale,struct('maximumConstraintRows',16));
catch e,rejected=strcmp(e.identifier,'TemperatureActive:Budget');end
assert(rejected,'TemperatureCapacity:Budget','用户显式设置的行上限必须保持拒绝。');
report=struct('status','passed-extended-seventeen-row-compact-projector-and-current-Q34-row-replay', ...
    'legacy',legacy,'realRowsExactlyReplayed',true,'realOriginalFaceCount',17,'explicitRowLimitStillRejected',rejected, ...
    'sourceStateSHA256',sourceStateSHA256,'geometrySHA256',geometrySHA256,'rankThresholdUnchanged',true, ...
    'constructionBudgetBytes',448*1024^2,'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_capacity_unit.mat'),'cases','realRows','report','-v7');disp(report);
end
