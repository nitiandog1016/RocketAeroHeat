function report=testRocketTemperatureDuplicateRows(folder)
%TESTROCKETTEMPERATUREDUPLICATEROWS 检查精确重复约束的等价表示及全部原行不变。
% 模块一：旧完整秩投影逐字段、逐值回放，重复行只复用完全相同的条件。
% 模块二：七个不同面标签对应三个条件；全七行均验收，不截断近似相关行。
% 模块三：两种CV排列下核对投影、伴随与独立稠密最终子空间解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'),fullfile(root,'tests','fixtures'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
s=load(fullfile(root,'outputs/v87_validation/unit01/temperature_active_unit.mat'),'cases');cases=cell(2,1);
duplicates=[1;2;1;3;2;3;1];
for j=1:2
    old=s.cases{j};factor=old.factor;n=factor.cellCount;N=5*n;
    % 原小域的全部三行本就含重复；旧算法当时只激活其中独立的两行。
    % 取这两行，再加一个不在其支撑内的独立代数条件，构成明确的制造三条件校核。
    baseRows=old.rows;baseRows.H=old.rows.H(old.detail.activeRowIndices,:);assert(size(baseRows.H,1)==2);
    column=find(full(sum(abs(baseRows.H),1))==0,1);assert(~isempty(column));
    baseRows.H=[baseRows.H;sparse(1,column,1,1,N)];baseRows.faceIDs=(1:3)';
    reference=rocketTemperatureProjectorV87Reference(factor,baseRows);
    current=buildRocketTemperatureProjector(factor,baseRows);
    for key=fieldnames(reference)',assert(isequaln(current.(key{1}),reference.(key{1})));end
    rows=struct('kind',old.rows.kind,'cellCount',n,'H',baseRows.H(duplicates,:),'faceIDs',(1:7)');
    projector=buildRocketTemperatureProjector(factor,rows);
    assert(projector.constraintCount==7&&projector.independentConstraintCount==3&& ...
        isequal(projector.originalToUniqueRow,duplicates)&&projector.allOriginalRowsExactlyRepresented&&~projector.rankTruncationUsed);
    assert(isequaln(projector.basis,reference.basis));
    P=zeros(N);Pi=P;
    for k=1:N
        e=zeros(n,5);e(k)=1;v=applyRocketAugmentedPatchQR(factor,e);P(:,k)=v(:);
        v=applyRocketTemperatureProjector(projector,e);Pi(:,k)=v(:);
    end
    assert(norm(Pi-Pi','fro')<1e-12&&norm(Pi*Pi-Pi,'fro')<1e-12);
    assert(norm(full(rows.H)*P*Pi,'fro')<=1e-10*norm(full(rows.H)*P,'fro'));
    rhs=old.source.base.rhs;A=[old.source.base.C;factor.penaltyInOriginalResidualUnits*old.source.W];
    Z=null(full(rows.H));dense=Z*((A*Z)\[rhs(:);zeros(N,1)]);
    [value,detail,coordinates]=solveRocketTemperatureProjectedDirection(old.packed,factor,rhs,old.source.metric,projector,rows,900);
    difference=norm(value(:)-dense)/max(1,norm(dense));assert(difference<1e-8&&detail.adjointIdentityRelativeDifference<1e-12);
    dependent=rows;dependent.H=[rows.H;rows.H(1,:)+rows.H(2,:)];dependent.faceIDs=(1:8)';
    rejected=false;try,buildRocketTemperatureProjector(factor,dependent);catch e,rejected=strcmp(e.identifier,'TemperatureProjection:Rank');end
    assert(rejected,'TemperatureDuplicate:Rank','非重复的相关条件仍须按原秩门槛拒绝。');
    cases{j}=struct('rows',rows,'source',old.source,'packed',old.packed,'factor',factor, ...
        'projector',projector,'reference',reference,'P',P,'Pi',Pi,'value',value,'detail',detail, ...
        'coordinates',coordinates,'dense',dense,'denseDifference',difference,'nonduplicateDependenceRejected',rejected);
end
report=struct('status','passed-exact-seven-row-three-condition-projection-and-legacy-full-rank-replay', ...
    'allOriginalConstraintRowsRetained',true,'onlyBitwiseIdenticalRowsShareRepresentation',true, ...
    'oldFullRankProjectorFieldsExactlyReplayed',true,'rankThresholdUnchanged',true,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_duplicate_unit.mat'),'cases','report','duplicates','-v7');disp(report);
end
