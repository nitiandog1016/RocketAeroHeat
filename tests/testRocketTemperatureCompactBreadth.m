function report=testRocketTemperatureCompactBreadth(folder)
%TESTROCKETTEMPERATURECOMPACTBREADTH 检查17及33个独立条件而非仅精确重复的扩展容量。
% 模块一：使用真实小域QR映射和原C/W，制造具有明确秩的代数H，不冒充物理温度行。
% 模块二：两种CV排列均核对完整支撑、域外坐标逐值不变及独立稠密受限解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
s=load(fullfile(root,'outputs/v89_validation/unit01/legacy/temperature_duplicate_unit.mat'),'cases');cases=cell(4,1);index=0;
for j=1:2
    old=s.cases{j};n=old.factor.cellCount;N=5*n;
    for q=[17 33]
        assert(N>q);index=index+1;columns=mod((0:q-1)*7,N)+1;assert(numel(unique(columns))==q);
        base=sparse((1:q)',columns',1+(1:q)'/q,q,N);mapping=[(1:q)';1;3;1;5];
        rows=struct('kind',old.rows.kind,'cellCount',n,'faceIDs',(1:numel(mapping))','H',base(mapping,:));
        projector=buildRocketTemperatureProjector(old.factor,rows);
        assert(projector.independentConstraintCount==q&&projector.constraintCount==numel(mapping));
        Pi=zeros(N);
        for k=1:N
            e=zeros(n,5);e(k)=1;v=applyRocketTemperatureProjector(projector,e);Pi(:,k)=v(:);
        end
        assert(norm(Pi-Pi','fro')<1e-12&&norm(Pi*Pi-Pi,'fro')<1e-12&&norm(full(rows.H)*old.P*Pi,'fro')<1e-9);
        input=reshape(sin((1:N)'/sqrt(2)),n,5);output=applyRocketTemperatureProjector(projector,input);
        outside=setdiff((1:N)',projector.basisCoordinateIndices);assert(isequal(output(outside),input(outside)));
        rhs=old.source.base.rhs;A=[old.source.base.C;old.factor.penaltyInOriginalResidualUnits*old.source.W];
        Z=null(full(rows.H));dense=Z*((A*Z)\[rhs(:);zeros(N,1)]);
        [value,detail,coordinates]=solveRocketTemperatureProjectedDirection(old.packed,old.factor,rhs,old.source.metric,projector,rows,900);
        difference=norm(value(:)-dense)/max(1,norm(dense));assert(difference<1e-8&&detail.adjointIdentityRelativeDifference<1e-12);
        dependent=rows;dependent.H=[rows.H;base(1,:)+base(2,:)];dependent.faceIDs=(1:size(dependent.H,1))';
        rejected=false;try,buildRocketTemperatureProjector(old.factor,dependent);catch e,rejected=strcmp(e.identifier,'TemperatureProjection:Rank');end
        assert(rejected);
        cases{index}=struct('source',old.source,'factor',old.factor,'rows',rows,'projector',projector,'P',old.P,'Pi',Pi, ...
            'value',value,'detail',detail,'coordinates',coordinates,'dense',dense,'denseDifference',difference, ...
            'nonduplicateDependenceRejected',rejected,'outsideCoordinatesExactlyUnchanged',true,'manufacturedAlgebraicRows',true);
    end
end
report=struct('status','passed-two-orderings-seventeen-and-thirty-three-independent-compact-conditions', ...
    'manufacturedAlgebraicRowsNotPhysicalTemperatures',true,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_capacity_unit.mat'),'cases','report','-v7');disp(report);
end
