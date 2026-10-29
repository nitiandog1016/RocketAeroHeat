function report=testRocketTemperatureRowConditioning(folder)
%TESTROCKETTEMPERATUREROWCONDITIONING 验证可逆差值表示保留全部条件及原秩/方向门槛。
% 模块一：二进制近邻制造行先重现旧QR拒绝，再核对新表示未删条件且可精确重构。
% 模块二：两种排列分别与全部等价条件的稠密零空间最小二乘比较，不冒充物理NS。
% 模块三：真正相关的倍数行仍拒绝；旧温域/完整NS/秩门槛不作放宽。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
p=load(rocketTemperatureResearchFixture('parent_least_squares_unit_evidence.mat'),'cases');
w=load(rocketTemperatureResearchFixture('peak_weighted_unit_evidence.mat'),'cases');cases=cell(2,1);
for j=1:2
    a=p.cases{j};b=w.cases{j};n=a.packed.cellCount;
    factor=buildRocketAugmentedPatchQR(a.packed,a.factor,b.metric,b.detail.penaltyInOriginalResidualUnits);
    first=sparse(1,1+(0:4)*n,[1 2 4 8 16],1,5*n);second=first;second(1,1+4*n)=16+eps(16);
    H=[first;second];original=H;
    rows=struct('kind','current-interior-selected-Riemann-sample-temperature-tangent-rows', ...
        'cellCount',n,'faceIDs',[1;2],'H',H);
    [working,conditioning]=conditionRocketTemperatureRows(H);
    assert(isequal(H,original)&&isequal(conditioning.transformedRowIndices,2)&& ...
        isequal(conditioning.referenceRowIndices,1)&&~conditioning.rankBasedRowSelection);
    restored=working(1,:)+pow2(full(working(2,:)),-conditioning.powerTwoExponents);
    assert(isequal(full(restored),full(H(2,:))));
    oldColumns=zeros(5*n,2);
    for k=1:2
        mapped=applyRocketAugmentedPatchQR(factor,reshape(full(H(k,:))',n,5),'transp');oldColumns(:,k)=mapped(:);
    end
    [~,oldR]=qr(oldColumns./vecnorm(oldColumns,2,1),0);oldPivots=abs(diag(oldR));
    assert(any(oldPivots<=1e-10),'TemperatureConditioning:Manufacture','制造行必须实际重现原秩门槛拒绝。');
    projector=buildRocketTemperatureProjector(factor,rows);
    assert(projector.constraintCount==2&&projector.independentConstraintCount==2&& ...
        all(projector.normalizedQRDiagonal>1e-10)&&~projector.rankTruncationUsed);
    assert(max(abs(projector.constraintColumnNorms-vecnorm(oldColumns,2,1))./vecnorm(oldColumns,2,1))<1e-12);
    [value,detail,coordinates]=solveRocketTemperatureProjectedDirection(a.packed,factor,b.base.rhs,b.metric,projector,rows,900);
    A=[b.base.C;factor.penaltyInOriginalResidualUnits*b.W];Z=null(full(working));
    dense=Z*((A*Z)\[b.base.rhs(:);zeros(5*n,1)]);
    difference=norm(value(:)-dense)/max(1,norm(dense));assert(difference<1e-8);
    dependent=rows;dependent.H=[first;2*first];rejected=false;
    try,buildRocketTemperatureProjector(factor,dependent);
    catch e,rejected=strcmp(e.identifier,'TemperatureProjection:Rank');end
    assert(rejected,'TemperatureConditioning:Rank','真正相关的不同原行仍由原1e-10门槛拒绝。');
    cases{j}=struct('H',H,'workingRows',working,'conditioning',conditioning,'oldNormalizedQRDiagonal',oldPivots, ...
        'projector',projector,'value',value,'detail',detail,'coordinates',coordinates,'C',b.base.C,'W',b.W, ...
        'penalty',factor.penaltyInOriginalResidualUnits,'rhs',b.base.rhs,'denseDifference',difference, ...
        'proportionalConditionsStillRejected',rejected);
end
report=struct('status','passed-two-orderings-exact-difference-all-conditions-dense-direction-and-unchanged-rank-refusal', ...
    'manufacturedRowsNotPhysicalRiemannDerivatives',true,'newActualNSEvaluations',0, ...
    'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_row_conditioning_unit.mat'),'cases','report','-v7');disp(report);
end
