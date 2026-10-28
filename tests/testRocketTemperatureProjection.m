function report=testRocketTemperatureProjection(folder)
%TESTROCKETTEMPERATUREPROJECTION 校核真实面解析温度导数、完整伴随和独立稠密受限解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'),fullfile(root,'tests','fixtures'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
s=load(fullfile(root,'outputs/v83_validation/riemannreservoirunit02/riemann_reservoir_unit_evidence.mat'),'result','cut','gas');
Q=s.result.Q;gas=s.gas;scale=s.result.residualScale;cut=s.cut;f=cut.openFaces;
face=find(f.area>0&f.cellMinus>0&f.cellPlus>0,1);rows=buildRocketInterfaceTemperatureRows(Q,cut,gas,scale,face);
% 原内核多面批次四个完整输出逐值回放；新单面输出应与重复批次首面完全相同。
ql=repmat(rows.QL,2,1);qr=repmat(rows.QR,2,1);nn=repmat(rows.normal,2,1);
[L,R,A,S]=cutCellExactEulerBranchJacobian(ql,qr,nn,gas);
[oldL,oldR,oldA,oldS]=rocketExactEulerBranchV85Reference(ql,qr,nn,gas);
assert(isequaln(L,oldL)&&isequaln(R,oldR)&&isequaln(A,oldA)&&isequaln(S,oldS));
[oneL,oneR,oneA,oneS]=cutCellExactEulerBranchJacobian(rows.QL,rows.QR,rows.normal,gas);
assert(isequaln(oneL,L(:,:,1))&&isequaln(oneR,R(:,:,1))&&isequaln(oneS.left,S.left(:,:,1))&&isequaln(oneS.right,S.right(:,:,1)));
assert(isequaln(oneA.baseFlux,A.baseFlux(1,:)));
fd=zeros(1,5*size(Q,1));
for side=1:2
    ids=[rows.ownerCV rows.neighbourCV];cellID=ids(side);
    for b=1:5
        h=1e-6*max(1,abs(Q(cellID,b)/scale(b)));qp=Q;qm=Q;
        qp(cellID,b)=qp(cellID,b)+h*scale(b);qm(cellID,b)=qm(cellID,b)-h*scale(b);
        [~,ap]=cutCellExactEulerFlux(qp(rows.ownerCV,:),qp(rows.neighbourCV,:),rows.normal,gas);
        [~,am]=cutCellExactEulerFlux(qm(rows.ownerCV,:),qm(rows.neighbourCV,:),rows.normal,gas);
        fd(cellID+(b-1)*size(Q,1))=(ap.samplePrimitive(5)-am.samplePrimitive(5))/(2*h);
    end
end
fdError=norm(fd-full(rows.H))/max(1,norm(full(rows.H)));assert(fdError<1e-7);
p=load(fullfile(root,'outputs/v83_validation/parentlsqrunit02/parent_least_squares_unit_evidence.mat'),'cases');
m=load(fullfile(root,'outputs/v83_validation/peakweightunit01/peak_weighted_unit_evidence.mat'),'cases');cases=cell(2,1);
for j=1:2
    a=p.cases{j};b=m.cases{j};n=a.packed.cellCount;N=5*n;assert(rows.cellCount==n);
    factor=buildRocketAugmentedPatchQR(a.packed,a.factor,b.metric,b.detail.penaltyInOriginalResidualUnits);
    projector=buildRocketTemperatureProjector(factor,rows);P=zeros(N);Pi=P;
    for k=1:N
        e=zeros(n,5);e(k)=1;v=applyRocketAugmentedPatchQR(factor,e);P(:,k)=v(:);
        v=applyRocketTemperatureProjector(projector,e);Pi(:,k)=v(:);
    end
    assert(norm(Pi-Pi','fro')<1e-12&&norm(Pi*Pi-Pi,'fro')<1e-12);
    assert(norm(full(rows.H)*P*Pi)<1e-10*norm(full(rows.H)*P));
    rhs=b.base.rhs;A=[b.base.C;factor.penaltyInOriginalResidualUnits*b.W];augmented=[rhs(:);zeros(N,1)];
    Z=null(full(rows.H));dense=Z*((A*Z)\augmented);
    [value,detail,coordinates]=solveRocketTemperatureProjectedDirection(a.packed,factor,rhs,b.metric,projector,rows,900);
    error=norm(value(:)-dense)/max(1,norm(dense));assert(error<1e-8&&detail.adjointIdentityRelativeDifference<1e-12);
    cases{j}=struct('source',b,'packed',a.packed,'factor',factor,'projector',projector,'rows',rows,'P',P,'Pi',Pi, ...
        'value',value,'detail',detail,'coordinates',coordinates,'dense',dense,'denseDifference',error);
end
report=struct('status','passed-real-interface-temperature-FD-two-orderings-symmetric-projector-and-dense-constrained-objective', ...
    'temperatureRowFDRelativeDifference',fdError,'oldCompleteMultiFaceBranchOutputsExactlyReplayed',true, ...
    'singleFaceExactlyMatchesDuplicateBatch',true,'physicalBoundaryConditionAdded',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_projection_unit.mat'),'cases','report','Q','gas','scale','rows','fd','-v7');disp(report);
end
