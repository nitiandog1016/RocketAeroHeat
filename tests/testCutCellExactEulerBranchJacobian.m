function report=testCutCellExactEulerBranchJacobian(folder)
%TESTCUTCELLEXACTEULERBRANCHJACOBIAN 当前波系解析导数对原通量中心差分及不变量校核。
% 模块一：覆盖六种采样区、激波/稀疏波/接触和实际近壁受限面，不以人工场替代真实点。
% 模块二：每侧五量独立扰动；跨分支点明确单列，不将分段导数称为全局光滑导数。
% 模块三：Euler齐次性、法向反转、三维旋转提供额外约束；不证明完整NS或物理精度。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'ExactBranchTest:Overwrite','不得覆盖实际导数校核。');mkdir(folder);
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
s=load(fullfile(root,'outputs','v83_validation','krylov03','wall_velocity_krylov_evidence.mat'),'result','gas');gas=s.gas;
PL=[];PR=[];
for shift=[0 100 700]
    PL=[PL;1 shift 50 -20 300;.125 -shift -30 70 240]; %#ok<AGROW>
    PR=[PR;.125 shift -30 70 240;1 -shift 50 -20 300]; %#ok<AGROW>
end
PL=[PL;1 670 15 -20 300;1 -520 20 40 300;1.2 1000 40 -30 400];
PR=[PR;1 -530 -30 50 300;1 680 -35 -10 300;1.2 1000 40 -30 400];
QL=cons(PL,gas);QR=cons(PR,gas);normal=repmat([1 0 0],size(PL,1),1);syntheticCount=size(PL,1);
actualFaces=[16397 16401 16403 19001 70995 19500]';f=d.domain.openFaces;
assert(all(f.cellPlus(actualFaces)>0));
QL=[QL;s.result.Q(f.cellMinus(actualFaces),:)];QR=[QR;s.result.Q(f.cellPlus(actualFaces),:)];normal=[normal;f.normal(actualFaces,:)];clear s d;
[left,right,audit]=cutCellExactEulerBranchJacobian(QL,QR,normal,gas);[~,base]=cutCellExactEulerFlux(QL,QR,normal,gas);
assert(isequal(unique(base.sampleRegion(1:syntheticCount)),(1:6)'),'ExactBranchTest:Coverage','六种采样波区未全部覆盖。');
N=size(QL,1);worst=0;tested=0;excluded=0;fdErrors=zeros(N,10);h=1e-6;
for side=1:2
    Q=QL;if side==2,Q=QR;end
    reference=left;if side==2,reference=right;end
    rho=Q(:,1);v=Q(:,2:4)./rho;T=(Q(:,5)-.5*rho.*sum(v.^2,2))./(rho*gas.R/(gas.gamma-1));
    speed=max(vecnorm(v,2,2),sqrt(gas.gamma*gas.R*T));basis=[rho,rho.*speed.*ones(1,3),rho.*speed.^2];
    for b=1:5
        step=h*basis(:,b);plus=Q;minus=Q;plus(:,b)=plus(:,b)+step;minus(:,b)=minus(:,b)-step;
        if side==1,[fp,ap]=cutCellExactEulerFlux(plus,QR,normal,gas);[fm,am]=cutCellExactEulerFlux(minus,QR,normal,gas);
        else,[fp,ap]=cutCellExactEulerFlux(QL,plus,normal,gas);[fm,am]=cutCellExactEulerFlux(QL,minus,normal,gas);end
        mask=ap.sampleRegion==base.sampleRegion&am.sampleRegion==base.sampleRegion;
        actual=(fp-fm)./(2*step);expected=reshape(reference(:,b,:),5,N)';
        error=max(abs(actual-expected),[],2)./max(1,max(abs(expected),[],2));
        fdErrors(:,(side-1)*5+b)=error;worst=max(worst,max(error(mask)));tested=tested+nnz(mask);excluded=excluded+nnz(~mask);
    end
end
assert(tested>=syntheticCount*10&&worst<5e-7,'ExactBranchTest:Derivative','原通量中心差分不符合当前分支导数。');
homogeneous=zeros(N,5);
for b=1:5,homogeneous=homogeneous+reshape(left(:,b,:),5,N)'.*QL(:,b)+reshape(right(:,b,:),5,N)'.*QR(:,b);end
homogeneity=max(abs(homogeneous-audit.baseFlux),[],'all')/max(1,max(abs(audit.baseFlux),[],'all'));
assert(homogeneity<2e-12,'ExactBranchTest:Homogeneity','Euler一阶齐次恒等式不符合。');
[reverseL,reverseR]=cutCellExactEulerBranchJacobian(QR,QL,-normal,gas);
mask=abs(audit.contactNormalVelocityMs)>1e-8;
reverse=max([max(abs(left(:,:,mask)+reverseR(:,:,mask)),[],'all'),max(abs(right(:,:,mask)+reverseL(:,:,mask)),[],'all')])/ ...
    max(1,max(abs(cat(2,left(:,:,mask),right(:,:,mask))),[],'all'));
assert(reverse<2e-12,'ExactBranchTest:Reverse','反法向/交换状态导数不符合。');
angle=.47;rotation=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];S=eye(5);S(2:4,2:4)=rotation;
ql=QL;qr=QR;ql(:,2:4)=QL(:,2:4)*rotation';qr(:,2:4)=QR(:,2:4)*rotation';
[rotL,rotR]=cutCellExactEulerBranchJacobian(ql,qr,normal*rotation',gas);rotationError=0;
for j=1:N
    expectedL=S*left(:,:,j)*S';expectedR=S*right(:,:,j)*S';
    rotationError=max(rotationError,max(abs([rotL(:,:,j)-expectedL,rotR(:,:,j)-expectedR]),[],'all')/max(1,max(abs([expectedL expectedR]),[],'all')));
end
assert(rotationError<3e-12,'ExactBranchTest:Rotation','三维旋转导数不符合。');
bad=load(fullfile(root,'outputs','v83_validation','rejected02','rejected_interface_evidence.mat'),'ql','qr','normal','gas');
rejected=false;try,cutCellExactEulerBranchJacobian(bad.ql,bad.qr,bad.normal,bad.gas);catch ex,rejected=strcmp(ex.identifier,'CutExactEuler:Temperature');end
assert(rejected,'ExactBranchTest:Guard','真实78.65K界面必须保持拒绝。');
report=struct('status','passed-six-sampling-regions-real-wall-faces-and-selected-branch-derivatives', ...
    'caseCount',N,'syntheticCases',syntheticCount,'actualFaces',actualFaces,'finiteDifferenceColumnsTested',tested, ...
    'crossBranchDifferenceColumnsExcluded',excluded,'maximumBranchFDRelativeDifference',worst, ...
    'homogeneityRelativeDifference',homogeneity,'reverseNormalRelativeDifference',reverse,'rotationRelativeDifference',rotationError, ...
    'originalRejectedInterfaceStillRejected',rejected,'globalSmoothJacobianClaimed',false,'fullNSJacobian',false, ...
    'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'exact_branch_jacobian_evidence.mat'),'report','QL','QR','normal','gas','left','right','fdErrors','audit','-v7');
fid=fopen(fullfile(folder,'exact_branch_jacobian_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function Q=cons(P,gas)
Q=[P(:,1),P(:,1).*P(:,2:4),P(:,1).*(gas.R/(gas.gamma-1)*P(:,5)+.5*sum(P(:,2:4).^2,2))];
end
