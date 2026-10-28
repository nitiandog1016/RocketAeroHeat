function report=testRiemannReservoirSelectedDirection(folder)
%TESTRIEMANNRESERVOIRSELECTEDDIRECTION 校核真实波系采样导数及完整新储库NS方向链。
% 模块一：旧三输出逐值保持；第四输出150列对原采样中心差分，六波区全部覆盖。
% 模块二：32CV新储库完整NS逐值回放；缓存/非缓存解析方向和两种原NS差分步核对。
% 模块三：只校核当前选中分支；不把方向一致性当作稳态、实验或合格面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
old=load(fullfile(root,'outputs','v83_validation','branchjac01','exact_branch_jacobian_evidence.mat'),'QL','QR','normal','gas','left','right','audit');
QL=old.QL;QR=old.QR;normal=old.normal;gas=old.gas;
[left,right,audit,sample]=cutCellExactEulerBranchJacobian(QL,QR,normal,gas);
assert(isequaln(left,old.left)&&isequaln(right,old.right)&&isequaln(audit,old.audit),'RiemannDirection:OldOutputs','新增采样导数改变了旧三输出。');
[l3,r3,a3]=cutCellExactEulerBranchJacobian(QL,QR,normal,gas);assert(isequal(l3,left)&&isequal(r3,right)&&isequaln(a3,audit));
N=size(QL,1);errors=zeros(N,10);tested=0;excluded=0;worst=0;h=1e-6;
for side=1:2
    q=QL;reference=sample.left;if side==2,q=QR;reference=sample.right;end
    rho=q(:,1);u=q(:,2:4)./rho;T=(q(:,5)-.5*rho.*sum(u.^2,2))./(rho*gas.R/(gas.gamma-1));
    speed=max(vecnorm(u,2,2),sqrt(gas.gamma*gas.R*T));basis=[rho,rho.*speed.*ones(1,3),rho.*speed.^2];
    for b=1:5
        step=h*basis(:,b);plus=q;minus=q;plus(:,b)=plus(:,b)+step;minus(:,b)=minus(:,b)-step;
        if side==1,[~,p]=cutCellExactEulerFlux(plus,QR,normal,gas);[~,m]=cutCellExactEulerFlux(minus,QR,normal,gas);
        else,[~,p]=cutCellExactEulerFlux(QL,plus,normal,gas);[~,m]=cutCellExactEulerFlux(QL,minus,normal,gas);end
        keep=p.sampleRegion==sample.sampleRegion&m.sampleRegion==sample.sampleRegion;
        expected=reshape(reference(:,b,:),5,N)';actual=(p.samplePrimitive-m.samplePrimitive)./(2*step);
        e=max(abs(actual-expected),[],2)./max(1,max(abs(expected),[],2));errors(:,(side-1)*5+b)=e;
        tested=tested+nnz(keep);excluded=excluded+nnz(~keep);worst=max(worst,max(e(keep)));
    end
end
assert(tested==150&&excluded==0&&worst<5e-7&&numel(unique(sample.sampleRegion))==6);
s=load(fullfile(root,'outputs','v83_validation','riemannreservoirunit02','riemann_reservoir_unit_evidence.mat'),'Q','rate','info','cut','gas','boundary','result');
Q=s.Q;cut=s.cut;gas=s.gas;boundary=s.boundary;scale=s.result.residualScale;
[rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);assert(isequaln(rate,s.rate)&&isequaln(info,s.info));
direction=reshape(.1*sin((1:numel(Q))*.37),size(Q));cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);
[value,detail]=applyCutCellSelectedDirection(cache,direction);
state=cutCellBoundaryStencilState(Q,cut,gas,boundary);bd=cutCellBoundaryDirectionalState(Q,direction,cut,gas,boundary,scale,state);
G=buildCutCellGradientOperator(cut,state);
uncached=cutCellWlsViscousDirection(Q,direction,info,cut,gas,boundary,scale,G,state,bd)+ ...
    cutCellSelectedEulerDirection(Q,direction,info,cut,gas,boundary,scale,state,bd);
cacheDifference=norm(value-uncached,'fro')/max(1,norm(value,'fro'));assert(cacheDifference<1e-12);
config=struct('differenceStep',1e-6,'scheme','forward','saveEvidence',true,'spatialConservationTolerance',1e-9);
[fd,fdEvidence]=rocketShellResidualDirectionalDerivative(Q,direction,rate,cut,gas,boundary,scale,config);
config.differenceStep=5e-7;[half,halfEvidence]=rocketShellResidualDirectionalDerivative(Q,direction,rate,cut,gas,boundary,scale,config);
differences=[norm(value-fd,'fro') norm(value-half,'fro')]/max(1,norm(value,'fro'));
assert(all(differences<2e-6),'RiemannDirection:FullNS','完整解析方向与原NS差分不一致。');
report=struct('status','passed-six-region-sample-Jacobians-unchanged-old-outputs-and-full-new-reservoir-NS-chain', ...
    'oldThreeOutputsExactlyReplayed',true,'sampleFDColumns',tested,'crossBranchColumnsExcluded',excluded,'maximumSampleFDDifference',worst, ...
    'cachedUncachedRelativeDifference',cacheDifference,'fullNSDirectionalRelativeDifferences',differences, ...
    'newActualNSEvaluations',3,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
save(fullfile(folder,'riemann_selected_direction_unit_evidence.mat'),'report','QL','QR','normal','left','right','audit','sample','errors', ...
    'Q','cut','gas','boundary','scale','rate','info','direction','value','uncached','detail','cache','fd','fdEvidence','half','halfEvidence','-v7');
fid=fopen(fullfile(folder,'riemann_selected_direction_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
