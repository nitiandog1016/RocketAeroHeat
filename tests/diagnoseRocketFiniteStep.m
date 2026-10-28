function report=diagnoseRocketFiniteStep(folder)
%DIAGNOSEROCKETFINITESTEP 对冻结Q22到Q24作逐面对流非线性分解，生成真实几何小复现。
% 全域分解只读历史，不重新授予方向资格；小域保留两圈状态以精确复现内部NS行。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);timer=tic;
stateFile=fullfile(root,'outputs/v83_validation/krylov24/wall_velocity_krylov_evidence.mat');
geometryFile=fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat');
s=load(stateFile,'result','gas','boundary');d=load(geometryFile,'domain');cut=d.domain;r=s.result;gas=s.gas;
Q0=r.initialQ;Q1=r.Q;dQ=Q1-Q0;scale=r.residualScale;f=cut.openFaces;n=size(Q0,1);
ids=find(f.cellMinus>0&f.cellPlus>0&f.area>0);m=f.cellMinus(ids);p=f.cellPlus(ids);
F0=r.initialInfo.convectiveFaceFlux(ids,:);F1=r.info.convectiveFaceFlux(ids,:);
% 分块限制临时Jacobian内存；导数作用采用实际已接受增量，单位不变。
dF=zeros(numel(ids),5);
for first=1:25000:numel(ids)
    k=first:min(first+24999,numel(ids));[L,R]=cutCellExactEulerBranchJacobian(Q0(m(k),:),Q0(p(k),:),f.normal(ids(k),:),gas);
    dF(k,:)=(reshape(pagemtimes(L,reshape(dQ(m(k),:)',5,1,[]))+pagemtimes(R,reshape(dQ(p(k),:)',5,1,[])),5,[]))';
end
remainder=F1-F0-dF;changed=r.initialInfo.exactEulerFluxAudit.sampleRegion~=r.info.exactEulerFluxAudit.sampleRegion;
div=@(F) divergence(F,ids,cut,n)./scale;
conv=div(remainder);same=div(remainder.*~changed);switching=div(remainder.*changed);
total=(r.rate-r.initialRate)./scale-r.trace{1}.alpha*r.linearLog{1}.steadyDirectionEvidence.value;
% 全域对流之外的剩余包括完整黏性及边界非线性和原NS差分舍入。
report=struct('status','diagnosed-saved-Q22-to-Q24','changedInteriorFaces',nnz(changed), ...
    'convectiveRemainderNorm',norm(conv,'fro'),'unchangedRegionRemainderNorm',norm(same,'fro'), ...
    'changedRegionRemainderNorm',norm(switching,'fro'),'totalRemainderNorm',norm(total,'fro'), ...
    'outsideInteriorConvectionNorm',norm(total-conv,'fro'), ...
    'unchangedChangedDot',sum(same.*switching,'all'),'newProductionSolvedAreaM2',0);
energy=sum(total.^2,2);[~,order]=sort(energy,'descend');seeds=order(1:3);facets=mod(seeds-1,19752)+1;
% 只取第3至6层为活动核心，避免把人工截断边界混同真实远场/墙条件。
core=[];
for facet=facets(:)'
    c=facet+3*19752;adj=unique([m(p==c);p(m==c);c]);adj=adj(floor((adj-1)/19752)==3);
    ff=mod(adj-1,19752)+1;core=[core;reshape(ff+(2:5)*19752,[],1)]; %#ok<AGROW>
end
core=unique(core);one=unique([core;m(ismember(p,core));p(ismember(m,core))]);
two=unique([one;m(ismember(p,one));p(ismember(m,one))]);map=zeros(n,1);map(two)=1:numel(two);
fluxFaces=find(ismember(f.cellMinus,core)|ismember(f.cellPlus,core));
gradientFaces=find(ismember(f.cellMinus,one)|ismember(f.cellPlus,one));
assert(all(f.cellMinus(gradientFaces)>0&f.cellPlus(gradientFaces)>0),'Patch:Boundary','此内部复现不能触及外域或墙。');
assert(~any(ismember(cut.walls.ownerCell,one)),'Patch:Wall','梯度模板不能漏墙。');
patch=struct('globalCellID',two,'core',map(core),'gradientCells',map(one), ...
    'centroid',cut.cells.centroid(two,:),'volume',cut.cells.volume(two), ...
    'fluxGlobalFaceID',fluxFaces,'gradientGlobalFaceID',gradientFaces, ...
    'fluxMinus',map(f.cellMinus(fluxFaces)),'fluxPlus',map(f.cellPlus(fluxFaces)), ...
    'fluxArea',f.area(fluxFaces),'fluxNormal',f.normal(fluxFaces,:), ...
    'gradientMinus',map(f.cellMinus(gradientFaces)),'gradientPlus',map(f.cellPlus(gradientFaces)), ...
    'gradientArea',f.area(gradientFaces),'Q0',Q0(two,:),'Q1',Q1(two,:), ...
    'direction',r.trace{1}.linearUpdateU(two,:),'scale',scale,'gas',gas, ...
    'expectedRate0',r.initialRate(core,:),'expectedRate1',r.rate(core,:), ...
    'expectedGradient0',r.initialInfo.gradient(:,:,one),'expectedGradient1',r.info.gradient(:,:,one), ...
    'expectedFlux0',r.initialInfo.faceFlux(fluxFaces,:),'expectedFlux1',r.info.faceFlux(fluxFaces,:), ...
    'sourceStateSHA256',rocketStreamingFileSHA256(stateFile),'sourceGeometrySHA256',rocketStreamingFileSHA256(geometryFile), ...
    'scope','exact-interior-NS-row-replay-with-two-hop-halo; frozen-exterior-only-for-local-experiments');
save(fullfile(folder,'finite_step_patch.mat'),'patch','-v7');
% 最坏核心的全部真实面，按同一方向扫描，分辨分支切换前的曲率。
worst=order(1);select=find(m==worst|p==worst);sweep=struct([]);
for fraction=[0 2.^(-12:-1) 1]
    qa=Q0+fraction*dQ;[fa,a]=cutCellExactEulerFlux(qa(m(select),:),qa(p(select),:),f.normal(ids(select),:),gas);
    row=struct('fractionOfAcceptedStep',fraction,'faceID',ids(select),'flux',fa,'region',a.sampleRegion, ...
        'pressureStar',a.pressureStarPa,'normalVelocityStar',a.normalVelocityStarMs, ...
        'nonlinearFluxRemainder',fa-F0(select,:)-fraction*dF(select,:));
    if isempty(sweep),sweep=row;else,sweep(end+1)=row;end %#ok<AGROW>
end
report.patchCells=numel(two);report.patchCoreCells=numel(core);report.patchFluxFaces=numel(fluxFaces);
report.elapsedSeconds=toc(timer);report.worstCell=worst;
save(fullfile(folder,'finite_step_decomposition.mat'),'report','sweep','conv','same','switching','total','-v7');
fid=fopen(fullfile(folder,'diagnosis.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
function rate=divergence(F,ids,cut,n)
% 共享面按原体积成对装配，不把各面范数误当单元残差范数。
f=cut.openFaces;m=f.cellMinus(ids);p=f.cellPlus(ids);rate=zeros(n,5);
for c=1:5
    value=f.area(ids).*F(:,c);rate(:,c)=(-accumarray(m,value,[n 1])+accumarray(p,value,[n 1]))./cut.cells.volume;
end
end
