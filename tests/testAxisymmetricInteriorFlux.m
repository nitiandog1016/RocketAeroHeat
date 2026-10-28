function report=testAxisymmetricInteriorFlux(folder,sourceFolder,trialFile)
%TESTAXISYMMETRICINTERIORFLUX 内部HLLE/原HLLC与实际NS+SA的质量一致耦合校核。
% 模块一：只改内部Euler通量，边界/黏性/几何源保持，SA使用同一实际质量通量。
% 模块二：四五量完整有色/独立逐列Jacobian、真实BE及全部共享面/源积分账。
% 模块三：逐值回放旧v72默认11端点，不重绑旧源；两种HLLE重算253真实候选。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v73_validation','interior_flux');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v72_validation','constant_uniform');end
if nargin<3,trialFile=fullfile(root,'outputs','v72_validation','continuation_trial','continuation_evidence.mat');end
assert(~exist(folder,'dir'),'InteriorFluxTest:Overwrite','不能覆盖真实证据。');
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[200 800], ...
    'transport',struct('model','sutherland','referenceViscosity',.05,'referenceTemperature',300,'sutherlandTemperature',110.4)));
[x,r]=ndgrid([0 .02 .05 .08],linspace(.03,.06,5));
mesh=buildAxisymmetricFVGeometry(x+.07*(r-.03),r,struct('periodicI',true));nc=mesh.cellCount;t=(1:nc)';
P=[1+.02*sin(t),-2+.2*cos(t),.02*sin(.8*t),500+60*cos(.4*t)];
Q4=axisymmetricFVConservative(P,gas);Q5=[Q4,P(:,1).*(.3+.02*sin(.6*t))];
per=struct('type','periodic');baseBC=struct('iMin',per,'iMax',per, ...
    'jMin',struct('type','wall','temperature',350,'velocity',[-2 0]), ...
    'jMax',struct('type','wall','temperature',650,'velocity',[0 0]),'pressureSourceQuadrature','cell-constant');
scale5=[1 20 20 2000 .5];records=struct([]);cases={};
for nq=[4 5]
    if nq==4,Q=Q4;else,Q=Q5;end
    scale=scale5(1:nq);U=Q./scale;
    [oldRate,oldInfo]=evaluate(Q,baseBC);hllc=baseBC;hllc.interiorConvectiveFlux='hllc';
    [same,sameInfo]=evaluate(Q,hllc);
    assert(isequal(same,oldRate)&&isequal(sameInfo.faceFluxIntegralPerRadian,oldInfo.faceFluxIntegralPerRadian)&& ...
        isequal(sameInfo.convectiveFrequency,oldInfo.convectiveFrequency));
    for policy={'hlle-davis','hlle-einfeldt'}
        name=policy{1};bc=baseBC;bc.interiorConvectiveFlux=name;
        [rate,info]=evaluate(Q,bc);audit=info.interiorConvectiveFlux;
        assert(contains(info.method,'first-order HLLE')&&~contains(info.method,'first-order HLLC'));
        ids=audit.faceID;bf=mesh.boundaryFaces;area=mesh.faces.areaPerRadian;
        assert(isequal(info.cellSourceIntegral,oldInfo.cellSourceIntegral)&&isequaln(info.wall,oldInfo.wall)&& ...
            isequal(info.faceFluxIntegralPerRadian(bf,:),oldInfo.faceFluxIntegralPerRadian(bf,:)));
        deltaFlux=info.faceFluxIntegralPerRadian-oldInfo.faceFluxIntegralPerRadian;
        expected=-mesh.incidence*deltaFlux./mesh.volumePerRadian;
        budget=2048*eps*max(1,max(abs(rate)+abs(oldRate),[],'all'));
        assert(max(abs(rate-oldRate-expected),[],'all')<=budget&&max(abs(rate-oldRate),[],'all')>1e-3);
        o=mesh.faces.owner(ids);j=mesh.faces.neighbour(ids);normal=mesh.faces.normal(ids,:);
        primitive=info.primitive;
        oldSpeed=max(abs(sum(primitive(o,2:3).*normal,2))+primitive(o,6), ...
            abs(sum(primitive(j,2:3).*normal,2))+primitive(j,6));
        speedIncrement=zeros(numel(area),1);speedIncrement(ids)=max(oldSpeed,max(abs(audit.waveSpeedsMPerS),[],2))-oldSpeed;
        predictedFrequency=oldInfo.convectiveFrequency+(mesh.cellFaceMultiplicity*(speedIncrement.*area))./mesh.volumePerRadian;
        assert(all(abs(info.convectiveFrequency-predictedFrequency)<=128*eps*(abs(info.convectiveFrequency)+abs(predictedFrequency))));
        if nq==5
            sa=info.sa;o=mesh.faces.owner(ids);j=mesh.faces.neighbour(ids);N=Q(:,5)./Q(:,1);
            mass=audit.fluxPerArea(:,1).*area(ids);upwind=N(o);upwind(mass<0)=N(j(mass<0));
            % (M*N)*A与(M*A)*N是同一三因子积；两条各两次乘法的差按舍入预算核对。
            expectedScalar=mass.*upwind;
            assert(isequal(sa.faceMassFluxPerRadian(ids),mass)&& ...
                all(abs(sa.faceAdvectiveFluxPerRadian(ids)-expectedScalar)<=8*eps*(abs(expectedScalar)+abs(sa.faceAdvectiveFluxPerRadian(ids))))&& ...
                isequal(sa.faceDiffusiveFluxPerRadian,oldInfo.sa.faceDiffusiveFluxPerRadian));
        end
        ledger=sum(rate.*mesh.volumePerRadian,1)+info.netBoundaryFlux-info.sourceIntegral;
        assert(max(abs(ledger))<1e-10);
        [J,details,pattern]=axisymmetricFVJacobian(Q,mesh,gas,bc,scale,struct('differenceStep',2e-6));assert(details.success);
        dense=zeros(nq*nc);
        for column=1:nq*nc
            h=2e-6*max(1,abs(U(column)));plus=U;minus=U;plus(column)=plus(column)+h;minus(column)=minus(column)-h;
            rp=evaluate(plus.*scale,bc)./scale;rm=evaluate(minus.*scale,bc)./scale;
            dense(:,column)=(rp(:)-rm(:))/(2*h);
        end
        relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');assert(relative<3e-7);
        outside=~full(kron(ones(nq),pattern.cellPattern));assert(max([0;abs(dense(outside))])<1e-6);
        options=struct('maxSteps',1,'maximumPseudoStep',1e-8,'residualScale',scale,'steadyResidualTolerance',0,'jacobianRefreshEvery',1);
        be=advanceAxisymmetricSteady(Q,mesh,gas,bc,options);assert(be.acceptedSteps==1&&~be.steadyConverged);
        [endRate,endInfo]=evaluate(be.Q,bc);G=be.Q./scale-U-be.history.dt(end)*(endRate./scale);
        assert(isequaln(endInfo,be.finalInfo)&&max(abs(G),[],'all')<=be.history.newtonTarget(end)+1024*eps&&be.globalBalanceRelative<5e-10);
        row=struct('quantities',nq,'policy',name,'jacobianColoredVsSerialRelative',relative, ...
            'conservationLedger',ledger,'residualCalls',be.residualEvaluations,'BEBalanceRelative',be.globalBalanceRelative);
        if isempty(records),records=row;else,records(end+1)=row;end %#ok<AGROW>
        cases{end+1}=struct('Q',Q,'mesh',mesh,'boundary',bc,'rate',rate,'info',info,'oldInfo',oldInfo,'J',J,'serial',dense,'be',be); %#ok<AGROW>
    end
end
refusals=struct('label',{},'identifier',{});
bad=baseBC;bad.interiorConvectiveFlux='automatic-clip';reject('unknown-policy',@()evaluate(Q5,bad));
bad=baseBC;bad.jMin.interiorConvectiveFlux='hlle-einfeldt';reject('silently-ignored-side-setting',@()evaluate(Q5,bad));

%% 三、旧完整默认路径数值兼容；新源不得冒充旧链执行身份
saved=load(fullfile(sourceFolder,'rocket_axisymmetric_evidence.mat'),'caseDefinition','durable','implicit');c=saved.caseDefinition;
assert(saved.implicit.acceptedSteps==11&&~isfield(c.boundary,'interiorConvectiveFlux'));
for k=1:11
    s=load(fullfile(sourceFolder,saved.durable.steps(k).relativePath),'segment');s=s.segment;
    [rate,info]=axisymmetricRANSResidual(s.Q,c.mesh,c.gas,c.boundary);
    assert(isequaln(info,s.finalInfo)&&max(abs(rate)./s.residualScale,[],'all')==s.residualNorm);
end
savedTrial=load(trialFile,'c','diagnostics');candidateRecords=struct([]);
for policy={'hlle-davis','hlle-einfeldt'}
    name=policy{1};bc=savedTrial.c.boundary;bc.interiorConvectiveFlux=name;count=0;minimumRho=Inf;minimumInternal=Inf;
    for k=1:numel(savedTrial.diagnostics.records)
        q=savedTrial.diagnostics.records(k).candidateQ;
        [rate,info]=axisymmetricRANSResidual(q,savedTrial.c.mesh,savedTrial.c.gas,bc);
        assert(all(isfinite(rate),'all')&&strcmp(info.interiorConvectiveFlux.policy,name));
        a=info.interiorConvectiveFlux;
        minimumRho=min(minimumRho,min(a.intermediateConservativeState(a.middleMask,1)));
        minimumInternal=min(minimumInternal,min(a.intermediateInternalEnergyPerVolume(a.middleMask)));count=count+1;
    end
    row=struct('policy',name,'completeRANSCandidatesAccepted',count,'minimumUsedIntermediateDensity',minimumRho, ...
        'minimumUsedIntermediateInternalEnergy',minimumInternal);
    if isempty(candidateRecords),candidateRecords=row;else,candidateRecords(end+1)=row;end %#ok<AGROW>
    assert(count==253);
end
report=struct('status','passed-explicit-interior-flux-full-NS-SA-Jacobian-BE-and-real-candidate-contract', ...
    'matlabEngineProcessID',feature('getpid'),'coupledCases',records,'legacyElevenInfoReplayed',true, ...
    'realCandidateRuns',candidateRecords,'refusals',refusals,'sharedNSandSAMassFluxPassed',true, ...
    'sourceOrBoundaryFluxChangedAtSameQ',false,'oldChainRebound',false, ...
    'physicalSolution',false,'experimentalValidation',false,'generalPositivityProved',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'interior_flux_evidence.mat'),'report','gas','cases','-v7');
fid=fopen(fullfile(folder,'interior_flux_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function [rate,info]=evaluate(state,boundary)
        if size(state,2)==5,[rate,info]=axisymmetricRANSResidual(state,mesh,gas,boundary);
        else,[rate,info]=axisymmetricNavierStokesResidual(state,mesh,gas,boundary);end
    end
    function reject(label,action)
        try,action();catch exception
            assert(strcmp(exception.identifier,'AxisFV:InteriorFluxPolicy'),'InteriorFluxTest:WrongRefusal','非预期%s。',exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('InteriorFluxTest:NoRefusal','非法输入%s没有被拒绝。',label);
    end
end
