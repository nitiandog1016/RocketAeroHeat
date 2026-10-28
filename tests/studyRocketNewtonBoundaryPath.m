function report=studyRocketNewtonBoundaryPath(folder)
%STUDYROCKETNEWTONBOUNDARYPATH 定位真实已保存Newton候选的开放特征边界限制。
% 模块一：只从krylov19原始方向重建实际候选，逐面调用原特征边界，不重新求NS。
% 模块二：保存原面/CV/层/位置、法向Mach和声学重构裕量，区分状态与分支拒绝。
% 模块三：对原分支拒绝面另试内部精确Riemann储库，仅为候选边界方法诊断。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'BoundaryPath:Overwrite','请使用未存在的诊断目录。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','krylov19','wall_velocity_krylov_evidence.mat');
s=load(sourceFile,'result','gas','boundary');r=s.result;gas=s.gas;boundary=s.boundary;
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;f=cut.openFaces;
Q=r.initialQ;direction=r.linearLog{1}.steadyDirectionEvidence.directionU;scale=r.residualScale;
assert(isequal(r.trace{1}.linearUpdateU,direction)&&isequal(r.trace{1}.initialQ,Q));
ids=find(f.area>0&(f.cellMinus==0|f.cellPlus==0));count=numel(ids);owner=f.cellMinus(ids);normal=f.normal(ids,:);side=2*ones(count,1);
left=owner==0;owner(left)=f.cellPlus(ids(left));normal(left,:)=-normal(left,:);side(left)=1;
names={'xmin','xmax';'ymin','ymax';'zmin','zmax'};labels=cell(count,1);reservoir=zeros(count,5);
for j=1:count
    labels{j}=names{f.axis(ids(j)),side(j)};b=boundary.(labels{j});
    assert(strcmp(b.type,'characteristic')&&strcmp(b.mode,'riemannFarfield')&&~isa(b.reservoirPrimitive,'function_handle'));
    reservoir(j,:)=b.reservoirPrimitive;
end
sourceFaces=size(cut.inputGeometry.F,1);map=struct('openFaceID',ids,'ownerCellID',owner,'normal',normal,'boundaryName',{labels}, ...
    'centroidM',f.centroid(ids,:),'sourceFacetID',mod(owner-1,sourceFaces)+1,'layer',floor((owner-1)/sourceFaces)+1);
cv=gas.R/(gas.gamma-1);trials=struct([]);lastAccepted=0;
for index=0:numel(r.lineLog)
    alpha=0;accepted=false;originalError='';if index>0,a=r.lineLog(index);alpha=a.alpha;accepted=a.accepted;originalError=a.errorIdentifier;end
    q=Q+alpha*direction.*scale;rho=q(:,1);energy=q(:,5)-.5*sum(q(:,2:4).^2,2)./rho;T=energy./(rho*cv);
    valid=all(rho>0)&&all(energy>0)&&all(T>=80&T<=2000);if ~valid,continue;end
    P=[rho,q(:,2:4)./rho,T];primitive=P(owner,:);facePrimitive=nan(count,5);errorID=cell(count,1);errorMessage=cell(count,1);branches=cell(count,1);
    cellNormal=sum(primitive(:,2:4).*normal,2);sound=sqrt(gas.gamma*gas.R*primitive(:,5));
    reservoirNormal=sum(reservoir(:,2:4).*normal,2);reservoirSound=sqrt(gas.gamma*gas.R*reservoir(:,5));
    Jplus=cellNormal+2*sound/(gas.gamma-1);Jminus=reservoirNormal-2*reservoirSound/(gas.gamma-1);
    acousticVelocity=.5*(Jplus+Jminus);acousticSound=(gas.gamma-1)*.25*(Jplus-Jminus);
    failed=false(count,1);
    for j=1:count
        try
            [facePrimitive(j,:),tag]=resolveNSCharacteristicBoundary3D(primitive(j,:),boundary.(labels{j}),normal(j,:),gas);branches{j}=tag.branch;
        catch ex
            assert(startsWith(ex.identifier,'RocketAeroHeat:NSDraft:Characteristic'));failed(j)=true;errorID{j}=ex.identifier;errorMessage{j}=ex.message;
        end
    end
    first=find(failed,1);originalGuardReproduced=true;
    if index==0||accepted,assert(~any(failed));
    elseif startsWith(originalError,'RocketAeroHeat:NSDraft:Characteristic'),originalGuardReproduced=strcmp(errorID{first},originalError);assert(originalGuardReproduced);end
    exact=struct('attemptedFaces',0,'succeeded',false,'failureIdentifier','','samplePrimitive',zeros(0,5),'sampleRegion',zeros(0,1), ...
        'vacuumVelocityMarginMs',zeros(0,1),'flux',zeros(0,5),'physicalBoundaryChanged',false);
    if any(failed)
        exact.attemptedFaces=nnz(failed);qr=reservoir(failed,:);ql=primitive(failed,:);nrm=normal(failed,:);
        L=conservative(ql);R=conservative(qr);
        try
            [exact.flux,audit]=cutCellExactEulerFlux(L,R,nrm,gas);exact.succeeded=true;
            for key={'samplePrimitive','sampleRegion','vacuumVelocityMarginMs'},exact.(key{1})=audit.(key{1});end
        catch ex,exact.failureIdentifier=ex.identifier;end
    end
    item=struct('sourceTrial',index,'alpha',alpha,'acceptedByOriginalSolver',accepted,'originalErrorIdentifier',originalError, ...
        'primitive',primitive,'facePrimitive',facePrimitive,'branch',{branches},'failed',failed,'errorIdentifier',{errorID},'errorMessage',{errorMessage}, ...
        'cellNormalMach',cellNormal./sound,'reservoirNormalMach',reservoirNormal./reservoirSound, ...
        'acousticNormalMach',acousticVelocity./acousticSound,'acousticSubsonicMarginMs',acousticSound-abs(acousticVelocity), ...
        'originalGuardReproduced',originalGuardReproduced,'firstFailedIndex',first,'exactReservoirDiagnostic',exact);
    if index>0&&accepted,lastAccepted=numel(trials)+1;end
    if isempty(trials),trials=item;else,trials(end+1)=item;end %#ok<AGROW>
end
assert(lastAccepted>0&&~any(trials(lastAccepted).failed));
details=struct([]);
for j=1:numel(trials)
    t=trials(j);i=t.firstFailedIndex;location=struct;
    if ~isempty(i),location=struct('openFaceID',ids(i),'ownerCellID',owner(i),'sourceFacetID',map.sourceFacetID(i),'layer',map.layer(i), ...
        'boundaryName',labels{i},'centroidM',map.centroidM(i,:),'cellNormalMach',t.cellNormalMach(i), ...
        'acousticNormalMach',t.acousticNormalMach(i),'acousticSubsonicMarginMs',t.acousticSubsonicMarginMs(i),'errorIdentifier',t.errorIdentifier{i});end
    row=struct('alpha',t.alpha,'sourceTrial',t.sourceTrial,'failedBoundaryFaces',nnz(t.failed),'firstFailure',location, ...
        'acceptedByOriginalSolver',t.acceptedByOriginalSolver,'exactDiagnosticAttemptedFaces',t.exactReservoirDiagnostic.attemptedFaces, ...
        'exactDiagnosticSucceeded',t.exactReservoirDiagnostic.succeeded,'exactDiagnosticFailure',t.exactReservoirDiagnostic.failureIdentifier);
    if isempty(details),details=row;else,details(end+1)=row;end %#ok<AGROW>
end
report=struct('status','tested-original-characteristic-boundary-path-with-separate-exact-reservoir-diagnostic', ...
    'sourceFile',sourceFile,'sourceSHA256',rocketStreamingFileSHA256(sourceFile),'boundaryFaces',count,'trialDetails',details, ...
    'originalBoundaryGuardKept',true,'allOriginalCellLegalCandidatesChecked',true,'newActualNSEvaluations',0, ...
    'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'newton_boundary_path_evidence.mat'),'Q','direction','scale','gas','boundary','map','reservoir','trials','report','-v7');
fid=fopen(fullfile(folder,'newton_boundary_path_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function q=conservative(p)
        q=[p(:,1),p(:,1).*p(:,2:4),p(:,1).*(cv*p(:,5)+.5*sum(p(:,2:4).^2,2))];
    end
end
