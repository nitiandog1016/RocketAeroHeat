function report=testRocketCellInterfaceLookahead(folder)
%TESTROCKETCELLINTERFACELOOKAHEAD 同时验证单元与星区界面的真实温域前瞻，不运行新NS。
% 模块一：真实Q37拒绝面仅作两CV夹具，并反转两侧/法向检查另一星区。
% 模块二：保留接触波零速度案例；另加整体平移案例供同波系公开差分验证。
% 模块三：旧三制造域核对单元阻挡面完整列出、原行不变，非温域错误必须传播。
% 模块四：多面、非连续面号和不同批量验证全部越温面及阻挡面均不遗漏。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
s=load(rocketTemperatureResearchFixture('interface_rejection.mat','outputs/v90_validation/interface01/interface_rejection.mat'),'baseRows','rejectedInputs');
gas=s.baseRows.gas;scale=s.baseRows.scale;cv=gas.R/(gas.gamma-1);pairs=cell(4,1);
for mode=1:4
    Q=[s.baseRows.QL;s.baseRows.QR];qc=[s.rejectedInputs.QL;s.rejectedInputs.QR];normal=s.baseRows.normal;
    if mode==2,Q=flipud(Q);qc=flipud(qc);normal=-normal;end
    if mode>=3
        normal=[1 0 0];rho=.1*ones(2,1);u=[-30 0 0;30 0 0];v=[-150 0 0;150 0 0];
        if mode==4,u(:,1)=u(:,1)+20;v(:,1)=v(:,1)+20;end
        Q=[rho rho.*u rho.*(cv*100+.5*sum(u.^2,2))];qc=[rho rho.*v rho.*(cv*100+.5*sum(v.^2,2))];
    end
    cut=struct('openFaces',struct('area',1,'cellMinus',1,'cellPlus',2,'normal',normal));
    [~,audit]=cutCellExactEulerFlux(Q(1,:),Q(2,:),normal,gas);audit.faceID=1;info=struct('exactEulerFluxAudit',audit);
    value=(qc-Q)./scale*8;candidate=rocketShellTrialState(Q,value,scale,1/8,gas,'conservative');
    base=selectRocketTemperatureBandRows(Q,info,cut,gas,scale);assert(isempty(base.faceIDs));
    [old,oldReport]=selectRocketTemperatureLookaheadRows(Q,info,cut,gas,scale,value,1/8,struct('selectionPolicy','original-and-predicted-required-faces'));
    assert(isempty(old.faceIDs)&&isempty(oldReport.predictedViolatingCellIDs));
    [rows,detail]=selectRocketTemperatureLookaheadRows(Q,info,cut,gas,scale,value,1/8, ...
        struct('selectionPolicy','original-and-predicted-cell-and-interface-faces'));
    p=detail.interfacePreview;one=findRocketCandidateTemperatureFaces(candidate,cut,gas,struct('maximumBatchFaces',1));
    assert(isequal(p.predictedViolatingInterfaceFaceIDs,1)&&p.allInteriorFacesAssessed&&p.allInputFacesAccountedFor&& ...
        isempty(p.blockedFaceIDs)&&isequal(rows.faceIDs,1)&&isequal(one.predictedViolatingInterfaceFaceIDs,p.predictedViolatingInterfaceFaceIDs));
    assert(abs(one.predictedInterfaceTemperatureK-p.predictedInterfaceTemperatureK)<=1e-10&&~p.temperatureDomainRelaxed);
    assert(ismember(rows.sampleRegion,[2 5])&&detail.allPredictedViolatingInterfacesRepresented);
    if mode<=2,assert(rows.temperatureBoundSigns==1&&p.predictedInterfaceTemperatureK>2000);
    else,assert(rows.temperatureBoundSigns==-1&&p.predictedInterfaceTemperatureK<80);end
    pairs{mode}=struct('Q',Q,'candidate',candidate,'cut',cut,'gas',gas,'scale',scale,'value',value,'info',info,'oldRows',old,'rows',rows,'detail',detail);
end

% 每个面使用独立单元对，加入面积为零及外边界槽，禁止把局部批次面号当成全局面号。
items=pairs([1 4 1 4 3 4]);candidate=vertcat(items{1}.candidate,items{2}.candidate, ...
    items{3}.Q,items{4}.candidate,items{5}.candidate,items{6}.Q);
candidate(7,5)=candidate(7,1)*(cv*2100+.5*sum((candidate(7,2:4)/candidate(7,1)).^2));
faceIDs=[1;3;4;6;7;8];minus=zeros(8,1);plus=zeros(8,1);area=ones(8,1);area(2)=0;
normal=repmat([1 0 0],8,1);
for k=1:6
    minus(faceIDs(k))=2*k-1;plus(faceIDs(k))=2*k;normal(faceIDs(k),:)=items{k}.cut.openFaces.normal;
end
cut=struct('openFaces',struct('area',area,'cellMinus',minus,'cellPlus',plus,'normal',normal));
batchReports=cell(3,1);batchSizes=[1 3 4096];
for k=1:numel(batchSizes)
    batchReports{k}=findRocketCandidateTemperatureFaces(candidate,cut,gas,struct('maximumBatchFaces',batchSizes(k)));
    assert(isequal(batchReports{k}.predictedViolatingInterfaceFaceIDs,[1;3;7])&& ...
        isequal(batchReports{k}.blockedFaceIDs,6)&&~batchReports{k}.allInteriorFacesAssessed&& ...
        batchReports{k}.allInputFacesAccountedFor);
end
batchCase=struct('candidate',candidate,'cut',cut,'gas',gas,'reports',{batchReports});
a=load(rocketTemperatureResearchFixture('temperature_lookahead_unit.mat','outputs/v90_validation/unit01/temperature_lookahead_unit.mat'),'cases');domains=cell(numel(a.cases),1);
for j=1:numel(a.cases)
    c=a.cases{j};
    [old,oldDetail]=selectRocketTemperatureLookaheadRows(c.Q,c.info,c.cut,c.gas,c.scale,c.value,1/8,struct('selectionPolicy','original-and-predicted-required-faces'));
    [rows,detail]=selectRocketTemperatureCellInterfaceRows(c.Q,c.info,c.cut,c.gas,c.scale,c.value,1/8);
    [found,slots]=ismember(old.faceIDs,rows.faceIDs);assert(all(found)&&isequal(old.H,rows.H(slots,:))&& ...
        isequal(old.temperatureBoundSigns,rows.temperatureBoundSigns(slots)));
    p=detail.interfacePreview;assert(~p.allInteriorFacesAssessed&&p.allInputFacesAccountedFor&& ...
        ~isempty(p.blockedFaceIDs)&&isequal(sort([p.assessedFaceIDs;p.blockedFaceIDs]),c.info.exactEulerFluxAudit.faceID));
    candidate=rocketShellTrialState(c.Q,c.value,c.scale,1/8,c.gas,'conservative');
    domains{j}=struct('source',c,'candidate',candidate,'oldRows',old,'oldDetail',oldDetail,'rows',rows,'detail',detail);
end
c=pairs{3};rho=c.Q(:,1);v=[-2000 0 0;2000 0 0];qc=[rho rho.*v rho.*(cv*100+.5*sum(v.^2,2))];value=(qc-c.Q)./scale*8;
vacuumRejected=false;identifier='';try
    selectRocketTemperatureCellInterfaceRows(c.Q,c.info,c.cut,c.gas,scale,value,1/8);
catch e,identifier=e.identifier;vacuumRejected=strcmp(identifier,'CutExactEuler:Vacuum');end
assert(vacuumRejected,'InterfacePreview:Suppression','非温域原内核错误不能当作已检查或静默忽略。');
report=struct('status','passed-upper-lower-star-state-preview-explicit-cell-blockage-original-row-retention-and-vacuum-propagation', ...
    'originalKernelTemperatureRangeUnchanged',true,'nonTemperatureFailurePropagated',true, ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'cell_interface_lookahead_unit.mat'),'pairs','domains','batchCase','report','identifier','-v7');disp(report);
end
