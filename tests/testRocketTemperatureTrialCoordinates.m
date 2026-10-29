function report=testRocketTemperatureTrialCoordinates(folder)
%TESTROCKETTEMPERATURETRIALCOORDINATES 验证实际候选路径和逐轮温域前瞻保持一致。
% 模块一：两条独立面记录检查同切向方向的非线性温度差，不冒充闭合NS体网格。
% 模块二：32CV人工面循环和耦合C/W制造系统检查两种排列的投影后前瞻；不是物理NS导数。
% 模块三：旧缺省输出另由兼容测试回放；原正性/界面温域/方向精度门槛不改。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
fixture=load(rocketTemperatureResearchFixture('interface_rejection.mat'),'baseRows');gas=fixture.baseRows.gas;
cv=gas.R/(gas.gamma-1);scale=ones(1,5);alpha=1/8;
Q=repmat([1 1000 0 0 cv*80.05+500000],2,1);value=repmat([0 2000 0 0 2000000],2,1);
cut=struct('openFaces',struct('area',ones(2,1),'cellMinus',[1;2],'cellPlus',[2;1], ...
    'normal',repmat([1 0 0],2,1)));
info=sample(Q,cut,gas);pairs=cell(2,1);names={'conservative','primitive-tangent'};
for mode=1:2
    coordinates=names{mode};candidate=rocketShellTrialState(Q,value,scale,alpha,gas,coordinates);
    [rows,detail]=selectRocketTemperatureCellInterfaceRows(Q,info,cut,gas,scale,value,alpha,coordinates);
    T=(candidate(:,5)-.5*sum(candidate(:,2:4).^2,2)./candidate(:,1))./(candidate(:,1)*cv);
    if mode==1
        assert(all(T<80)&&isequal(detail.predictedViolatingCellIDs,[1;2])&& ...
            isequal(detail.interfacePreview.blockedFaceIDs,[1;2]));
    else
        assert(max(abs(T-80.05))<1e-10&&isempty(detail.predictedViolatingCellIDs)&& ...
            detail.interfacePreview.allInteriorFacesAssessed&&isempty(detail.predictedViolatingInterfaceFaceIDs)&& ...
            strcmp(detail.trialStateCoordinates,'primitive-tangent'));
    end
    pairs{mode}=struct('Q',Q,'value',value,'scale',scale,'gas',gas,'cut',cut,'info',info, ...
        'coordinates',coordinates,'candidate',candidate,'rows',rows,'detail',detail);
end
% 保留已有耦合C/W，改变的只是明确标为制造的Q/右端和独立面记录。
p=load(rocketTemperatureResearchFixture('parent_least_squares_unit_evidence.mat'),'cases');
w=load(rocketTemperatureResearchFixture('peak_weighted_unit_evidence.mat'),'cases');cases=cell(2,1);
for j=1:2
    a=p.cases{j};b=w.cases{j};n=a.packed.cellCount;
    q=repmat([1 1000 0 0 cv*80.9+500000],n,1);
    f=struct('area',ones(n,1),'cellMinus',(1:n)','cellPlus',[2:n 1]', ...
        'normal',repmat([1 0 0],n,1));domain=struct('openFaces',f);state=sample(q,domain,gas);
    base=selectRocketTemperatureBandRows(q,state,domain,gas,scale);
    factor=buildRocketAugmentedPatchQR(a.packed,a.factor,b.metric,b.detail.penaltyInOriginalResidualUnits);
    A=[b.base.C;factor.penaltyInOriginalResidualUnits*b.W];desired=zeros(n,5);desired(1,5)=-cv*8;
    rhs=reshape((b.base.C')\((A'*A)*desired(:)),n,5);
    preview=struct('Q',q,'info',state,'cut',domain,'gas',gas,'scale',scale,'targetAlpha',alpha, ...
        'selectionPolicy','original-and-iteratively-predicted-cell-and-interface-faces', ...
        'trialStateCoordinates','primitive-tangent');
    [v,detail,evidence]=solveRocketTemperatureActiveDirection(a.packed,factor,rhs,b.metric,base, ...
        struct('maximumIterations',900,'lookaheadPreview',preview));
    assert(detail.activeSolveCount>=1&&detail.repreview.previewCount>=2&& ...
        strcmp(detail.repreview.trialStateCoordinates,'primitive-tangent'));
    values=cell(1,1+numel(evidence.cycles));values{1}=evidence.free;
    for k=1:numel(evidence.cycles),values{k+1}=evidence.cycles{k}.value;end
    for k=1:numel(values)
        selection=evidence.previewSelections{k};qc=rocketShellTrialState(q,values{k},scale,alpha,gas,'primitive-tangent');
        rho=qc(:,1);u=qc(:,2:4)./rho;temp=(qc(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);
        assert(isequal(find(temp<80|temp>2000),selection.predictedViolatingCellIDs)&& ...
            strcmp(selection.trialStateCoordinates,'primitive-tangent'));
    end
    cases{j}=struct('Q',q,'cut',domain,'gas',gas,'scale',scale,'info',state,'baseRows',base, ...
        'C',b.base.C,'W',b.W,'factor',factor,'metric',b.metric,'rhs',rhs,'value',v,'detail',detail,'evidence',evidence);
end
unknownRejected=false;
try,selectRocketTemperatureCellInterfaceRows(Q,info,cut,gas,scale,value,alpha,'unknown');
catch e,unknownRejected=strcmp(e.identifier,'TemperatureLookahead:Coordinates');end
assert(unknownRejected);
% 模块四：必须使用真正由选择器返回的空集合，不用手工清空已有行替代。
emptyQ=repmat([1 1000 0 0 cv*100+500000],2,1);emptyInfo=sample(emptyQ,cut,gas);
actualEmpty=selectRocketTemperatureBandRows(emptyQ,emptyInfo,cut,gas,scale);
wider=selectRocketTemperatureBandRows(emptyQ,emptyInfo,cut,gas,scale,struct('temperatureBandK',30));
assert(isempty(actualEmpty.faceIDs)&&isequal(actualEmpty.scale,scale)&&isequal(actualEmpty.gas,gas));
mergedEmpty=mergeRocketTemperatureRequiredRows(actualEmpty,wider);
assert(isequal(mergeRocketTemperatureRequiredRows(wider,actualEmpty),wider)&& ...
    isequal(mergedEmpty.H,wider.H)&& ...
    isequal(mergeRocketTemperatureRequiredRows(actualEmpty,actualEmpty),actualEmpty));
emptyRejects=false(3,2);
for field=1:3
    wrong=actualEmpty;
    if field==1,wrong.gas.R=wrong.gas.R*1.01;
    elseif field==2,wrong.scale(1)=wrong.scale(1)*1.01;
    else,wrong=rmfield(wrong,'gas');end
    for order=1:2
        try
            if order==1,mergeRocketTemperatureRequiredRows(wrong,wider);else,mergeRocketTemperatureRequiredRows(wider,wrong);end
        catch e,emptyRejects(field,order)=strcmp(e.identifier,'TemperatureRepreview:Identity');end
    end
end
assert(all(emptyRejects,'all'));
report=struct('status','passed-coordinate-consistent-pair-preview-and-two-ordering-projected-rechecks', ...
    'manufacturedFaceRecordsNotClosedPhysicalMesh',true,'manufacturedCNotPhysicalNSJacobian',true, ...
    'unknownCoordinatesRejected',unknownRejected,'actualEmptySelectorMergeAndIdentityPassed',true, ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_trial_coordinates_unit.mat'),'pairs','cases','report','actualEmpty','wider','emptyRejects','-v7');disp(report);
end
function info=sample(Q,cut,gas)
f=cut.openFaces;ids=find(f.area>0&f.cellMinus>0&f.cellPlus>0);
[~,audit]=cutCellExactEulerFlux(Q(f.cellMinus(ids),:),Q(f.cellPlus(ids),:),f.normal(ids,:),gas);
audit.faceID=ids;info=struct('exactEulerFluxAudit',audit);
end
