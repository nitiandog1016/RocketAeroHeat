function report=testRocketTemperatureLookahead(folder)
%TESTROCKETTEMPERATURELOOKAHEAD 检查上下温域的前瞻选面、原行保持及显式拒绝。
% 模块一：制造均匀高速小域，预览当前方向在1/8处的冷热越域；不当作物理解。
% 模块二：原10K面、原H及正负号均保留，新增原态面与单元EOS导数独立比较。
% 模块三：以已保存Q35→Q36作为真实夹具，只复核选面，不沿用其方向推进Q37。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
s=load(fullfile(root,'outputs/v83_validation/riemannreservoirunit02/riemann_reservoir_unit_evidence.mat'),'result','cut','gas');
cut=s.cut;gas=s.gas;scale=s.result.residualScale;n=size(s.result.Q,1);cv=gas.R/(gas.gamma-1);
inside=find(cut.openFaces.area>0&cut.openFaces.cellMinus>0&cut.openFaces.cellPlus>0);cases=cell(3,1);
for mode=1:3
    rho=.1*ones(n,1);u=repmat([3000 5100 8200],n,1);T=100*ones(n,1);
    if mode==2,T(:)=1980;end
    if mode==3,T(:)=300;T(1)=100;T(end-1)=85;T(end)=1995;end
    Q=[rho rho.*u rho.*(cv*T+.5*sum(u.^2,2))];
    [~,audit]=cutCellExactEulerFlux(Q(cut.openFaces.cellMinus(inside),:),Q(cut.openFaces.cellPlus(inside),:),cut.openFaces.normal(inside,:),gas);audit.faceID=inside;
    pure=find(audit.sampleRegion==1|audit.sampleRegion==4);assert(~isempty(pure));
    if mode==3
        owner=double(cut.openFaces.cellMinus(inside(pure)));neighbour=double(cut.openFaces.cellPlus(inside(pure)));
        pureCV=owner;choose=audit.sampleRegion(pure)==4;pureCV(choose)=neighbour(choose);
        selected=find(pureCV==1,1);assert(~isempty(selected));slot=pure(selected);cellID=1;
    else
        slot=pure(1);cellID=double(cut.openFaces.cellMinus(inside(slot)));
        if audit.sampleRegion(slot)==4,cellID=double(cut.openFaces.cellPlus(inside(slot)));end
    end
    value=zeros(n,5);change=-200;if mode==2,change=200;end
    value(cellID,5)=rho(cellID)*cv*change/scale(5);info=struct('exactEulerFluxAudit',audit);
    base=selectRocketTemperatureBandRows(Q,info,cut,gas,scale);
    [rows,detail]=selectRocketTemperatureLookaheadRows(Q,info,cut,gas,scale,value,1/8);
    assert(detail.previewCellStateEvaluable&&isequal(detail.predictedViolatingCellIDs,cellID)&& ...
        detail.selectedTemperatureBandK>10&&ismember(inside(slot),rows.faceIDs)&&detail.originalRowsExactlyRetained);
    if mode==3,assert(~isempty(base.faceIDs));end
    density=value;density(:,1)=-16*rho/scale(1);
    [densityRows,densityReport]=selectRocketTemperatureLookaheadRows(Q,info,cut,gas,scale,density,1/8);
    assert(~densityReport.previewCellStateEvaluable&&isequaln(densityRows,base));
    zeroRows=selectRocketTemperatureLookaheadRows(Q,info,cut,gas,scale,zeros(n,5),1/8);assert(isequaln(zeroRows,base));
    cases{mode}=struct('Q',Q,'cut',cut,'gas',gas,'scale',scale,'value',value,'info',info,'base',base,'rows',rows,'detail',detail, ...
        'negativeDensityPreviewRetainsOriginalRows',true,'zeroDirectionRetainsOriginalRows',true,'manufacturedState',true);
end
source=fullfile(root,'outputs/v83_validation/krylov35/wall_velocity_krylov_evidence.mat');
next=fullfile(root,'outputs/v83_validation/krylov36/wall_velocity_krylov_evidence.mat');
geometry=fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat');
before=rocketStreamingFileSHA256(source);nextSHA=rocketStreamingFileSHA256(next);geometrySHA=rocketStreamingFileSHA256(geometry);
a=load(source,'result','gas');b=load(next,'result');d=load(geometry,'domain');trace=b.result.trace{1};
assert(isequal(a.result.Q,trace.initialQ));
[realRows,realReport]=selectRocketTemperatureLookaheadRows(a.result.Q,a.result.info,d.domain,a.gas,b.result.residualScale,trace.linearUpdateU,1/8);
assert(numel(realReport.predictedViolatingCellIDs)==3&&realReport.selectedTemperatureBandK>10);
assert(strcmp(before,rocketStreamingFileSHA256(source))&&strcmp(nextSHA,rocketStreamingFileSHA256(next)));
report=struct('status','passed-both-temperature-ends-original-row-retention-and-real-saved-Q35-selection-fixture', ...
    'sourceStateSHA256',before,'directionEndpointSHA256',nextSHA,'geometrySHA256',geometrySHA, ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'oldDirectionUsedToAdvanceNewState',false, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_lookahead_unit.mat'),'cases','realRows','realReport','report','-v7');disp(report);
end
