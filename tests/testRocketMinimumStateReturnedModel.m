function report=testRocketMinimumStateReturnedModel(folder)
%TESTROCKETMINIMUMSTATERETURNEDMODEL 用已保存Q24核对可行/不可行时返回方向与报告像一致。
% 只重新作用当前模型，不运行LSQR、原NS或候选推进；冻结原整箭试验证据不改。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
p=jsondecode(fileread(fullfile(root,'outputs/v84_validation/q24system/peak_system.json')));
c=load(p.artifacts.cache.file,'cache','rate');rows=cell(1,2);
for k=1:2
    segment=k+1;s=load(fullfile(root,sprintf('outputs/v84_validation/q24direction%d/peak_zero_candidate_evidence.mat',segment)),'value');
    [selected,r]=selectRocketMinimumStateDirection(c.cache,s.value,c.rate,.0999);
    actual=applyCutCellSelectedDirection(c.cache,selected);
    difference=norm(actual-r.selectedModel,'fro')/max(norm(actual,'fro'),realmin);
    assert(difference<1e-12);
    if segment==2
        assert(~r.applied&&isequal(selected,s.value)&&isequal(r.selectedModel,r.baselineModel));
        assert(isequal(r.returnedCoefficients,r.normalization));
    else
        saved=load(fullfile(root,'outputs/v84_validation/full02/baseline/step1/selected_direction.mat'),'selected');
        assert(r.applied&&isequal(selected,saved.selected));
        assert(isequal(r.returnedCoefficients,r.coefficients));
    end
    rows{k}=struct('segment',segment,'applied',r.applied,'modelRelativeDifference',difference, ...
        'returnedModelDefect',r.returnedModelDefect,'minimumModelDefectInSpan',r.minimumModelDefectInSpan);
end
report=struct('status','passed-returned-direction-model-consistency-and-feasible-direction-bitwise-replay', ...
    'cases',{rows},'newOriginalNSCalls',0,'newLSQRIterations',0,'newFlowAcceptances',0);
fid=fopen(fullfile(folder,'returned_model.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
