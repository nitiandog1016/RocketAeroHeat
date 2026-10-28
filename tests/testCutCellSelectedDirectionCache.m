function report=testCutCellSelectedDirectionCache(folder,kind)
%TESTCUTCELLSELECTEDDIRECTIONCACHE 缓存作用对已独立核对的冷作用及线性叠加，不重复NS试算。
% 模块一：单位两方向或真实失败方向使用冻结证据，缓存本身不产生新NS调用。
% 模块二：变换系数后核对固定线性作用，未知物理资格不从缓存测试继承。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'DirectionCacheTest:Overwrite','不得覆盖方向缓存证据。');
if strcmp(kind,'unit'),name='wlsunit01';file='wls_direction_unit_evidence.mat';else,assert(strcmp(kind,'real'));name='wlsrocket01';file='wls_direction_real_evidence.mat';end
sourceFile=fullfile(root,'outputs','v83_validation',name,file);s=load(sourceFile,'Q','baseInfo','cut','gas','boundary','scale','allEvidence');
Q=s.Q;baseInfo=s.baseInfo;cut=s.cut;gas=s.gas;boundary=s.boundary;scale=s.scale;
directions=cell(numel(s.allEvidence),1);expected=directions;
for j=1:numel(directions),directions{j}=s.allEvidence{j}.directionU;expected{j}=s.allEvidence{j}.viscous+s.allEvidence{j}.euler;end
clear s;t=tic;cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);buildSeconds=toc(t);
errors=zeros(numel(directions),1);applySeconds=errors;values=cell(numel(directions),1);
for j=1:numel(directions)
    t=tic;values{j}=applyCutCellSelectedDirection(cache,directions{j});applySeconds(j)=toc(t);
    errors(j)=norm(values{j}(:)-expected{j}(:))/max(1,norm(expected{j}(:)));
end
assert(max(errors)<1e-10,'DirectionCacheTest:Cold','缓存没有复现已独立核对的冷作用。');
other=directions{end};if numel(directions)==1,other=other.*reshape(cos((1:numel(other))*.37),size(other));end
second=applyCutCellSelectedDirection(cache,other);combined=applyCutCellSelectedDirection(cache,.7*directions{1}-1.3*other);
reference=.7*values{1}-1.3*second;linearity=norm(combined(:)-reference(:))/max(1,norm(reference(:)));
assert(linearity<1e-10,'DirectionCacheTest:Linearity','缓存不是固定基态的线性作用。');
report=struct('status','passed-selected-complete-direction-cache-cold-reference-and-linearity','kind',kind, ...
    'cellCount',size(Q,1),'buildSeconds',buildSeconds,'applySeconds',applySeconds,'coldActionRelativeDifferences',errors, ...
    'linearSuperpositionRelativeDifference',linearity,'newOriginalNSEvaluations',0,'newAcceptedIterations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'selected_direction_cache_evidence.mat'),'report','sourceFile','directions','expected','values','other','second','combined','-v7');
fid=fopen(fullfile(folder,'selected_direction_cache_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
