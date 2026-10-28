function report=probeWallVelocityRejectedInterface(folder,sourceFolder,segmentIteration,alpha)
%PROBEWALLVELOCITYREJECTEDINTERFACE 还原已保存真实方向的界面温域拒绝，保留原护栏。
% 模块一：读取真实保存方向，确认单元rho/内能/温度仍合法，不把单元合法当界面合法。
% 模块二：原精确Riemann批量调用给出真实失败面，映回原内部面/owner/层数。
% 模块三：仅诊断拒绝，不裁温、换通量或将候选作为接受状态。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
% 模块四：缺省保留旧第五次半步诊断；新接续须显式指定实际来源及候选。
if nargin<2,sourceFolder=fullfile(root,'outputs','v83_validation','krylov02');end
if nargin<3,segmentIteration=2;end
if nargin<4,alpha=.5;end
validateattributes(segmentIteration,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(alpha,{'numeric'},{'scalar','real','positive','finite','<=',1});
assert(~exist(folder,'dir'),'WallReject:Overwrite','不得覆盖真实诊断。');mkdir(folder);
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
s=load(fullfile(sourceFolder,'wall_velocity_krylov_evidence.mat'),'result','report','gas');
assert(segmentIteration<=numel(s.result.trace),'WallReject:Source','指定方向没有保存的接受trace。');
t=s.result.trace{segmentIteration};rows=s.result.lineLog;
recorded=[rows.iteration]==segmentIteration&[rows.pseudoDt]==t.pseudoDt&[rows.alpha]==alpha;
assert(nnz(recorded)==1&&~rows(recorded).accepted&&strcmp(rows(recorded).errorIdentifier,'CutExactEuler:Temperature'), ...
    'WallReject:Source','指定候选必须是原记录的未接受界面温域拒绝。');
Q=t.initialQ+alpha*t.linearUpdateU.*s.result.residualScale;gas=s.gas;
rho=Q(:,1);u=Q(:,2:4)./rho;T=(Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*gas.R/(gas.gamma-1));
assert(all(rho>0)&&all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)));
f=d.domain.openFaces;ids=find(f.area>0&f.cellMinus>0&f.cellPlus>0);identifier='';message='';
try
    cutCellExactEulerFlux(Q(f.cellMinus(ids),:),Q(f.cellPlus(ids),:),f.normal(ids,:),gas);
catch ex
    identifier=ex.identifier;message=ex.message;
    if ~strcmp(identifier,'CutExactEuler:Temperature'),rethrow(ex);end
end
assert(strcmp(identifier,'CutExactEuler:Temperature'),'WallReject:Source','没有复现记录中的原界面温域拒绝。');
token=regexp(message,'批次面(\d+)','tokens','once');assert(~isempty(token));batch=str2double(token{1});faceID=ids(batch);
nf=size(d.domain.inputGeometry.F,1);owner=f.cellMinus(faceID);neighbour=f.cellPlus(faceID);
ql=Q(owner,:);qr=Q(neighbour,:);normal=f.normal(faceID,:);
report=struct('status','observed-real-rejected-interface-temperature-with-valid-cell-state', ...
    'newEquationIteration',segmentIteration,'segmentIteration',segmentIteration,'alpha',alpha,'errorIdentifier',identifier,'message',message, ...
    'sourceEvidenceFolder',sourceFolder,'actualOriginalInterfaceCallCount',1,'fullNSResidualRecomputed',false, ...
    'globalInteriorFaceID',faceID,'batchFaceIndex',batch,'centroidM',f.centroid(faceID,:), ...
    'ownerCV',owner,'neighbourCV',neighbour,'ownerLayer',floor((owner-1)/nf)+1,'neighbourLayer',floor((neighbour-1)/nf)+1, ...
    'ownerSourceFacet',mod(owner-1,nf)+1,'neighbourSourceFacet',mod(neighbour-1,nf)+1, ...
    'candidateCellTemperatureRangeK',[min(T),max(T)],'candidateAccepted',false, ...
    'originalTemperatureRangeKept',true,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
if isfield(s.report,'priorAcceptedNonlinearIterations'),report.newEquationIteration=segmentIteration+s.report.priorAcceptedNonlinearIterations;end
save(fullfile(folder,'rejected_interface_evidence.mat'),'report','ql','qr','normal','gas','-v7');
fid=fopen(fullfile(folder,'rejected_interface_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
