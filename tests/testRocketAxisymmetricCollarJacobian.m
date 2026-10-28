function report=testRocketAxisymmetricCollarJacobian(folder,sourceFolder)
%TESTROCKETAXISYMMETRICCOLLARJACOBIAN 核对新真实近壁网格已接受端点的完整五量Jacobian。
% 模块一：只读已保存七步前缀，重建工况及终态空间方程，不修改真实链。
% 模块二：按生产默认2e-6构造完整有色矩阵，独立逐列复算鼻部/后体五量的25列。
% 模块三：检查这些列的数值及保守稀疏支持；这是选列核对，不冒称全矩阵串行验证。
% 中央割线跨非光滑分支时不等于唯一解析导数，本测试不宣称分支光滑或物理精度。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v70_validation','collar_jacobian');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v70_validation','normal_collar');end
path=fullfile(sourceFolder,'rocket_axisymmetric_evidence.mat');before=fileSHA(path);
s=load(path,'caseDefinition','implicit');c=s.caseDefinition;last=s.implicit;
assert(isequaln(c,createRocketAxisymmetricCollarCase(c.options))&&last.acceptedSteps>=1);
[r,i]=axisymmetricRANSResidual(last.Q,c.mesh,c.gas,c.boundary);assert(isequaln(i,last.finalInfo));
Q=last.Q;U=Q./c.residualScale;nc=c.mesh.cellCount;
[J,details,pattern]=axisymmetricFVJacobian(Q,c.mesh,c.gas,c.boundary,c.residualScale,struct);
assert(details.success&&details.columnsDifferentiated==numel(Q)&&details.differenceStep==2e-6&& ...
    details.retries==0,'RocketCollarJacobian:Build','需要无重试的默认中央差分完整矩阵作为明确对照。');
cells=[1 2 105 106 107];columns=reshape(cells(:)+(0:4)*nc,[],1);serial=zeros(numel(Q),numel(columns));
errors=zeros(numel(columns),1);outside=zeros(numel(columns),1);h=2e-6*max(1,abs(U(:)));
for k=1:numel(columns)
    col=columns(k);a=U;b=U;a(col)=a(col)+h(col);b(col)=b(col)-h(col);
    rp=axisymmetricRANSResidual(a.*c.residualScale,c.mesh,c.gas,c.boundary)./c.residualScale;
    rm=axisymmetricRANSResidual(b.*c.residualScale,c.mesh,c.gas,c.boundary)./c.residualScale;
    serial(:,k)=(rp(:)-rm(:))/(2*h(col));
    errors(k)=max(abs(full(J(:,col))-serial(:,k)))/max(max(abs(serial(:,k))),realmin);
    cellID=mod(col-1,nc)+1;allowed=repmat(full(pattern.cellPattern(:,cellID))~=0,5,1);
    outside(k)=max([0;abs(serial(~allowed,k))]);
end
assert(all(errors<2e-6)&&all(outside<1e-6), ...
    'RocketCollarJacobian:SelectedColumns','真实鼻/后体五量选列或保守支持与独立串行不符。');
assert(strcmp(before,fileSHA(path))&&max(abs(r)./c.residualScale,[],'all')==last.residualNorm);
report=struct('status','passed-real-collar-five-quantity-default-Jacobian-selected-serial-columns', ...
    'matlabEngineProcessID',feature('getpid'),'sourceEvidenceSHA256',before, ...
    'acceptedPrefixSteps',last.acceptedSteps,'fullMatrixSize',size(J),'selectedCells',cells, ...
    'selectedColumnCount',numel(columns),'relativeColumnError',errors,'outsidePatternMaximum',outside, ...
    'differenceStep',details.differenceStep,'entireMatrixSerialVerified',false, ...
    'analyticSmoothDerivativeProved',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'collar_jacobian_evidence.mat'),'report','details','pattern','columns','serial','J');
fid=fopen(fullfile(folder,'collar_jacobian_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Real collar full five-quantity Jacobian: %d selected columns, maximum relative error %.6g.\n',numel(columns),max(errors));
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));v=fread(fid,Inf,'*uint8');clear guard;
h=java.security.MessageDigest.getInstance('SHA-256');h.update(v);digest=lower(reshape(dec2hex(typecast(h.digest(),'uint8'),2)',1,[]));
end
