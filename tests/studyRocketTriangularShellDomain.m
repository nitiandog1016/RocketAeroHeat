function report=studyRocketTriangularShellDomain(folder,options)
%STUDYROCKETTRIANGULARSHELLDOMAIN 在真实19752原面箭体上构造完整三维外域并保存证据。
% 模块一：只读原几何和原工况，未读取旧压力/热流作为三维体场入口。
% 模块二：近壁/过渡/外域均为真实体积，保留原面编号及给定壁温身份。
% 模块三：本次仅验收几何和实际代价，不将其称作20度流场或物理解完成。
if nargin<2,options=struct;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(~exist(folder,'dir'),'RocketShellStudy:Overwrite','不能覆盖已有试验。');
source=fullfile(root,'outputs','v56_validation','state_validity','state_validity_actual_rocket.mat');
if ~exist(source,'file'),source=fullfile(root,'validation','inherited_v56','v56','state_validity','state_validity_actual_rocket.mat');end
old=load(source,'result','fs','cfg');V=old.result.V;F=old.result.F;
assert(size(F,1)==19752&&old.fs.alpha==20&&old.fs.Mach==5);
t=tic;domain=buildRocketTriangularShellDomain(V,F,options);elapsed=toc(t);
assert(isequal(V,domain.inputGeometry.V)&&isequal(F,domain.inputGeometry.F));
assert(isequal(domain.walls.originalFaceID,(1:size(F,1))')&&isequal(domain.walls.ownerCell,(1:size(F,1))'));
assert(max(abs(domain.walls.area-old.result.area)./old.result.area)<3e-13);
fs=old.fs;cfg=old.cfg;report=domain.quality;
report.status='constructed-real-original-triangle-three-dimensional-domain-not-flow-solved';
report.sourceSHA256=fileSHA(source);report.sourceFlowVelocityVector=fs.flowVelocityVector;
report.angleOfAttackDeg=20;report.fullThreeDimensionalMomentumRequired=true;
report.oldEngineeringPressureOrHeatUsedAsVolumeInlet=false;report.elapsedSeconds=elapsed;
report.engineProcessID=feature('getpid');report.flowSolved=false;
mkdir(folder);save(fullfile(folder,'shell_domain_evidence.mat'),'domain','report','fs','cfg','-v7');
fid=fopen(fullfile(folder,'shell_domain_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function value=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
