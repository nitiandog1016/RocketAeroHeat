function report=verifyRocketWlsNormalBlocks(folder,resultCase,parentCase)
%VERIFYROCKETWLSNORMALBLOCKS 在真实原CV单位方向上复核五带提取，不更新流场。
% 模块一：修正矩阵属于新试算初始Q，原图来自同Q父端点，不用新接受态混验旧块。
% 模块二：两个真实CV的全部五量方向，逐个比较同原面完整WLS差与矩阵补丁。
% 模块三：十个方向只作模型验证，原NS/主链计数不增加。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'),fullfile(root,'tests'));setupCutCellFlow;
assert(~exist(folder,'dir'),'WlsNormalVerify:Overwrite','不得覆盖真实块验证。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation',resultCase,'wall_velocity_krylov_evidence.mat');
parentFile=fullfile(root,'outputs','v83_validation',parentCase,'wall_velocity_krylov_evidence.mat');
s=load(sourceFile,'result','gas','boundary');r=s.result;gas=s.gas;boundary=s.boundary;Q=r.initialQ;scale=r.residualScale;operator=r.firstOperator;
baseInfo=struct;for key={'primitive','gradient','viscousFaceFlux'},baseInfo.(key{1})=r.initialInfo.(key{1});end;clear r s;
p=load(parentFile,'result');assert(isequaln(p.result.Q,Q));graph=p.result.firstOperator;clear p;
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;clear d;
cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);nf=operator.sourceFaces;cells=[19623 78879];errors=zeros(2,5);evidence=cell(2,5);
for a=1:numel(cells),for b=1:5
    directionU=zeros(size(Q));directionU(cells(a),b)=1;[full,approx]=splitRocketGraphDirection(cache,graph,directionU);
    fullValue=full.momentumStress+full.heatConduction+full.workGradient+full.workVelocity;
    approxValue=approx.momentumStress+approx.heatConduction+approx.workGradient+approx.workVelocity;
    same=mod((1:size(Q,1))'-1,nf)==mod(cells(a)-1,nf);expected=fullValue-approxValue;expected(~same,:)=0;
    actual=applyRocketShellGraphOperator(graph,directionU)-applyRocketWlsNormalOperator(operator,directionU);
    errors(a,b)=norm(actual-expected,'fro')/max(1,norm(expected,'fro'));
    assert(errors(a,b)<2e-9,'WlsNormalVerify:CurrentAction','真实CV五量块不等于原完整WLS同线投影。');
    evidence{a,b}=struct('cellID',cells(a),'equationID',b,'actual',actual,'fullViscousAction',fullValue,'approxViscousAction',approxValue);
    fprintf('WLS_NORMAL_PROBE cell=%d equation=%d error=%.6g\n',cells(a),b,errors(a,b));
end,end
report=struct('status','passed-ten-real-original-CV-conservative-directions-WLS-same-line-block-extraction', ...
    'cellIDs',cells,'relativeDifferences',errors,'maximumRelativeDifference',max(errors,[],'all'), ...
    'completeSelectedReferenceApplications',10,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'wls_normal_real_blocks.mat'),'report','sourceFile','parentFile','cells','errors','evidence','-v7');
fid=fopen(fullfile(folder,'wls_normal_real_blocks_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
