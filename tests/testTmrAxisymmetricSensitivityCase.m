function report=testTmrAxisymmetricSensitivityCase(folder,parentEvidenceFile)
%TESTTMRAXISYMMETRICSENSITIVITYCASE 原节点四档工况、逻辑初猜和独立初始账的定向校核。
% 模块一：逐节点核对官方level3抽样、各向加密和真实墙面，不以单元计数代替几何证据。
% 模块二：独立双循环核对逻辑父复制全部五量，并重新积分新网格初始守恒量。
% 模块三：固定气体/边界/SA/尺度与原稳态门槛，核对初始空间残差和五量边界/源账。
% 模块四：非法档位、路径、父状态/策略/尺度/旧结果必须拒绝；不改冻结v66证据。
% 本测试没有推进或求解新稳态，不是网格收敛、实验或生产面积验收。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v67_validation','sensitivity_cases');end
if nargin<2
    parentEvidenceFile=fullfile(root,'outputs','v66_validation','tmr_wall_thermal','tmr_wall_thermal_evidence.mat');
    if ~exist(parentEvidenceFile,'file'),parentEvidenceFile=fullfile(root,'validation','v66','tmr_wall_thermal','tmr_wall_thermal_evidence.mat');end
end
assert(exist(parentEvidenceFile,'file')==2,'AxisFVTest:TmrSensitivitySource','缺少冻结v66父检查点。');
parentSHA=fileSHA(parentEvidenceFile);parent=load(parentEvidenceFile,'mesh','gas','boundary','implicit','report');
dataRoot=fullfile(root,'tests','data','tmr_mach7');fine=loadTmrMach7Grid(3,dataRoot);coarse=loadTmrMach7Grid(4,dataRoot);
assert(isequal(fine.X(1:2:end,1:2:end),coarse.X)&&isequal(fine.R(1:2:end,1:2:end),coarse.R));
names={'baseline','streamwise','normal','both'};factorI=[1 2 1 2];factorJ=[1 1 2 2];
expectedNodes=[41 26;81 26;41 51;81 51];expectedCells=[1000 2000 2000 4000];cases=cell(1,4);records=struct([]);
refusals=struct('label',{},'identifier',{});

%% 一、真实原节点和墙面保留；法向加密不偷偷改变流向，反之亦然
for k=1:4
    item=createTmrAxisymmetricSensitivityCase(names{k},parentEvidenceFile,dataRoot);mesh=item.mesh;g=item.gridStudy;
    assert(isequal(sort(fieldnames(item)),sort({'mesh';'gas';'boundary';'Q';'residualScale';'parentBinding';'gridStudy'})));
    is=1:(2/factorI(k)):81;js=1:(2/factorJ(k)):51;
    assert(isequal(size(mesh.X),expectedNodes(k,:))&&mesh.cellCount==expectedCells(k)&& ...
        isequal(mesh.X,fine.X(is,js))&&isequal(mesh.R,fine.R(is,js))&& ...
        isequal(g.officialLevel3IIndices,is)&&isequal(g.officialLevel3JIndices,js)&& ...
        g.factorI==factorI(k)&&g.factorJ==factorJ(k), ...
        'AxisFVTest:TmrSensitivityNodes','%s不是指定方向的原节点抽样。',names{k});
    assert(all(mesh.planarArea>0)&&all(mesh.volumePerRadian>0)&&all(mesh.faces.normalDistance>0)&& ...
        all(diff(mesh.X,1,1)<0,'all')&&all(diff(mesh.R,1,2)>0,'all')&& ...
        isequal(mesh.X(1:factorI(k):end,1:factorJ(k):end),coarse.X)&& ...
        isequal(mesh.R(1:factorI(k):end,1:factorJ(k):end),coarse.R)&& ...
        isequal(mesh.X(:,1),fine.X(is,1))&&isequal(mesh.R(:,1),fine.R(is,1))&& ...
        g.wallFaceCount==expectedNodes(k,1)-1&&g.parentWallNodesPreservedExactly, ...
        'AxisFVTest:TmrSensitivityWall','墙面或几何正性未保留。');

    %% 二、独立逻辑父双循环；只作初猜，新的体积分可以与父网格不同
    copied=zeros(mesh.cellCount,5);parentIDs=zeros(mesh.cellCount,1);
    for j=1:mesh.shape(2)
        for i=1:mesh.shape(1)
            row=(j-1)*mesh.shape(1)+i;
            parentI=1+floor((i-1)/factorI(k));parentJ=1+floor((j-1)/factorJ(k));parentRow=(parentJ-1)*40+parentI;
            copied(row,:)=parent.implicit.Q(parentRow,:);parentIDs(row)=parentRow;
        end
    end
    assert(isequal(item.Q,copied)&&isequal(g.parentCellID,parentIDs)&&all(item.Q(:,5)>=0)&& ...
        all(accumarray(parentIDs,1,[1000 1])==factorI(k)*factorJ(k))&& ...
        ~g.conservativeRemapPerformed&&g.initializationCopiesAllFiveConservedComponents, ...
        'AxisFVTest:TmrSensitivityCopy','逻辑父单元或全部五量复制不正确。');
    initialIntegral=zeros(1,5);parentIntegral=zeros(1,5);
    for c=1:mesh.cellCount,initialIntegral=initialIntegral+item.Q(c,:)*mesh.volumePerRadian(c);end
    for c=1:1000,parentIntegral=parentIntegral+parent.implicit.Q(c,:)*parent.mesh.volumePerRadian(c);end
    integralBudget=64*eps*(sum(abs(item.Q.*mesh.volumePerRadian),1)+sum(abs(parent.implicit.Q.*parent.mesh.volumePerRadian),1));
    assert(all(abs(initialIntegral-g.initialConservedIntegralPerRadian)<=integralBudget)&& ...
        all(abs(parentIntegral-g.parentConservedIntegralPerRadian)<=integralBudget)&& ...
        all(abs(initialIntegral-parentIntegral-g.initialIntegralDifferencePerRadian)<=integralBudget), ...
        'AxisFVTest:TmrSensitivityIntegral','新网格初始账未从真实新体积与复制Q重建。');
    grouped=zeros(1000,1);for c=1:mesh.cellCount,grouped(parentIDs(c))=grouped(parentIDs(c))+mesh.volumePerRadian(c);end
    assert(isequal(grouped-parent.mesh.volumePerRadian,g.childGroupedVolumeDifferencePerParentM3));

    %% 三、固定合同、原始残差和独立有单位的五量空间账；不推进时间
    assert(isequaln(item.gas,parent.gas)&&isequaln(item.boundary,parent.boundary)&& ...
        isequal(item.residualScale,parent.implicit.residualScale)&& ...
        strcmp(item.parentBinding.sha256,parentSHA)&&item.parentBinding.completeSpatialReplayPassed&& ...
        isequaln(item.parentBinding.gas,parent.gas)&&isequaln(item.parentBinding.boundary,parent.boundary)&& ...
        isequal(item.parentBinding.residualScale,item.residualScale)&& ...
        g.originalSteadyCriteria.steadyResidualTolerance==2e-6&&g.originalSteadyCriteria.steadyRequiredSteps==3&& ...
        g.originalSteadyCriteria.minimumSteadySteps==3&&g.originalSteadyCriteria.conservationTolerance==5e-10);
    [rate,info]=axisymmetricRANSResidual(item.Q,mesh,item.gas,item.boundary);spatialLedger=ledger(rate,info,mesh);
    assert(strcmp(info.sa.modelID,'SA-noft2-nonnegative-AJS2012-v1')&&info.sa.wallMuTIsZero&& ...
        all(info.wallThermalGradient.tangentialGradientInputKm==0)&& ...
        ~info.physicalSolution&&~info.eligibleForSurfaceCoverage&& ...
        g.timeStepsTaken==0&&g.physicalTimeAdvancedS==0&&~g.newGridSteadySolved&& ...
        ~g.gridConvergenceEstablished&&~g.physicalSolution&&~g.experimentalValidation&&g.newProductionSolvedAreaM2==0);
    if k==1
        assert(isequaln(mesh,parent.mesh)&&isequal(item.Q,parent.implicit.Q)&&isequaln(info,parent.implicit.finalInfo)&& ...
            max(abs(rate)./item.residualScale,[],'all')==parent.implicit.residualNorm&&all(g.initialIntegralDifferencePerRadian==0));
    end
    record=struct('caseID',names{k},'nodeShape',size(mesh.X),'cellCount',mesh.cellCount, ...
        'factorI',factorI(k),'factorJ',factorJ(k),'wallFaces',g.wallFaceCount, ...
        'initialIntegralDifferencePerRadian',g.initialIntegralDifferencePerRadian, ...
        'volumeDifferencePerRadianM3',g.volumeDifferencePerRadianM3, ...
        'maximumGroupedVolumeRelativeDifference',g.maximumGroupedVolumeRelativeDifference, ...
        'initialOriginalResidual',max(abs(rate)./item.residualScale,[],'all'),'initialSpatialLedger',spatialLedger, ...
        'newSteadySolved',false,'timeStepsTaken',0,'conservativeRemapPerformed',false);
    if isempty(records),records=record;else,records(end+1)=record;end %#ok<AGROW>
    cases{k}=item;
end
% 两个单方向对照必须相互不同，不能只改标签/单元数或在两方向同时重采样。
assert(isequal(cases{2}.mesh.X(1:2:end,:),cases{1}.mesh.X)&& ...
    isequal(cases{3}.mesh.X(:,1:2:end),cases{1}.mesh.X)&& ...
    isequal(cases{4}.mesh.X(:,1:2:end),cases{2}.mesh.X)&& ...
    isequal(cases{4}.mesh.X(1:2:end,:),cases{3}.mesh.X));

%% 四、非法输入和陈旧父数据必须先拒绝；临时副本与旧v66证据分开
for value={42,true,NaN,{'normal'},['normal';'normal'],"bad",'stream','normal '}
    reject('invalid-case-name',@()createTmrAxisymmetricSensitivityCase(value{1},parentEvidenceFile,dataRoot),'AxisFV:TmrSensitivityInput');
end
reject('invalid-parent-path-type',@()createTmrAxisymmetricSensitivityCase('normal',42,dataRoot),'AxisFV:TmrSensitivityInput');
reject('missing-parent',@()createTmrAxisymmetricSensitivityCase('normal',[tempname '.mat'],dataRoot),'AxisFV:TmrSensitivitySource');
reject('invalid-grid-root-type',@()createTmrAxisymmetricSensitivityCase('normal',parentEvidenceFile,42),'AxisFV:TmrSensitivityInput');
temporary=[tempname '_RocketAeroHeat_TmrSensitivityTest'];mkdir(temporary);cleanup=onCleanup(@()removeTemporary(temporary)); %#ok<NASGU>
badPath=fullfile(temporary,'bad_parent.mat');
bad=parent;bad.implicit.steadyConverged=false;saveBad(badPath,bad);
reject('nonsteady-parent',@()createTmrAxisymmetricSensitivityCase('normal',badPath,dataRoot),'AxisFV:TmrSensitivityParent');
bad=parent;bad.implicit.Q(1,5)=-1;saveBad(badPath,bad);
reject('negative-SA-parent',@()createTmrAxisymmetricSensitivityCase('normal',badPath,dataRoot),'AxisFV:TmrSensitivityParent');
bad=parent;bad.boundary.jMin.thermalGradient='legacy';saveBad(badPath,bad);
reject('legacy-wall-parent',@()createTmrAxisymmetricSensitivityCase('normal',badPath,dataRoot),'AxisFV:TmrSensitivityParent');
bad=parent;bad.implicit.options.steadyResidualTolerance=1e-4;saveBad(badPath,bad);
reject('loosened-steady-parent',@()createTmrAxisymmetricSensitivityCase('normal',badPath,dataRoot),'AxisFV:TmrSensitivityParent');
bad=parent;bad.implicit.residualScale(4)=2*bad.implicit.residualScale(4);saveBad(badPath,bad);
reject('changed-scale-parent',@()createTmrAxisymmetricSensitivityCase('normal',badPath,dataRoot),'AxisFV:TmrSensitivityParent');
bad=parent;bad.implicit.Q(1,4)=1.001*bad.implicit.Q(1,4);saveBad(badPath,bad);
reject('stale-state-vs-result-parent',@()createTmrAxisymmetricSensitivityCase('normal',badPath,dataRoot),'AxisFV:TmrSensitivityParent');
assert(strcmp(parentSHA,fileSHA(parentEvidenceFile)),'AxisFVTest:TmrSensitivityParentChanged','测试改写了冻结v66父证据。');
report=struct('status','passed-four-directional-original-node-grid-cases-and-initialization-contract', ...
    'records',records,'refusals',refusals,'parentEvidenceSHA256',parentSHA,'parentEvidenceUnchanged',true, ...
    'sourceBindings',[binding(fullfile(root,'regionalFlow','axisymmetric','createTmrAxisymmetricSensitivityCase.m'), ...
        'regionalFlow/axisymmetric/createTmrAxisymmetricSensitivityCase.m'),binding([mfilename('fullpath') '.m'], ...
        'tests/testTmrAxisymmetricSensitivityCase.m')], ...
    'timeStepsTaken',0,'newSteadySolved',false,'gridConvergenceEstablished',false, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'sensitivity_cases_evidence.mat'),'report','cases','-v7');
fid=fopen(fullfile(folder,'sensitivity_cases_summary.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear closer;
fprintf('TMR sensitivity cases passed: original nodes 1000/2000/2000/4000, five-component logical initialization, fixed contracts, %d refusals; no new steady or grid-convergence claim.\n',numel(refusals));
    function reject(label,action,expected)
        try,action();catch exception
            assert(strcmp(exception.identifier,expected),'AxisFVTest:TmrSensitivityUnexpectedRefusal','%s期望%s，实际%s。',label,expected,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFVTest:TmrSensitivityMissingRefusal','非法合同%s没有拒绝。',label);
    end
end

function audit=ledger(rate,info,mesh)
% 模块用途：每个不同单位分量分别核对边界/源账，包含共享通量相消的舍入预算。
volumeTerms=rate.*mesh.volumePerRadian;flux=info.faceFluxIntegralPerRadian;source=info.cellSourceIntegral;
boundary=flux(mesh.boundaryFaces,:);defect=sum(volumeTerms,1)+sum(boundary,1)-sum(source,1);
physical=sum(abs(volumeTerms),1)+sum(abs(boundary),1)+sum(abs(source),1);
count=full(sum(mesh.cellFaceMultiplicity,1))';allAbsolute=physical+sum(abs(flux).*count,1);
k=2*mesh.cellCount+numel(mesh.boundaryFaces)+16;gamma=k*eps/(1-k*eps);budget=2e-12*physical+gamma*allAbsolute;
assert(all(abs(defect)<=budget),'AxisFVTest:TmrSensitivityLedger','新初始状态五量账不闭合。');
audit=struct('defect',defect,'physicalScale',physical,'combinedBudget',budget,'perRadian',true,'unitSeparated',true,'passed',true);
end

function saveBad(path,item)
% 模块用途：只写本测试私有临时副本，不覆盖原父证据。
save(path,'-struct','item','-v7');
end

function record=binding(path,relative)
% 模块用途：绑定实际完整MATLAB文件字节，避免仅以文件名声称测试版本。
record=struct('path',relative,'sha256',fileSHA(path));
end

function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function removeTemporary(path)
% 模块用途：只删除已确认位于系统临时根目录的本测试私有路径。
absolute=char(java.io.File(path).getCanonicalPath());temporary=char(java.io.File(tempdir).getCanonicalPath());
if startsWith(lower(absolute),[lower(temporary) filesep])&&endsWith(absolute,'_RocketAeroHeat_TmrSensitivityTest')&&exist(absolute,'dir')
    rmdir(absolute,'s');
end
end
