function study=loadNonsimilarSphericalStudy(studyFile,inletFile)
%LOADNONSIMILARSPHERICALSTUDY 读取冻结球鼻条件试验并复核每个保存站。
% 模块用途：将原outputs草稿的数值证据变为可移植的研究输入，不重跑长段。
% 结果不是实验验证或生产解。源码与旧运行的绑定是事后重算，不能伪称
% 旧运行开始前已锁定源码。初始局部相似入口不算非相似求解站。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
folder=fullfile(root,'tests','data','nonsimilar_sphere');
if nargin<1,studyFile=fullfile(folder,'sphere_n121_40_60.mat');end
if nargin<2,inletFile=fullfile(folder,'sphere_local_similarity_inlet.mat');end
% 冻结证据可搬到别处，但内容必须与显式清单一致；不暗中接受改写结果。
manifest=jsondecode(fileread(fullfile(folder,'source_manifest.json')));
if ~strcmp(shaFile(studyFile),manifest.files(1).sha256)|| ...
        ~strcmp(shaFile(inletFile),manifest.files(2).sha256)
    error('RocketAeroHeat:NonSimilarStudy:FrozenInput','输入不是清单绑定的冻结条件试验。');
end
raw=load(studyFile,'report','profiles','steps','gas','o');
seed=load(inletFile,'fs','result');
required={'report','profiles','steps','gas','o'};
if ~all(isfield(raw,required))||~all(isfield(seed,{'fs','result'}))
    error('RocketAeroHeat:NonSimilarStudy:Schema','缺少原始剖面、步骤、气体或入口数据。');
end
report=raw.report;profiles=raw.profiles;steps=raw.steps;fs=seed.fs;
for key={'Mach','pressure','temperature','speed','R','gamma','alpha','beta'}
    validateattributes(fs.(key{1}),{'numeric'},{'scalar','real','finite'});
end
validateattributes(fs.unitFlowDirection,{'numeric'},{'real','finite','numel',3});
if fs.Mach<=1||fs.pressure<=0||fs.temperature<=0||fs.R<=0|| ...
        raw.gas.R~=fs.R||raw.gas.gamma~=fs.gamma
    error('RocketAeroHeat:NonSimilarStudy:Condition','入口气体和自由流不一致或不是本超声速范围。');
end
if ~flagIs(report.stoppedEarly,false)||~iscell(profiles)||~iscell(steps)|| ...
        numel(profiles)~=numel(steps)+1||numel(steps)<2|| ...
        any(cellfun(@isempty,profiles))||any(cellfun(@isempty,steps))|| ...
        report.acceptedNumericalSteps~=numel(steps)
    error('RocketAeroHeat:NonSimilarStudy:Incomplete','仅接受完整、已结束的条件试验。');
end
inletSHA=shaFile(inletFile);
if ~strcmp(inletSHA,report.inlet.sourceSha256)
    error('RocketAeroHeat:NonSimilarStudy:InletBinding','入口源文件与试验记录的SHA256不符。');
end
if ~flagIs(report.uniformWallHistory,true)||~strcmp(raw.gas.thermalModel,'thermallyPerfect')|| ...
        fs.alpha~=0||fs.beta~=0||norm(fs.unitFlowDirection-[-1 0 0])>1e-12
    error('RocketAeroHeat:NonSimilarStudy:Scope','本数据接口只支持同一均匀壁温历史、NASA7和零攻角球鼻。');
end
R=seed.result.noseRadius;validateattributes(R,{'numeric'},{'scalar','real','finite','positive'});
Tw=report.wallTemperatureInputK;station=struct([]);
validateattributes(Tw,{'numeric'},{'scalar','real','finite','>=',80,'<=',2500});
validateattributes(raw.o.angularStepDeg,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(raw.o.startAngleDeg,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(raw.o.endAngleDeg,{'numeric'},{'scalar','real','finite','>',raw.o.startAngleDeg});
for j=1:numel(steps)
    step=steps{j};c=profiles{j+1};history=profiles{j};
    if j>1,history=[history,profiles{j-1}];end
    if ~flagIs(step.available,true)||~flagIs(step.converged,true)||~isequaln(c,step.profile)|| ...
            c.Tw~=Tw||~isequal(c.n,profiles{1}.n)
        error('RocketAeroHeat:NonSimilarStudy:Profile','保存站没有一致的已收敛剖面及壁温历史。');
    end
    d=nonsimilarLaminarResidual(c,history,raw.gas);
    angle=rad2deg(c.s/R);row=report.rows(j);
    validateattributes(step.wallHeatFluxIntoWall,{'numeric'},{'scalar','real','finite'});
    validateattributes(row.angleDeg,{'numeric'},{'scalar','real','finite'});
    validateattributes(row.deltaOverNoseRadius,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(row.ReynoldsS,{'numeric'},{'scalar','real','finite','positive'});
    if ~flagIs(row.numericalThinLayerChecks,true)
        error('RocketAeroHeat:NonSimilarStudy:Applicability','原数值薄层筛选未通过或状态非法。');
    end
    scale=max([abs(d.wallHeatFluxIntoWall),abs(step.wallHeatFluxIntoWall),1]);
    expectedAngle=min(raw.o.startAngleDeg+j*raw.o.angularStepDeg,raw.o.endAngleDeg);
    if abs(angle-row.angleDeg)>1e-9||abs(angle-expectedAngle)>1e-9|| ...
            abs(c.r-R*sind(angle))>1e-10*R|| ...
            abs(d.wallHeatFluxIntoWall-step.wallHeatFluxIntoWall)>1e-9*scale|| ...
            norm(d.scaledResidual,inf)>1e-7|| ...
            max(abs(d.wallBoundaryResidual))>1e-10||max(abs(d.edgeBoundaryResidual))>1e-8
        error('RocketAeroHeat:NonSimilarStudy:Replay','当前内核与保存的剖面/壁导数/残差不一致。');
    end
    normalPower=c.r*d.wallHeatFluxIntoWall;
    item=struct('angleDeg',angle,'heatFlux',d.wallHeatFluxIntoWall, ...
        'wallShearStress',d.wallShearStress,'edgePressure',c.p, ...
        'edgeTemperature',c.Te,'edgeSpeed',c.Ue,'wallTemperatureInput',c.Tw, ...
        'normalizedResidual',norm(d.scaledResidual,inf), ...
        'energyDefectWPerM',d.conservativeEnergyIntegral.defectWPerM, ...
        'energyDefectOverWallPower',d.conservativeEnergyIntegral.defectWPerM/max(abs(normalPower),realmin), ...
        'thinLayerChecksPassed',row.numericalThinLayerChecks, ...
        'deltaOverNoseRadius',row.deltaOverNoseRadius,'ReynoldsS',row.ReynoldsS, ...
        'positiveWallShear',d.wallShearStress>0,'profileIndex',j+1);
    if isempty(station),station=item;else,station(end+1)=item;end %#ok<AGROW>
end
if any(diff([station.angleDeg])<=0)|| ...
        ~all([station.thinLayerChecksPassed])||~all([station.positiveWallShear])
    error('RocketAeroHeat:NonSimilarStudy:Applicability','保存站没有通过原条件数值筛选。');
end
study=struct('schemaVersion',1,'kind','conditional-spherical-nonsimilar-study', ...
    'freestream',fs,'noseRadius',R,'gas',raw.gas,'uniformWallTemperature',Tw, ...
    'inletAngleDeg',raw.o.startAngleDeg,'maximumAngularStepDeg',raw.o.angularStepDeg, ...
    'acceptedIntervalDeg',[station(1).angleDeg,station(end).angleDeg], ...
    'station',station,'inlet',report.inlet,'profiles',{profiles}, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false, ...
    'newSolvedAreaM2',0,'historicalPreRunSourceBinding',false, ...
    'currentSourceReplayPassed',true,'recoveryTemperatureSolved',false, ...
    'wallTemperatureSolved',false,'externalCFDRequired',false, ...
    'source',struct('studyFile',char(studyFile),'studySHA256',shaFile(studyFile), ...
        'inletFile',char(inletFile),'inletSHA256',inletSHA), ...
    'method','TP-3271物理s/n非相似层流，真实上游u/T，BDF1启动后BDF2，Newton隐式推进', ...
    'scope','给定局部相似有限角入口与固定熵工程外缘的条件研究；未验证真实激波层、区域接口或实验精度');
currentPaths={'regionalFlow/nonsimilar/nonsimilarLaminarResidual.m','airThermalProperties.m'};
study.currentSourceHashes=struct([]);
for k=1:numel(currentPaths)
    study.currentSourceHashes(k).path=currentPaths{k};
    study.currentSourceHashes(k).sha256=shaFile(fullfile(root,currentPaths{k}));
end
study.replayBindingSHA256=nonsimilarStudyDigest(study);
end

function yes=flagIs(value,expected)
yes=(islogical(value)||isnumeric(value))&&isreal(value)&&isscalar(value)&& ...
    isfinite(value)&&ismember(value,[0 1])&&logical(value)==expected;
end

function hash=shaFile(path)
% 文件绑定使用原始二进制；文件路径变化不改变证据身份。
fid=fopen(path,'rb');assert(fid>=0,'RocketAeroHeat:NonSimilarStudy:File','文件不能读取。');
clean=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
