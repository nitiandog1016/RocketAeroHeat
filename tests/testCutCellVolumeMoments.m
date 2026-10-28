function report=testCutCellVolumeMoments(folder)
%TESTCUTCELLVOLUMEMOMENTS 四个冻结几何对独立oracle及局部闭式矩的验证。
% 模块用途：只测试新增MATLAB二阶矩；不重构旧四个几何、不运行流场推进。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,folder=fullfile(root,'outputs','v60_validation','library_moments');end
if ~isfolder(folder),mkdir(folder);end
maxNumCompThreads(1);addpath(fullfile(root,'internalFlow'));setupCutCellFlow;started=tic;
v58=validationFolder(root,'v58');v59=validationFolder(root,'v59');
names={'actual','sphere_grid8','sphere_grid12','sphere_grid16'};cases=cell(4,1);checks=cell(4,1);
for j=1:4
    if j==1
        d=load(fullfile(v58,'cut_cell_draft/cut_cell_geometry_results.mat'),'actual');source=d.actual;
        d=load(fullfile(v58,'cut_cell_transport_draft/actual_aggregation_results.mat'),'group');group=d.group;
    else
        d=load(fullfile(v59,'geometry',[names{j},'.mat']),'source','group');source=d.source;group=d.group;
    end
    d=load(fullfile(v59,'curved_accuracy_oracle',[names{j},'_fv_volume_moments.mat']),'moments');oracle=d.moments;
    fprintf('MATLAB_SECOND_MOMENT_START %s\n',names{j});
    [moments,details]=buildCutCellVolumeMoments(source,group);
    covError=max(abs(moments.covariance-oracle.covariance),[],'all');
    centralError=max(abs(moments.centralSecondMoment-oracle.centralSecondMoment),[],'all');
    volumeError=max(abs(moments.volume-oracle.volume(:)));centroidError=max(abs(moments.centroid-oracle.centroid),[],'all');
    assert(covError<3e-12 && centralError<3e-15 && volumeError<3e-15 && centroidError<3e-12, ...
        'CutMomTest:IndependentOracle','MATLAB prefix帽二阶矩与独立截面环oracle不一致。');
    full=group.cells.state==2;fullError=0;
    for id=find(full)'
        fullError=max(fullError,max(abs(reshape(moments.covariance(id,:,:),3,3)-diag(source.grid.spacing.^2/12)),[],'all'));
    end
    assert(fullError<2e-15);
    checks{j}=struct('name',names{j},'groups',numel(moments.volume),'covarianceOracleMaximumDifference',covError, ...
        'centralSecondOracleMaximumDifference',centralError,'volumeOracleMaximumDifference',volumeError, ...
        'centroidOracleMaximumDifference',centroidError,'fullCubeCovarianceMaximumDifference',fullError, ...
        'minimumCovarianceEigenvalue',details.minimumCovarianceEigenvalue,'elapsedSeconds',details.elapsedSeconds);
    cases{j}=struct('name',names{j},'moments',moments,'details',details);
    fprintf('MATLAB_SECOND_MOMENT_PASS %s groups=%d CovError=%.9g elapsed=%.6g\n', ...
        names{j},numel(moments.volume),covError,details.elapsedSeconds);
end

%% 闭式体积交集：盒、斜切simplex、共面/接触域、单格包裹实体
d=load(fullfile(v58,'cut_cell_draft/analytic_cut_cell_results.mat'),'analyticCases','fixtures');
selected=[1 2 3 6 7];analytic=cell(numel(selected),1);analyticCases=cell(numel(selected),1);
for j=1:numel(selected)
    index=selected(j);source=d.analyticCases{index};[moments,details]=buildCutCellVolumeMoments(source);
    covariance=exactFluidCovariance(source,index==3);valid=source.cells.volume>0;
    error=max(abs(moments.covariance(valid,:,:)-covariance(valid,:,:)),[],'all');
    assert(error<3e-11,'CutMomTest:ClosedForm','闭式流体二阶矩不一致。');
    analytic{j}=struct('name',d.fixtures{index}.name,'covarianceClosedFormMaximumDifference',error, ...
        'minimumPositivePatchArea',min(source.walls.area),'positivePatchCount',numel(source.walls.area), ...
        'wallRelativeAreaError',details.maximumWallRelativeAreaError);
    analyticCases{j}=struct('name',d.fixtures{index}.name,'moments',moments,'details',details);
end

%% 大坐标平移：几何局部坐标不变，不以全局二阶原矩相减制造Cov
base=d.analyticCases{1};[before,beforeDetails]=buildCutCellVolumeMoments(base); %#ok<ASGLU>
shift=[1024 -2048 4096];translated=translateCut(base,shift);[after,afterDetails]=buildCutCellVolumeMoments(translated); %#ok<ASGLU>
valid=base.cells.volume>0;translationError=max(abs(before.covariance(valid,:,:)-after.covariance(valid,:,:)),[],'all');
assert(translationError<3e-12,'CutMomTest:Translation','局部Cov大坐标平移敏感。');

%% 拒绝合同：离域/错误方向/损坏原法向、面面积、组体积及无流体
rejections=cell(6,1);
bad=base;bad.inputGeometry.V(1,1)=base.grid.origin(1)-base.grid.spacing(1);
rejections{1}=reject(@()buildCutCellVolumeMoments(bad),'CutMom:OutsideDomain');
bad=base;bad.inputGeometry.F=bad.inputGeometry.F(:,[1 3 2]);
rejections{2}=reject(@()buildCutCellVolumeMoments(bad),'CutMom:Orientation');
bad=base;bad.walls.normal(1,:)=2*bad.walls.normal(1,:);
rejections{3}=reject(@()buildCutCellVolumeMoments(bad),'CutMom:WallNormal');
bad=base;id=find(bad.openFaces.area>0,1);bad.openFaces.area(id)=.5*bad.openFaces.area(id);
rejections{4}=reject(@()buildCutCellVolumeMoments(bad),'CutMom:ApertureMismatch');
s=load(fullfile(v59,'geometry/sphere_grid8.mat'),'source','group');badGroup=s.group;badGroup.cells.volume(1)=1.01*badGroup.cells.volume(1);
rejections{5}=reject(@()buildCutCellVolumeMoments(s.source,badGroup),'CutMom:GroupVolume');
rejections{6}=reject(@()buildCutCellVolumeMoments(d.analyticCases{8}),'CutMom:NoFluid');
report=struct('status','passed-tested-MATLAB-prefix-second-volume-moments','actualAndSphere',{checks}, ...
    'analyticClosedForm',{analytic},'translationCovarianceMaximumDifference',translationError, ...
    'translation',shift,'rejections',{rejections},'elapsedSeconds',toc(started), ...
    'independentPythonOracleUsedOnlyInValidation',true,'runtimeRequiresExternalOracle',false, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false);
save(fullfile(folder,'matlab_volume_moment_evidence.mat'),'report','cases','analyticCases','before','after','-v7');
fid=fopen(fullfile(folder,'matlab_volume_moment_summary.json'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));fprintf('MATLAB_VOLUME_MOMENT_TEST_PASS elapsed=%.6g\n',toc(started));
end

function p=validationFolder(root,version)
p=fullfile(root,'outputs',[version,'_validation']);
if strcmp(version,'v58'),required=fullfile('cut_cell_draft','cut_cell_geometry_results.mat');
else,required=fullfile('geometry','sphere_grid8.mat');end
if ~isfile(fullfile(p,required)),p=fullfile(root,'validation',version);end
end

function output=exactFluidCovariance(c,oblique)
% 模块用途：独立盒/截simplex体积积分；不使用待测平面边界恒等式。
V=c.inputGeometry.V;bodyLo=min(V,[],1);bodyHi=max(V,[],1);threshold=max(sum(V,2));
h=c.grid.spacing;o=c.grid.origin;index=c.cells.index;N=numel(c.cells.volume);output=nan(N,3,3);
for id=find(c.cells.volume>0)'
    lo=o+(index(id,:)-1).*h;hi=lo+h;reference=c.cells.centroid(id,:);
    [vf,f1,f2]=boxMoment(lo,hi,reference);a=max(lo,bodyLo);b=min(hi,bodyHi);
    if oblique,[vs,s1,s2]=simplexMoment(a,b,threshold,reference);else,[vs,s1,s2]=boxMoment(a,b,reference);end
    volume=vf-vs;mean=(f1-s1)/volume;C=(f2-s2)/volume-mean'*mean;output(id,:,:)=reshape(C,1,3,3);
end
end
function [v,m1,m2]=boxMoment(lo,hi,q)
if any(hi<=lo),v=0;m1=zeros(1,3);m2=zeros(3);return;end
v=prod(hi-lo);r=(lo+hi)/2-q;m1=v*r;m2=v*(r'*r+diag((hi-lo).^2/12));
end
function [v,m1,m2]=simplexMoment(lo,hi,threshold,q)
v=0;m1=zeros(1,3);m2=zeros(3);if any(hi<=lo) || threshold<=sum(lo),return;end
if threshold>=sum(hi),[v,m1,m2]=boxMoment(lo,hi,q);return;end
width=hi-lo;
for mask=0:7
    active=double(bitget(mask,1:3));t=max(0,threshold-sum(lo+active.*width));sgn=(-1)^sum(active);r=lo+active.*width-q;
    a=t^3/6;b=ones(1,3)*t^4/24;S=(ones(3)+eye(3))*t^5/120;
    v=v+sgn*a;m1=m1+sgn*(b+a*r);m2=m2+sgn*(S+r'*b+b'*r+a*(r'*r));
end
end
function c=translateCut(c,shift)
c.grid.origin=c.grid.origin+shift;c.inputGeometry.V=c.inputGeometry.V+shift;
for key={'center','centroid','rawCentroid'},if isfield(c.cells,key{1}),c.cells.(key{1})=c.cells.(key{1})+shift;end,end
for key={'centroid','rawCentroid'},if isfield(c.openFaces,key{1}),c.openFaces.(key{1})=c.openFaces.(key{1})+shift;end,end
c.walls.centroid=c.walls.centroid+shift;c.walls.grid.origin=c.walls.grid.origin+shift;
for k=1:numel(c.walls.polygonVertices),c.walls.polygonVertices{k}=c.walls.polygonVertices{k}+shift;end
end
function label=reject(fn,expected)
got='';try,fn();catch err,got=err.identifier;end
assert(strcmp(got,expected),'CutMomTest:Rejection','期待%s，实际%s。',expected,got);label=got;
end
