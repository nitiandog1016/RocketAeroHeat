function report = testRocketRemeshFeatures(outputDirectory)
%TESTROCKETREMESHFEATURES 核对完整箭体突起在既有重网格中的实际保留程度。
% 模块用途：以真实三角形射线求交比较同一 (X,theta) 位置的径向尺寸，
% 并检查闭合、有向边、连通分量和 Euler 示性数；不以最大半径替代特征保真。
% 默认参数与显式加密参数分别记录，只有显式加密结果参加本测试的通过门槛。
% 本测试只验证几何和网格，不作遮挡、分离或气动热模型的物理验证。
if nargin<1,outputDirectory='';end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if ~isempty(outputDirectory) && ~exist(outputDirectory,'dir'),mkdir(outputDirectory);end

%% 1. 建立光滑整箭及同一母线上的两个不相交突起
% 阈值先于计算固定：RMS <= 10% 高度，最大误差 <= 20%，中心高度保留 >= 85%。
features(1)=struct('type','fairing','height',0.0192,'length',0.24, ...
    'centerFraction',0.38,'azimuthDeg',35,'halfWidthDeg',22);
features(2)=struct('type','ridge','height',0.0192,'length',0.50, ...
    'centerFraction',0.63,'azimuthDeg',210,'halfWidthDeg',12);
options=struct('resolution',72,'protrusions',features);
[rawFeatured,truth]=buildRocketTestBody(options);
options.protrusionEnabled=false;[rawBaseline,baselineTruth]=buildRocketTestBody(options);
report=struct('geometryOnly',true,'physicalValidation',false, ...
    'newPhysicalSolvedArea',0,'coordinateConvention','SI; nose +X; flow -X', ...
    'rawBaseline',topologyAudit(rawBaseline),'rawFeatured',topologyAudit(rawFeatured), ...
    'thresholds',struct('relativeRms',0.10,'relativeMaximum',0.20, ...
    'minimumCenterRetention',0.85,'maximumCenterRetention',1.10), ...
    'truth',truth.features,'baselineBodyRadius',baselineTruth.bodyRadius);
requireClosed(report.rawBaseline);requireClosed(report.rawFeatured);

%% 2. 独立探针：沿突起的 X/周向支持域采样，先确认源几何确实含目标形状
% patch 含支持域以外的点，可发现特征偏移、错转和突起外泄。
patch=cell(numel(features),1);
for k=1:numel(features)
    f=truth.features(k);z=linspace(-1.1,1.1,13);
    xx=f.centerX+f.length/2*z;tt=deg2rad(f.azimuthDeg+f.halfWidthDeg*z);
    original=radialSamples(rawFeatured,xx,tt)-radialSamples(rawBaseline,xx,tt);
    analytic=analyticFeature(f,xx,tt);
    sourceError=original-analytic;
    assert(sqrt(mean(sourceError(:).^2))<0.025*f.height, ...
        '源三角网格未充分解析突起，不能把重网格结果当成保真证据。');
    patch{k}=struct('x',xx,'thetaRad',tt,'rawHeight',original, ...
        'analyticHeight',analytic,'sourceRelativeRms',sqrt(mean(sourceError(:).^2))/f.height);
end

%% 3. 保持主配置与重网格程序原样，分别实测默认及显式加密设置
baseCfg=createCaseConfig();fineCfg=baseCfg;
fineCfg.remesh.hMin=0.002;fineCfg.remesh.hMax=0.008;
fineCfg.remesh.initialAxialSpacing=0.008;
fineCfg.remesh.maxCircumPoints=360;fineCfg.remesh.maxFaces=200000;
configs={baseCfg,fineCfg};labels={'default','featureResolved'};
for j=1:numel(configs)
    timer=tic;
    remeshedBaseline=autoRemeshSurface(rawBaseline,configs{j});
    remeshedFeatured=autoRemeshSurface(rawFeatured,configs{j});
    result=struct('label',labels{j},'config',configs{j}.remesh, ...
        'baselineTopology',topologyAudit(remeshedBaseline), ...
        'featuredTopology',topologyAudit(remeshedFeatured), ...
        'baselineQuality',meshQuality(remeshedBaseline,configs{j},struct([])), ...
        'featuredQuality',meshQuality(remeshedFeatured,configs{j},truth.features),'elapsedSeconds',toc(timer));
    requireClosed(result.baselineTopology);requireClosed(result.featuredTopology);
    for k=1:numel(features)
        p=patch{k};f=truth.features(k);
        h=radialSamples(remeshedFeatured,p.x,p.thetaRad)- ...
            radialSamples(remeshedBaseline,p.x,p.thetaRad);
        delta=h-p.rawHeight;middle=(size(h,1)+1)/2;
        metrics=struct('type',f.type,'relativeRms',sqrt(mean(delta(:).^2))/f.height, ...
            'relativeMaximum',max(abs(delta(:)))/f.height, ...
            'sourceRelativeRms',p.sourceRelativeRms, ...
            'centerRetention',h(middle,middle)/p.rawHeight(middle,middle), ...
            'x',p.x,'thetaRad',p.thetaRad,'rawHeight',p.rawHeight,'remeshedHeight',h);
        metrics.passed=metrics.relativeRms<=report.thresholds.relativeRms && ...
            metrics.relativeMaximum<=report.thresholds.relativeMaximum && ...
            metrics.centerRetention>=report.thresholds.minimumCenterRetention && ...
            metrics.centerRetention<=report.thresholds.maximumCenterRetention;
        result.features(k)=metrics;
        fprintf('%s %s: RMS/H %.6f, max/H %.6f, center retained %.6f, pass %d\n', ...
            labels{j},f.type,metrics.relativeRms,metrics.relativeMaximum,metrics.centerRetention,metrics.passed);
    end
    result.allFeatureGatesPassed=all([result.features.passed]);
    if ~isempty(outputDirectory)
        result.meshEvidenceFile=['remesh_',labels{j},'_meshes.mat'];
        save(fullfile(outputDirectory,result.meshEvidenceFile),'remeshedBaseline','remeshedFeatured');
    end
    report.(labels{j})=result;
end

%% 4. 记录明确边界并交付原始采样矩阵，避免把几何测试误称物理解面积
report.limitations={ ...
    'Only smooth cylinder-mounted fairings/ridges represented by a single radial graph were tested.', ...
    'No independent fins, undercuts, intersecting shells or multiple-loop X sections are certified.', ...
    'Topology and geometric retention do not validate shadowing, separation, pressure or heat flux.', ...
    'Geometric retention does not certify triangle quality; actual aspect ratio and angle failures are recorded.', ...
    'Default mesh settings must be assessed per feature; explicit fine settings do not alter defaults.'};
report.passed=report.featureResolved.allFeatureGatesPassed;
if ~isempty(outputDirectory)
    if ~exist(outputDirectory,'dir'),mkdir(outputDirectory);end
    save(fullfile(outputDirectory,'remesh_features.mat'),'report');
    fid=fopen(fullfile(outputDirectory,'remesh_features.json'),'w','n','UTF-8');
    assert(fid>=0,'无法写出几何试验记录。');guard=onCleanup(@()fclose(fid));
    fprintf(fid,'%s',jsonencode(report));clear guard;
end
assert(report.passed,'显式加密后的突起未满足预先规定的几何保真门槛。');
fprintf('Rocket remesh geometry tests passed; new physical solved area = 0.\n');
end

function audit=topologyAudit(mesh)
%% 有向拓扑审计：共享边恰两面且方向相反，单连通闭合箭体 Euler=2
V=mesh.V;F=mesh.F;E=[F(:,[1,2]);F(:,[2,3]);F(:,[3,1])];
[~,~,eid]=unique(sort(E,2),'rows');count=accumarray(eid,1);
direction=accumarray(eid,sign(E(:,2)-E(:,1)));
G=graph(E(:,1),E(:,2),[],size(V,1));component=conncomp(G);
A=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2)/2;
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
audit=struct('numberVertices',size(V,1),'numberFaces',size(F,1), ...
    'boundaryEdges',sum(count==1),'nonManifoldEdges',sum(count>2), ...
    'inconsistentEdges',sum(direction~=0),'connectedComponents',max(component), ...
    'eulerCharacteristic',size(V,1)-numel(count)+size(F,1), ...
    'minimumTriangleArea',min(vecnorm(A,2,2)),'volume',volume);
end

function q=meshQuality(mesh,cfg,features)
%% 原样记录全局面元质量，不因拓扑/突起保真通过而宣称网格达到工程质量
V=mesh.V;F=mesh.F;P=V(F(:,1),:);Q=V(F(:,2),:);S=V(F(:,3),:);
a=vecnorm(Q-S,2,2);b=vecnorm(S-P,2,2);c=vecnorm(P-Q,2,2);
area=vecnorm(cross(Q-P,S-P,2),2,2)/2;
rin=2*area./max(a+b+c,eps);ar=max([a,b,c],[],2)./max(2*sqrt(3)*rin,eps);
aa=acosd(max(-1,min(1,(b.^2+c.^2-a.^2)./max(2*b.*c,eps))));
ab=acosd(max(-1,min(1,(c.^2+a.^2-b.^2)./max(2*c.*a,eps))));
angle=min([aa,ab,180-aa-ab],[],2);center=(P+Q+S)/3;
q=struct('meanAspectRatio',mesh.remeshInfo.meanAspectRatio, ...
    'maximumAspectRatio',mesh.remeshInfo.maxAspectRatio, ...
    'meanMinimumAngleDeg',mesh.remeshInfo.meanMinimumAngle, ...
    'minimumAngleDeg',mesh.remeshInfo.minimumAngle, ...
    'requestedMaxAspectRatio',cfg.remesh.maxAspectRatio, ...
    'requestedMinAngleDeg',cfg.remesh.minAngle);
q.requestedQualityGatesPassed=q.maximumAspectRatio<=q.requestedMaxAspectRatio && ...
    q.minimumAngleDeg>=q.requestedMinAngleDeg;
q.qualityConfiguredMeets=q.requestedQualityGatesPassed;
bad=ar>q.requestedMaxAspectRatio | angle<q.requestedMinAngleDeg;
q.numberFailingFaces=sum(bad);q.failingAreaFraction=sum(area(bad))/sum(area);
[~,i]=max(ar);q.worstAspectFace=faceRecord(i);
[~,i]=min(angle);q.worstAngleFace=faceRecord(i);
for k=1:numel(features)
    f=features(k);t=atan2(center(:,3),center(:,2));
    d=atan2(sin(t-deg2rad(f.azimuthDeg)),cos(t-deg2rad(f.azimuthDeg)));
    selected=abs(center(:,1)-f.centerX)<f.length/2 & abs(d)<deg2rad(f.halfWidthDeg);
    r=struct('type',f.type,'selection','face centroid inside exact X-theta support', ...
        'numberFaces',sum(selected),'maximumAspectRatio',max(ar(selected)), ...
        'minimumAngleDeg',min(angle(selected)),'numberFailingFaces',sum(bad & selected), ...
        'failingAreaFraction',sum(area(bad & selected))/sum(area(selected)));
    q.featureRegions(k)=r;
end
    function r=faceRecord(i)
        r=struct('faceIndex',i,'sectionBand',mesh.param.faceBandID(i), ...
            'centerXYZ',center(i,:),'verticesXYZ',V(F(i,:),:), ...
            'aspectRatio',ar(i),'minimumAngleDeg',angle(i));
    end
end

function requireClosed(a)
assert(a.boundaryEdges==0 && a.nonManifoldEdges==0 && a.inconsistentEdges==0 && ...
    a.connectedComponents==1 && a.eulerCharacteristic==2 && ...
    a.minimumTriangleArea>0 && a.volume>0,'箭体不满足单连通有向闭合曲面约定。');
end

function radius=radialSamples(mesh,xx,tt)
%% 实际三角面 Moller-Trumbore 射线求交，射线起点为轴线上的 [X,0,0]
% 仅检查当前 X 所跨三角形，重复共边命中合并；多次不同半径命中明确拒绝。
V=mesh.V;F=mesh.F;p=V(F(:,1),:);a=V(F(:,2),:)-p;b=V(F(:,3),:)-p;
lo=min(reshape(V(F(:),1),size(F)),[],2);hi=max(reshape(V(F(:),1),size(F)),[],2);
radius=nan(numel(xx),numel(tt));scale=max(max(V)-min(V));
for i=1:numel(xx)
    selected=lo<=xx(i)+1e-12*scale & hi>=xx(i)-1e-12*scale;
    P=p(selected,:);A=a(selected,:);B=b(selected,:);S=[xx(i),0,0]-P;
    Q=cross(S,A,2);
    for j=1:numel(tt)
        D=repmat([0,cos(tt(j)),sin(tt(j))],size(P,1),1);
        H=cross(D,B,2);det=sum(A.*H,2);valid=abs(det)>1e-14*scale^2;
        U=sum(S.*H,2)./det;W=sum(D.*Q,2)./det;R=sum(B.*Q,2)./det;
        hit=valid & U>=-1e-10 & W>=-1e-10 & U+W<=1+1e-10 & R>0;
        values=sort(R(hit));assert(~isempty(values),'径向射线未命中箭体。');
        assert(all(abs(values-values(1))<1e-8*scale), ...
            '射线命中多个不同半径，不能按单径向图解释该模型。');
        radius(i,j)=mean(values);
    end
end
end

function h=analyticFeature(f,xx,tt)
%% 突起解析高程，仅用于确认源三角网格自身解析了指定的形状
[X,T]=ndgrid(xx,tt);xi=2*(X-f.centerX)/f.length;
eta=atan2(sin(T-deg2rad(f.azimuthDeg)),cos(T-deg2rad(f.azimuthDeg)))/deg2rad(f.halfWidthDeg);
inside=abs(xi)<1 & abs(eta)<1;h=zeros(size(X));
if strcmp(f.type,'fairing'),fx=0.5*(1+cos(pi*xi(inside)));
else,q=abs(xi(inside));fx=ones(size(q));ramp=q>0.5;fx(ramp)=0.5*(1+cos(2*pi*(q(ramp)-0.5)));end
h(inside)=f.height*fx.*(0.5*(1+cos(pi*eta(inside))));
end
