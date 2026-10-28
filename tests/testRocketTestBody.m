function report = testRocketTestBody()
%TESTROCKETTESTBODY 核对箭体夹具的闭合拓扑、几何精度、突起真值和拒绝合同。
% 这里是几何验证；不等同于激波、分离、热流的物理或实验验证。
addpath(fileparts(fileparts(mfilename('fullpath'))));
cases={'none','fairing','ridge'};checks=0;records=cell(0,1);
for j=1:numel(cases)
    err=zeros(2,1);areaErr=zeros(2,1);
    for k=1:2
        cfg=struct('resolution',24*2^(k-1),'protrusionType',cases{j});
        [m,t]=buildRocketTestBody(cfg);checkTopology(m,t);checks=checks+1;
        err(k)=abs(t.geometryAudit.enclosedVolume/t.exactVolume-1);
        if j==1,areaErr(k)=abs(sum(m.area)/t.exactArea-1);end
        assert(all(isfinite(m.normal(:)))&&all(abs(vecnorm(m.normal,2,2)-1)<1e-12));
        assert(t.geometryOnly && ~t.physicalValidation);
        assert(isequal(size(t.faceFeatureMask),[size(m.F,1),numel(t.features)]));
        assert(all(diff(m.param.sectionX)>0)&&all(m.param.sectionPerimeter>0));
        assert(numel(m.param.faceBandID)==size(m.F,1));
        if j>1
            f=t.features(1);peak=max(t.nodeFeatureHeight(:,1));
            assert(abs(peak-f.height)<1e-13,'没有采到突起峰值。');
            assert(sum(t.faceFeatureMask(:,1))>30,'突起没有保留为真实三角面。');
            ids=find(t.faceFeatureMask(:,1));
            angular=abs(atan2(sin(t.faceParamTheta(ids)-deg2rad(f.azimuthDeg)),cos(t.faceParamTheta(ids)-deg2rad(f.azimuthDeg))));
            assert(all(abs(t.faceParamX(ids)-f.centerX)<f.length/2));
            assert(all(angular<deg2rad(f.halfWidthDeg)));
            assert(all(strcmp(t.faceBaseRegion(ids),'cylinder')));
            assert(all(strcmp(t.faceRegion(ids),['protrusion_',f.type])));
        else
            assert(all(t.nodeRadialOffset==0)&&isempty(t.features));
        end
        records{end+1}=struct('type',cases{j},'resolution',cfg.resolution, ...
            'vertices',size(m.V,1),'faces',size(m.F,1),'relativeVolumeError',err(k), ...
            'peakHeights',t.geometryAudit.featurePeakHeight,'audit',t.geometryAudit); %#ok<AGROW>
    end
    assert(err(2)<err(1)&&err(2)<0.012,'细化未改善解析体积误差。');checks=checks+1;
    if j==1,assert(areaErr(2)<areaErr(1)&&areaErr(2)<0.012);checks=checks+1;end
    fprintf('%s: volume relative error %.6g -> %.6g\n',cases{j},err(1),err(2));
end

%% 2. 两类不重叠突起、跨周期接缝和重叠支持域分别核查
f=struct('type','fairing','height',0.018,'length',0.24,'centerFraction',0.35,'azimuthDeg',358,'halfWidthDeg',18);
f(2)=struct('type','ridge','height',0.012,'length',0.30,'centerFraction',0.68,'azimuthDeg',175,'halfWidthDeg',10);
[multi,mt]=buildRocketTestBody(struct('resolution',48,'protrusions',f));checkTopology(multi,mt);
assert(all(sum(mt.faceFeatureMask,1)>0)&&~mt.featuresOverlap&&isfinite(mt.exactVolume));
assert(all(abs(max(mt.nodeFeatureHeight,[],1)-[f.height])<1e-13));checks=checks+1;
assert(any(mt.faceFeatureMask(mt.faceParamTheta<0.05,1))&&any(mt.faceFeatureMask(mt.faceParamTheta>6.2,1)));checks=checks+1;
f(2)=f(1);f(2).azimuthDeg=4;
[overlap,ot]=buildRocketTestBody(struct('resolution',24,'protrusions',f));checkTopology(overlap,ot);
assert(ot.featuresOverlap&&isnan(ot.exactVolume)&&any(sum(ot.faceFeatureMask,2)>1));checks=checks+1;
[disabled,dt]=buildRocketTestBody(struct('resolution',24,'protrusions',f,'protrusionEnabled',false));
[baseline,bt]=buildRocketTestBody(struct('resolution',24,'protrusionType','none'));
assert(isequal(disabled.V,baseline.V)&&isequal(disabled.F,baseline.F)&&isempty(dt.features));checks=checks+1;
assert(mt.geometryAudit.enclosedVolume>bt.exactBaselineVolume*0.97); % 两网格分开，不冒称逐节点相同。

%% 3. 对不支持的输入明确拒绝，避免配置错误静默改变测试几何
bad={struct('resolution',15),struct('resolution',24.5),struct('bodyRadius',NaN), ...
    struct('bodyRadius',0.01,'baseRadius',0.08),struct('noseRadius',1), ...
    struct('coneAngleDeg',90),struct('protrusionEnabled',1),struct('typo',1), ...
    struct('protrusions',struct('type','fin')), ...
    struct('protrusions',struct('height',-1)), ...
    struct('protrusions',struct('length',2)), ...
    struct('protrusions',struct('halfWidthDeg',90)), ...
    struct('protrusions',struct('azimuthDeg',1i)), ...
    struct('protrusions',struct('centerFraction',0)), ...
    struct('protrusions',struct('unknown',1))};
for j=1:numel(bad)
    rejected=false;try,buildRocketTestBody(bad{j});catch,rejected=true;end
    assert(rejected,'无效输入 %d 未拒绝。',j);checks=checks+1;
end
report=struct('status','passed','checks',checks,'cases',[records{:}],'physicalValidation',false, ...
    'scope','Topology, analytic geometry convergence and exact protrusion support only.');
fprintf('Rocket geometry fixture: %d checks passed; physical validation not performed.\n',checks);
end

function checkTopology(m,t)
% 独立重算边与有向体积，不把生成器自报的审核当作唯一证据。
F=m.F;V=m.V;edges=[F(:,[1,2]);F(:,[2,3]);F(:,[3,1])];
[~,~,id]=unique(sort(edges,2),'rows');
assert(all(accumarray(id,1)==2),'开口或非流形边。');
assert(all(accumarray(id,sign(edges(:,2)-edges(:,1)))==0),'边方向不一致。');
assert(size(V,1)-max(id)+size(F,1)==2,'Euler 示性数不是2。');
assert(numel(unique(F))==size(V,1),'存在未使用顶点。');
normalArea=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
assert(all(vecnorm(normalArea,2,2)>0));
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
assert(volume>0&&abs(volume-t.geometryAudit.enclosedVolume)<1e-13);
base=strcmp(t.faceBaseRegion,'base');assert(all(m.normal(base,1)<-1+1e-12));
side=~base;radial=sum(m.normal(side,2:3).*m.center(side,2:3),2);
assert(all(radial>0),'侧壁法向未朝外。');
assert(t.geometryAudit.connectedComponents==1);
assert(all(diff(m.param.sectionX)>1e-13*max(1,t.totalLength)),'近重复轴向节点产生了数值退化细片。');
end
