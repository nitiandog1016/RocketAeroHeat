function report=testParabolicShockLayer(folder,sourceFile)
%TESTPARABOLICSHOCKLAYER 跨层质量、同流线熵、非正交度量及真实20度箭体校核。
% 模块一：同一共同激波和同一真实自由流求层内压力/温度/速度/厚度。
% 模块二：检查流线来源、全列能量/熵、积分加密及刚体变换，不套旧壁面场。
% 模块三：在原19752面上实际迭代候选形状，匹配失败不能冒充物理解完成。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupEngineeringShockLayer;
assert(~exist(folder,'dir'),'ShockLayer:Overwrite','跨层研究须写新目录。');
d=load(sourceFile,'result','fs');geometry=d.result;fs=d.fs;
assert(size(geometry.F,1)==19752&&fs.alpha==20&&fs.Mach==5);
upstream=struct('pressure',double(fs.pressure),'temperature',double(fs.temperature), ...
    'density',double(fs.density),'velocityVector',double(fs.flowVelocityVector(:)), ...
    'method','原19752面箭体M5/30km/alpha20均匀自由流；未使用旧压力/外缘/NS状态');
ex=upstream.velocityVector/norm(upstream.velocityVector);ey=[0;1;0];ey=ey-ex*dot(ey,ex);ey=ey/norm(ey);
frame=[ex ey cross(ex,ey)];nose=estimateNoseRadius(geometry);assert(nose.available);
depth=max(geometry.V(:,1))-geometry.V(:,1);points=double(geometry.V(depth<=nose.outerDepthRadiusFraction*nose.radius,:));
center=mean(points,1);scale=max(vecnorm(points-center,2,2));z=(points-center)/scale;
fit=[2*z ones(size(z,1),1)]\sum(z.^2,2);noseCenter=(center+scale*fit(1:3)')';
radius=nose.radius;shape=struct('origin',noseCenter-(radius+.15*radius)*ex, ...
    'frame',frame,'K',[.72 .08;.08 .95]/radius);
eta=linspace(1,0,11)';q=[.022;.013];modes={'caloricallyPerfect','thermallyPerfect'};cases=cell(2,1);
for j=1:2
    gas=createSurfaceFlowGas(modes{j},double(fs.R),1.4);
    c=solveParabolicShockLayerColumn(upstream,gas,shape,q,eta,struct('quadratureOrder',16));
    assert(all(diff(c.normalDistanceM)>0)&&all(c.jacobianRatio>0));
    assert(max(abs(c.state.energyResidual))<1e-10&&max(abs(c.state.entropyResidual))<1e-8);
    assert(max(c.state.sourceEntropy)-min(c.state.sourceEntropy)>1);
    assert(max(abs(c.offDiagonalMetricRatio))>1e-5,'ShockLayer:TestMetric','本例必须检验真实非正交项。');
    assert(~c.bodyMatched&&~c.boundaryLayerEdgeAvailable&&~c.productionQualified);
    fine=solveParabolicShockLayerColumn(upstream,gas,shape,q,eta,struct('quadratureOrder',32));
    assert(max(abs(c.normalDistanceM-fine.normalDistanceM))<1e-8);
    % 同一物理流线在不同列有不同eta，但共同激波来源和熵必须相同。
    [Q,D]=eig(shape.K);tau=.009;downstreamQ=Q*diag(exp(diag(D)*tau))*Q'*q;
    firstEta=.6;secondEta=firstEta*exp(-trace(shape.K)*tau);
    a=solveParabolicShockLayerColumn(upstream,gas,shape,q,[1;firstEta;0]);
    b=solveParabolicShockLayerColumn(upstream,gas,shape,downstreamQ,[1;secondEta;0]);
    assert(norm(a.state.sourceShockCoordinate(2,:)-b.state.sourceShockCoordinate(2,:))<1e-14);
    assert(abs(a.state.entropy(2)-b.state.entropy(2))<1e-9);
    % 刚体变换必须只改变三维位置和速度方向，不改变层内标量或厚度。
    angle=.47;rotation=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];shift=[.8;-.3;1.1];
    moved=shape;moved.origin=rotation*shape.origin+shift;moved.frame=rotation*shape.frame;
    changed=upstream;changed.velocityVector=rotation*upstream.velocityVector;
    transformed=solveParabolicShockLayerColumn(changed,gas,moved,q,eta);
    assert(max(abs(transformed.normalDistanceM-c.normalDistanceM))<1e-10);
    assert(max(vecnorm(transformed.velocityVector-c.velocityVector*rotation',2,2))<1e-7);
    assert(max(vecnorm(transformed.position-(c.position*rotation'+shift'),2,2))<1e-11);
    cases{j}=struct('column',c,'refined',fine,'sameStreamlineUpstream',a,'sameStreamlineDownstream',b, ...
        'downstreamParameterM',tau,'transformed',transformed,'rotation',rotation,'translation',shift);
end

% 拟合前保存跨层计算检查点；后续报告/拟合问题不得抹掉已算出的列。
mkdir(folder);save(fullfile(folder,'columns_before_fit.mat'),'cases','geometry','fs','nose','-v7');

%% 二、原几何上的真正形状反算，不拟合给定压力，不以几何闭合代替流动匹配
gas=createSurfaceFlowGas('thermallyPerfect',double(fs.R),1.4);
seed=struct('noseCenter',noseCenter,'referenceRadiusM',radius,'frame',frame, ...
    'standOffM',.15*radius,'K',.8/radius*eye(2));
radii=[.012 .022 .032];angles=[0 2*pi/3 4*pi/3];coordinates=zeros(9,2);k=0;
for r=radii
    for angle=angles,k=k+1;coordinates(k,:)=[r*cos(angle) r*sin(angle)];end
end
timer=tic;bodyFit=fitParabolicShockToBody(double(geometry.V),double(geometry.F),upstream,gas,seed,coordinates, ...
    struct('maximumIterations',8,'bodyMatchToleranceM',1e-5,'quadratureOrder',8));elapsed=toc(timer);
assert(norm(bodyFit.result.residualM)<norm(bodyFit.initial.residualM));
assert(~bodyFit.boundaryLayerEdgeAvailable&&~bodyFit.productionQualified&&bodyFit.newProductionSolvedAreaM2==0);

%% 三、输入和极限范围必须显式拒绝，不能让未解状态进入体表结果
rejections={};bad=shape;bad.K(1,1)=-1;
rejections{end+1}=refused(@()solveParabolicShockLayerColumn(upstream,gas,bad,q,eta),'ShockLayer:Shape');
rejections{end+1}=refused(@()solveParabolicShockLayerColumn(upstream,gas,shape,[0;0],eta),'ShockLayer:StagnationLimit');
rejections{end+1}=refused(@()solveParabolicShockLayerColumn(upstream,gas,shape,q,[1;.5;.6;0]),'ShockLayer:MassCoordinate');
report=struct('status','tested-internal-parabolic-shock-layer-columns-and-real-body-shape-iteration-prototype', ...
    'actualFaces',19752,'actualAngleDeg',20,'thermalModes',2,'sameStreamlineEntropyPreserved',true, ...
    'realBodyMatchingConverged',bodyFit.shapeMatchingConverged,'bodyMatchToleranceM',bodyFit.bodyMatchToleranceM, ...
    'initialBodyRMSM',norm(bodyFit.initial.residualM)/3,'finalBodyRMSM',norm(bodyFit.result.residualM)/3, ...
    'maximumBodyResidualM',max(abs(bodyFit.result.residualM)), ...
    'bodyFitEvaluations',bodyFit.evaluations,'bodyFitSeconds',elapsed,'rejections',{rejections}, ...
    'stagnationGlobalMatchingSolved',false,'boundaryLayerEdgeAvailable',false,'heatFluxSolved',false, ...
    'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'parabolic_shock_layer_evidence.mat'), ...
    'cases','bodyFit','seed','coordinates','geometry','fs','nose','report','-v7');
fid=fopen(fullfile(folder,'parabolic_shock_layer_report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(report);
end
function name=refused(action,expected)
name='';try,action();catch ex,name=ex.identifier;end
assert(strcmp(name,expected),'ShockLayer:TestRejection','实际拒绝%s，期望%s。',name,expected);
end
