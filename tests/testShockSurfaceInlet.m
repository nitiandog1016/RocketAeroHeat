function report=testShockSurfaceInlet(folder,sourceFile)
%TESTSHOCKSURFACEINLET 共同三维激波的守恒、物理梯度及真实20度上游入口校核。
% 模块一：从真实箭体快照只读V/F和自由来流；旧外缘/热流不作为入口。
% 模块二：同一非轴对称候选曲面采样，测试热力模式、坐标变换与刚体旋转。
% 模块三：有限差分只校核隐式梯度；所有形状/壁面/热流资格仍为未解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupEngineeringShockLayer;
assert(~exist(folder,'dir'),'ShockLayer:Overwrite','验证必须写新目录。');
old=load(sourceFile,'result','fs');mesh=old.result;fs=old.fs;
assert(size(mesh.F,1)==19752&&fs.alpha==20&&fs.Mach==5&&fs.altitude==30000);
geometry=struct('V',double(mesh.V),'F',double(mesh.F),'area',double(mesh.area));clear mesh old;
upstream=struct('pressure',double(fs.pressure),'temperature',double(fs.temperature), ...
    'density',double(fs.density),'velocityVector',double(fs.flowVelocityVector(:)), ...
    'method','真实19752面箭体M5/30km/alpha20/beta0自由来流；未读旧外缘作为入口');
K=[4 1.1;1.1 2.8];tip=[max(geometry.V(:,1))+.02;0;0];
[y,z]=ndgrid([-.06 0 .06],[-.05 0 .05]);coordinates=[y(:) z(:)];
cases=cell(2,1);differences=cell(2,1);
for j=1:2
    modes={'caloricallyPerfect','thermallyPerfect'};gas=createSurfaceFlowGas(modes{j},double(fs.R),1.4);
    p=struct('sample',@(q)surface(q,K,tip,eye(2),eye(3)), ...
        'upstream',upstream,'gas',gas,'normalOrientation',-1,'sources',struct( ...
        'geometry','一张给定非轴对称抛物候选激波；K=[4,1.1;1.1,2.8]/m；未拟合真实箭体', ...
        'upstream',upstream.method,'jumpRelations','三维定常Rankine-Hugoniot；理想固定组分气体，同cp积分'));
    inlet=sampleShockSurfaceInlet(p,coordinates);
    assert(all(cellfun(@(s)s.available,inlet.states))&&~inlet.bodyMatched&&~inlet.boundaryLayerEdgeAvailable);
    assert(max(cellfun(@(s)abs(s.residual.totalEnthalpy),inlet.states))<1e-10);
    entropy=cellfun(@(s)s.entropy,inlet.states);assert(max(entropy)-min(entropy)>1);
    % 对真实20度上游检查非零横向速度；禁止把速度投影回自由流切向。
    assert(any(cellfun(@(s)norm(s.partialDerivatives.velocityVector,'fro')>1,inlet.states)));
    q=[.017;-.023];h=[2e-4 1e-4 5e-5];checks=cell(3,1);
    center=shockSurfaceJumpState(upstream,gas,p.sample(q),-1);
    for k=1:3
        numerical=zeros(5,2);
        for d=1:2
            shift=zeros(2,1);shift(d)=h(k);
            plus=shockSurfaceJumpState(upstream,gas,p.sample(q+shift),-1);
            minus=shockSurfaceJumpState(upstream,gas,p.sample(q-shift),-1);
            numerical(:,d)=(values(plus)-values(minus))/(2*h(k));
        end
        exact=[center.partialDerivatives.pressure;center.partialDerivatives.temperature; ...
            center.partialDerivatives.density;center.partialDerivatives.speed;center.partialDerivatives.entropy];
        checks{k}=struct('stepM',h(k),'numerical',numerical,'exact',exact, ...
            'relativeDifference',vecnorm(numerical-exact,2,2)./max(1,vecnorm(exact,2,2)));
    end
    first=checks{1}.relativeDifference;last=checks{3}.relativeDifference;
    assert(all(first./last>10)&&max(last)<2e-6,'ShockLayer:TestDerivative','两向导数未随网格收敛。');
    % 非正交/有量纲坐标变换及刚体旋转，物理梯度应随旋转而非标签缩放。
    map=[2 .3;-.2 .7];angle=.63;rotation=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
    transformedUpstream=upstream;transformedUpstream.velocityVector=rotation*upstream.velocityVector;
    s=shockSurfaceJumpState(transformedUpstream,gas,surface(map\q,K,tip,map,rotation),-1);
    assert(max(abs(values(s)-values(center))./max(1,abs(values(center))))<1e-10);
    assert(norm(s.velocityVector-rotation*center.velocityVector)<1e-9);
    for name=fieldnames(center.physicalGradients)'
        key=name{1};assert(norm(s.physicalGradients.(key)-rotation*center.physicalGradients.(key)) ...
            <1e-8*max(1,norm(center.physicalGradients.(key))));
    end
    cases{j}=inlet;differences{j}=checks;
end

%% 四、标准平面激波/近弱分支供独立闭式和Brent审核，不只核对自有残差
planeCases=cell(4,1);gas=createSurfaceFlowGas('caloricallyPerfect',287.05,1.4);
up=struct('pressure',1200,'temperature',220,'density',1200/(287.05*220), ...
    'velocityVector',[-5*sqrt(1.4*287.05*220);0;0],'method','标准M5平面激波解析校核上游');
for j=1:4
    normalMach=[1.01 1.5 3 5];beta=asin(normalMach(j)/5);
    n=[-sin(beta);cos(beta);0];a=[cos(beta);sin(beta);0];b=cross(n,a);
    sample=struct('position',zeros(3,1),'firstDerivatives',[a b],'secondDerivatives',zeros(3,3));
    state=shockSurfaceJumpState(up,gas,sample,1);
    r=2.4*normalMach(j)^2/(.4*normalMach(j)^2+2);
    assert(abs(state.densityRatio-r)<1e-9);
    planeCases{j}=struct('upstream',up,'geometry',sample,'gas',gas,'state',state);
end

%% 五、状态不一致、翻向、无激波、未知来源和退化几何均须拒绝
rejections={};gas=cases{2}.gas;sample=surface([0;0],K,tip,eye(2),eye(3));
bad=upstream;bad.density=bad.density*1.01;
rejections{end+1}=refused(@()shockSurfaceJumpState(bad,gas,sample,-1),'ShockLayer:EOS');
rejections{end+1}=refused(@()shockSurfaceJumpState(upstream,gas,sample,1),'ShockLayer:NormalMach');
bad=upstream;bad.velocityVector=bad.velocityVector*.1;
rejections{end+1}=refused(@()shockSurfaceJumpState(bad,gas,sample,-1),'ShockLayer:NormalMach');
badGeometry=sample;badGeometry.firstDerivatives(:,2)=badGeometry.firstDerivatives(:,1);
rejections{end+1}=refused(@()shockSurfaceJumpState(upstream,gas,badGeometry,-1),'ShockLayer:Geometry');
bad=p;bad.sources=rmfield(bad.sources,'geometry');
rejections{end+1}=refused(@()sampleShockSurfaceInlet(bad,coordinates),'ShockLayer:Source');
report=struct('status','tested-common-3D-candidate-shock-jumps-and-inlet-gradients-not-body-matched', ...
    'actualGeometryFaces',19752,'actualFreestreamAngleDeg',20,'commonSurfacePointsPerMode',9, ...
    'planeCases',4,'rejections',{rejections},'sameShockGeometryForAllSamples',true, ...
    'oldExternalStateUsedAsInlet',false,'shapeSolved',false,'bodyMatched',false, ...
    'boundaryLayerEdgeAvailable',false,'pressureOnBodySolved',false,'heatFluxSolved',false, ...
    'experimentalValidation',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'shock_surface_inlet_evidence.mat'), ...
    'cases','differences','planeCases','report','geometry','fs','K','tip','coordinates','-v7');
fid=fopen(fullfile(folder,'shock_surface_inlet_report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(report);
end

function sample=surface(q,K,tip,map,rotation)
% 模块用途：同一非轴对称候选曲面；仅几何给定，流动状态由跃迁方程求解。
v=map*q;position=tip+[-.5*v'*K*v;v];base=[-(K*v)';eye(2)];
A=base*map;H=-map'*K*map;B=[H(1,1) H(1,2) H(2,2);zeros(2,3)];
sample=struct('position',rotation*position,'firstDerivatives',rotation*A,'secondDerivatives',rotation*B);
end
function v=values(s)
v=[s.pressure;s.temperature;s.density;s.speed;s.entropy];
end
function name=refused(action,expected)
name='';try,action();catch ex,name=ex.identifier;end
assert(strcmp(name,expected),'ShockLayer:TestRejection','实际拒绝%s，期望%s。',name,expected);
end
