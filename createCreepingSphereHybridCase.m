function cfg=createCreepingSphereHybridCase(mesh,c,wallTemperatureK,options)
%CREATECREEPINGSPHEREHYBRIDCASE 均匀低速来流中等温实心球的黏性外流与热浴传热。
% 仅适用于恒物性、无体力的小 Reynolds 数球体。动量取 Stokes 极限，
% 能量保留有限 Pe 的对流扩散，外侧热浴恒温。接口内外都有法向流动。
if nargin<4,options=struct;end
defaults=struct('cutRatio',2,'bathRatio',8,'innerCells',24,'outerCells',24,'angularCells',32,'saveFolder','');
assert(isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'RAH:Sphere:Options','未知网格选项。');
for key=fieldnames(defaults)'
    if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end
end
c=conditionForSphere(c);
for key={'innerCells','outerCells','angularCells'}
    validateattributes(options.(key{1}),{'numeric'},{'scalar','integer','>=',4});
end
validateattributes(options.cutRatio,{'numeric'},{'real','finite','scalar','>',1});
validateattributes(options.bathRatio,{'numeric'},{'real','finite','scalar','>',options.cutRatio});
n=size(mesh.F,1);
if isscalar(wallTemperatureK),wallTemperatureK=repmat(wallTemperatureK,n,1);end
validateattributes(wallTemperatureK,{'numeric'},{'vector','numel',n,'positive','real','finite'});
assert(all(wallTemperatureK==wallTemperatureK(1)),'RAH:Sphere:Wall','本方法要求等温球壁。');
% 热浴半径是热边界的位置。切面和离散参数只控制域分解与网格。
c.thermalBathRadiusM=options.bathRatio*c.radiusM;
cfg=struct('analysisMode','hybridSurface','hybridSurface',struct('condition',c, ...
    'region',repmat("sphere-wall",n,1),'wallTemperatureK',wallTemperatureK(:), ...
    'engineeringSolver',@outerFlow,'localSolver',@nearWall));

    function wall=outerFlow(g,condition,tw)
        % 一、外域衰减模态：由切面速度确定远场扰动，再向内域返回应力条件。
        condition=conditionForSphere(condition);
        assert(all(tw==tw(1)),'RAH:Sphere:Wall','本方法要求等温球壁。');
        assert(all(g.region=="sphere-wall"),'RAH:Sphere:Geometry','当前适配器求解单一球壁。');
        validateattributes(condition.thermalBathRadiusM,{'numeric'},{'finite','scalar','real','>',condition.radiusM*options.cutRatio});
        unit=(g.V-condition.centerM)/condition.radiusM;
        assert(all(abs(vecnorm(unit,2,2)-1)<1e-8),'RAH:Sphere:Geometry','顶点须在指定球面上。');
        edges=[g.F(:,[1 2]);g.F(:,[2 3]);g.F(:,[3 1])];
        [~,~,group]=unique(sort(edges,2),'rows');
        assert(all(accumarray(group,1)==2) && all(accumarray(group,sign(edges(:,2)-edges(:,1)))==0), ...
            'RAH:Sphere:Geometry','球壁须为方向一致的闭合三角网格。');
        a=unit(g.F(:,1),:);b=unit(g.F(:,2),:);d=unit(g.F(:,3),:);
        angles=2*atan2(dot(a,cross(b,d,2),2),1+dot(a,b,2)+dot(a,d,2)+dot(b,d,2));
        assert(all(angles>0) && abs(sum(angles)-4*pi)<1e-8,'RAH:Sphere:Geometry','球壁方向或球面覆盖无效。');
        s=options.cutRatio;
        basis=[s 1/s;1 -1/s^2];uniform=[s^2;2*s];
        response=[0 2/s^3;0 -6/s^4]/basis;
        wall=table(zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
            'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
        wall.Properties.UserData=struct('basis',basis,'uniform',uniform,'response',response, ...
            'offset',[2;0]-response*uniform,'cutRatio',s,'condition',condition,'wallTemperatureK',tw(1));
    end

    function local=nearWall(g,condition,tw,regionIDs,targets,engineering)
        % 二、用 MATLAB 配点法求近区 Stokes 流函数。速度、剪切应力与法向应力在切面连续。
        cut=engineering.Properties.UserData;s=cut.cutRatio;
        [condition,speed,re,pe]=conditionForSphere(condition);
        assert(isequaln(condition,cut.condition) && all(tw==cut.wallTemperatureK), ...
            'RAH:Sphere:Interface','切面来源与当前工况或壁温不匹配。');
        bc=@(left,right)[left(1:2);right(3:4)-cut.response*right(1:2)-cut.offset];
        initial=bvpinit(linspace(1,s,12),@(r)[0;0;0;0]);
        folder=options.saveFolder;
        if ~isempty(folder)
            assert(isfolder(folder)&&~isfile(fullfile(folder,'initial.mat'))&&~isfile(fullfile(folder,'terminal.mat')), ...
                'RAH:Sphere:Output','保存目录须存在且不得覆盖原运行。');
            save(fullfile(folder,'initial.mat'),'initial','cut','condition','options','tw');
        end
        flow=bvp4c(@stokesODE,bc,initial,bvpset('RelTol',1e-10,'AbsTol',1e-12,'NMax',4000));
        inner=deval(flow,s);coef=cut.basis\(inner(1:2)-cut.uniform);
        outer=exterior(s,coef);
        pts=linspace(1,s,1001);[state,derivative]=deval(flow,pts);rhs=stokesODE(pts,state);
        residual=max(abs(derivative-rhs)./(1+abs(rhs)),[],'all');
        stressJump=max(abs(traction(s,inner)-traction(s,outer)));
        traceJump=max(abs(inner(1:2)-outer(1:2)));
        % 三、近区和外侧热浴分别使用各自的速度，联立求解切面两侧的热流与界面温度。
        thermalGrid=options;thermalGrid.bathRatio=condition.thermalBathRadiusM/condition.radiusM;
        thermal=solveSphereHeatBath(flow,coef,pe,thermalGrid);
        thermal.temperatureK=condition.ambientTemperatureK+(tw(1)-condition.ambientTemperatureK)*thermal.normalizedTemperature;
        fluxScale=condition.conductivityWmK*(tw(1)-condition.ambientTemperatureK)/condition.radiusM;
        wallQ=-fluxScale*thermal.wallFluxDensity;
        wallState=deval(flow,1);
        dpScale=condition.viscosityPaS*speed/condition.radiusM;
        wallPressureMode=dpScale*pressureMode(1,wallState);
        motion=[-1 0 0];if speed>0,motion=condition.velocityMps/speed;end
        centers=(g.V(g.F(:,1),:)+g.V(g.F(:,2),:)+g.V(g.F(:,3),:))/3-condition.centerM;
        mu=(centers*motion')./vecnorm(centers,2,2);
        % 理论范围端点可因点积舍入偏离 1 个 ulp；按最外控制带处理两个极点。
        band=discretize(mu,[-Inf;thermal.muEdges(2:end-1);Inf]);
        wallScale=2*pi*condition.radiusM^2;
        mscale=wallScale*condition.densityKgM3*speed;
        mInner=mscale*inner(1)*thermal.muMoments;mOuter=mscale*outer(1)*thermal.muMoments;
        massJump=max(abs(mInner-mOuter))/max(mscale,realmin);
        % 按物理半球积分。即使某个角带跨赤道，两个方向也分别计量。
        south=diff(min(thermal.muEdges,0).^2)/2;north=diff(max(thermal.muEdges,0).^2)/2;
        massIn=-mscale*inner(1)*sum(south);massOut=mscale*inner(1)*sum(north);
        enerScale=2*pi*condition.conductivityWmK*(tw(1)-condition.ambientTemperatureK)*condition.radiusM;
        advected=enerScale*pe*inner(1)*thermal.muMoments.*thermal.cutNormalizedTemperature;
        audit=struct('reynoldsDiameter',re,'pecletRadius',pe,'flowResidual',residual, ...
            'interfaceVelocityModeJump',traceJump,'interfaceStressModeJump',stressJump, ...
            'interfaceMassInKgS',massIn,'interfaceMassOutKgS',massOut,'interfaceScaledMassJump',massJump, ...
            'interfaceNetMassOutKgS',sum(mInner),'interfaceHeatFluxJump',thermal.interfaceDefect, ...
            'interfaceAdvectedHeatAbsoluteW',sum(abs(enerScale*pe*inner(1)*thermal.cutNormalizedTemperature).*(north-south)), ...
            'interfaceAdvectedHeatNetOutW',sum(advected), ...
            'interfaceConductedHeatNetOutW',sum(enerScale*thermal.innerCutFlux-advected), ...
            'energyBalanceRelative',thermal.balance,'thermalResidual',thermal.residual, ...
            'curvedWallPowerW',-enerScale*sum(thermal.wallFlux),'surfacePressureModePa',wallPressureMode, ...
            'cutCountedAsWall',false,'flowNodes',numel(flow.x),'thermalUnknowns',numel(thermal.temperatureK));
        ok=residual<1e-8 && traceJump<1e-8 && stressJump<1e-8 && massJump<1e-8 && ...
            (speed==0 || (massIn>0 && massOut>0)) && thermal.interfaceDefect<1e-8 && ...
            thermal.balance<1e-9 && thermal.residual<1e-10 && thermal.bounded;
        if ~isempty(folder)
            save(fullfile(folder,'terminal.mat'),'flow','coef','thermal','audit','ok','wallQ','mInner','mOuter','advected');
        end
        local=struct('accepted',ok,'reason',"球外穿流接口或离散收支未通过",'wall',table(),'diagnostics',audit);
        if ~ok,return;end
        assert(all(ismember(targets,regionIDs)));
        ids=targets+1;
        pressure=condition.ambientPressurePa+wallPressureMode*mu(ids);
        local.wall=table(targets,pressure,wallQ(band(ids)), ...
            repmat("壁压：球外 Stokes 配点解；入壁热流：有限 Pe 对流扩散有限体积解",numel(ids),1), ...
            'VariableNames',engineering.Properties.VariableNames);
        local.wall.pressureMethod=repmat("球外 Stokes 配点解，由动量方程恢复压力",numel(ids),1);
        local.wall.heatFluxMethod=repmat("有限 Pe 球坐标对流扩散有限体积解，取球壁 Fourier 通量",numel(ids),1);
        adjacent=thermal.temperatureK(2,:)';
        local.wall.nearWallCellTemperatureK=adjacent(band(ids));
        local.wall.nearWallTemperatureMethod=repmat("球坐标对流扩散有限体积解，取所属角带首层流体单元温度",numel(ids),1);
        local.reason="内外域速度、应力、双向质量交换与热通量已对照";
    end
end

function [c,speed,re,pe]=conditionForSphere(c)
% 首次配置和每次运行均以当次输入计算无量纲参数，工况编辑后重新检查。
for key={'radiusM','densityKgM3','viscosityPaS','specificHeatJkgK','conductivityWmK','ambientTemperatureK','ambientPressurePa'}
    validateattributes(c.(key{1}),{'numeric'},{'real','finite','scalar','positive'});
end
for key={'centerM','velocityMps'}
    validateattributes(c.(key{1}),{'numeric'},{'real','finite','vector','numel',3});
    c.(key{1})=reshape(c.(key{1}),1,3);
end
speed=norm(c.velocityMps);
re=2*c.densityKgM3*speed*c.radiusM/c.viscosityPaS;
pe=c.densityKgM3*c.specificHeatJkgK*speed*c.radiusM/c.conductivityWmK;
assert(re<=0.01 && pe<=1,'RAH:Sphere:Regime','超出本低速、弱对流球体适配器的范围。');
end

function rate=stokesODE(r,y)
% 无旋流分量的轴对称 Stokes 流函数满足 (D²−2/r²)² f=0。
rate=[y(2:4,:);4*y(3,:)./r.^2-8*y(2,:)./r.^3+8*y(1,:)./r.^4];
end

function y=exterior(s,c)
% Stokes 远域的均匀流、Stokeslet 和偶极模态；系数由切面解返回。
y=[s.^2+c(1)*s+c(2)./s;2*s+c(1)-c(2)./s.^2; ...
    2+2*c(2)./s.^3;-6*c(2)./s.^4];
end

function p=pressureMode(s,y)
p=(y(4)-2*y(2)/s^2+4*y(1)/s^3)/2;
end

function stress=traction(s,y)
stress=[-pressureMode(s,y)+2*(y(2)/s^2-2*y(1)/s^3);-y(3)/(2*s)+y(2)/s^2-y(1)/s^3];
end

function heat=solveSphereHeatBath(flow,coef,pe,opt)
% 四、轴对称球壳的稳态有限体积能量方程。长度以球半径计，温度以壁浴温差计。
% 交界温度是一个零容积未知量。内外控制体分别贡献通量，界面方程要求二者相等。
cut=opt.cutRatio;
ri=exp(linspace(0,log(cut),opt.innerCells+1))';ri([1 end])=[1;cut];
ro=exp(linspace(log(cut),log(opt.bathRatio),opt.outerCells+1))';ro([1 end])=[cut;opt.bathRatio];
radius=[1;(ri(1:end-1)+ri(2:end))/2;cut;(ro(1:end-1)+ro(2:end))/2;opt.bathRatio];
faceRadius=[ri;ro];width=[0;diff(ri);0;diff(ro);0];
muEdges=linspace(-1,1,opt.angularCells+1)';dmu=diff(muEdges);
muMoments=diff(muEdges.^2)/2;
innerValues=deval(flow,ri');outerValues=exterior(ro',coef);
fi=innerValues(1,:)';fo=outerValues(1,:)';
faceF=[fi;fo];deltaF=[0;diff(fi);0;diff(fo);0];
index=reshape(1:numel(radius)*opt.angularCells,numel(radius),opt.angularCells);

% 径向导热系数使用半径倒数之差，静止流体的球面导热阻力可逐层相加。
% 对流通量由流函数的带端差积分，温度以同一界面值传递。
radialResistance=1./radius(1:end-1)-1./radius(2:end);
rightWeight=(1./radius(1:end-1)-1./faceRadius)./radialResistance;
conductance=(1./radialResistance)*dmu';
transport=pe*faceF*muMoments';
left=conductance+transport.*(1-rightWeight);right=-conductance+transport.*rightWeight;
% 余弦方向正向为极角减小。两端轴上的面积为零，无需虚构极点温度。
angularFactor=1-muEdges(2:end-1)'.^2;
angularDiffusion=width*angularFactor/(2/opt.angularCells);
angularTransport=(pe/2)*deltaF*angularFactor;
leftAngular=angularDiffusion+angularTransport/2;
rightAngular=-angularDiffusion+angularTransport/2;
lo=[reshape(index(1:end-1,:),[],1);reshape(index(:,1:end-1),[],1)];
hi=[reshape(index(2:end,:),[],1);reshape(index(:,2:end),[],1)];
cl=[left(:);leftAngular(:)];cr=[right(:);rightAngular(:)];
n=numel(index);
matrix=sparse([lo;lo;hi;hi],[lo;hi;lo;hi],[cl;cr;-cl;-cr],n,n);
fixed=[index(1,:) index(end,:)]';matrix(fixed,:)=0;
matrix=matrix+sparse(fixed,fixed,ones(size(fixed)),n,n);
rhs=zeros(n,1);rhs(index(1,:))=1;
theta=matrix\rhs;
residual=max(abs(matrix*theta-rhs)./max(1,abs(matrix)*abs(theta)+abs(rhs)));
t=reshape(theta,size(index));
flux=left.*t(1:end-1,:)+right.*t(2:end,:);
cutRow=opt.innerCells+2;wall=flux(1,:)';
heat=struct('radiusRatio',radius,'muEdges',muEdges,'muMoments',muMoments,'normalizedTemperature',t, ...
    'wallFlux',wall,'wallFluxDensity',wall./dmu,'bathFlux',flux(end,:)', ...
    'cutNormalizedTemperature',t(cutRow,:)','innerCutFlux',flux(cutRow-1,:)','outerCutFlux',flux(cutRow,:)', ...
    'interfaceDefect',max(abs(flux(cutRow-1,:)-flux(cutRow,:)))/max(abs(wall)), ...
    'balance',abs(sum(wall)-sum(flux(end,:)))/sum(abs(wall)), ...
    'residual',residual,'bounded',all(isfinite(theta) & theta>=-1e-10 & theta<=1+1e-10));
end
