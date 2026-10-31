function result=solveVariableMaterialConduction(h,Taw,time,material,initial,law)
%SOLVEVARIABLEMATERIALCONDUCTION 单层变物性材料的一维守恒导热。
% 储热：rho * 控制体厚度 * integral(cp(T),T_old,T_new)，不能用
% 当前 cp 乘温差代替跨温度区间的真实储热。密度保持常数。
% 导热：相邻节点的 Kirchhoff 势差 integral(k(T),T_i,T_j)/dx，
% 每个界面通量成对抵消。隐式 Euler + 精确三对角 Newton 雅可比，
% 用表内限幅和残差线搜索保证每次试探都有明确材料物性。
% 外侧对流/辐射，背面默认绝热；非绝热背面交给同一逐层守恒核。
% 无相变、烧蚀或横向导热。
if nargin<6, law=[]; end
backBoundary=wallBackBoundary(material);
if backBoundary.active
    result=calcWallConduction1D(h,Taw,time,material,initial,law);
    return;
end
nodes=12; if isfield(material,'thermalNodes'), nodes=material.thermalNodes; end
validateattributes(nodes,{'numeric'},{'real','finite','scalar','integer','>=',3,'<=',100});
validateattributes(time,{'numeric'},{'real','finite','vector','nonempty'});
time=time(:)';
if time(1)~=0 || any(diff(time)<=0)
    error('RocketAeroHeat:MaterialConduction:Time','时间须从零开始并严格递增。');
end
h=h(:); Taw=Taw(:); count=numel(h);
if numel(Taw)~=count, error('RocketAeroHeat:MaterialConduction:Size','h/Taw 尺寸不符。'); end
dx=material.thickness/(nodes-1); volume=dx*ones(1,nodes); volume([1 end])=dx/2;
mass=material.density*volume;
rad=material.emissivity*5.670374419e-8; Tr=material.radiationTemperature;
if isscalar(initial), old=initial*ones(count,nodes);
elseif isequal(size(initial),[count nodes]), old=initial;
else, error('RocketAeroHeat:MaterialConduction:Initial','初温须为标量或完整厚度剖面。'); end
validateattributes(old,{'numeric'},{'real','finite','positive'});
props=wallMaterialProperties(material,old);
range=props.temperatureRange;
surface=nan(count,numel(time)); surface(:,1)=old(:,1);
back=nan(count,numel(time)); back(:,1)=old(:,end);
for step=2:numel(time)
    dt=time(step)-time(step-1); value=old;
    for iteration=1:60
        [res,diagonal,below,above,q]=equations(value,old,dt);
        tolerance=1e-8+1e-10*max(1,abs(q));
        errorByFace=sum(abs(res),2);
        if all(errorByFace<=tolerance), break; end
        change=triSolve(diagonal,below,above,-res);
        fraction=ones(count,1);
        % 保持全部厚度节点在物性表内；不能把越界物性固定在端点继续算。
        for n=1:nodes
            rising=change(:,n)>0; falling=change(:,n)<0;
            fraction(rising)=min(fraction(rising), ...
                .99*(range(2)-value(rising,n))./change(rising,n));
            fraction(falling)=min(fraction(falling), ...
                .99*(range(1)-value(falling,n))./change(falling,n));
        end
        pending=errorByFace>tolerance;
        trial=value;
        for search=1:35
            candidate=value+fraction.*change;
            candidate(~pending,:)=trial(~pending,:);
            candidateRes=equations(candidate,old,dt);
            improved=sum(abs(candidateRes),2)<errorByFace | ...
                sum(abs(candidateRes),2)<=tolerance;
            accepted=pending & improved;
            trial(accepted,:)=candidate(accepted,:); pending(accepted)=false;
            if ~any(pending), break; end
            fraction(pending)=fraction(pending)/2;
        end
        if any(pending)
            error('RocketAeroHeat:MaterialConduction:RangeOrConvergence', ...
                '变物性导热无法在表内收敛；检查温度表范围、材料数据和时间步长。');
        end
        value=trial;
    end
    [res,~,~,~,q]=equations(value,old,dt);
    tolerance=1e-8+1e-10*max(1,abs(q));
    if any(sum(abs(res),2)>tolerance)
        error('RocketAeroHeat:MaterialConduction:Convergence','变物性导热残差未收敛。');
    end
    old=value; surface(:,step)=value(:,1);
    back(:,step)=value(:,end);
end
result=struct('method','材料比热精确积分储热 + Kirchhoff 导热通量 + 隐式 Euler/Newton', ...
    'scope','单层常密度材料，用户 cp(T)/k(T) 分段线性表；绝热背面，无相变/烧蚀，禁止表外外推', ...
    'surfaceTemperature',surface,'backTemperature',back,'finalDepthTemperature',old, ...
    'depth',linspace(0,material.thickness,nodes),'thermalNodes',nodes);

    function [res,d,lower,upper,q]=equations(T,previous,dt)
        p=wallMaterialProperties(material,T);
        storage=integrateWallMaterialProperty(material,'specificHeat',previous,T).*mass/dt;
        rightFlux=integrateWallMaterialProperty(material,'thermalConductivity', ...
            T(:,1:end-1),T(:,2:end))/dx;
        [localH,dh]=evaluateWallConvection(h,Taw,T(:,1),law);
        q=localH.*(Taw-T(:,1))-rad*(T(:,1).^4-Tr^4);
        res=storage;
        res(:,1:end-1)=res(:,1:end-1)-rightFlux;
        res(:,2:end)=res(:,2:end)+rightFlux;
        res(:,1)=res(:,1)-q;
        d=p.specificHeat.*mass/dt;
        links=[1 2*ones(1,nodes-2) 1];
        d=d+p.thermalConductivity.*links/dx;
        d(:,1)=d(:,1)+localH-dh.*(Taw-T(:,1))+4*rad*T(:,1).^3;
        lower=-p.thermalConductivity(:,1:end-1)/dx;
        upper=-p.thermalConductivity(:,2:end)/dx;
    end
end

function x=triSolve(d,lo,up,rhs)
% 批量 Thomas 法：变导热系数使上下对角不再是同一个常数。
[nf,n]=size(d); u=zeros(nf,n); b=u;
u(:,1)=up(:,1)./d(:,1); b(:,1)=rhs(:,1)./d(:,1);
for j=2:n
    pivot=d(:,j)-lo(:,j-1).*u(:,j-1);
    if j<n, u(:,j)=up(:,j)./pivot; end
    b(:,j)=(rhs(:,j)-lo(:,j-1).*b(:,j-1))./pivot;
end
x=b;
for j=n-1:-1:1, x(:,j)=b(:,j)-u(:,j).*x(:,j+1); end
end
