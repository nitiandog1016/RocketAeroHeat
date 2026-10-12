function result=solveLayeredWallConduction(h,Taw,time,material,initial,law)
%SOLVELAYEREDWALLCONDUCTION 多层、逐层变物性和接触热阻的一维导热。
% 模块一：逐层半控制体储热，理想界面共温度，热阻界面允许温度跳变。
% 模块二：隐式 Euler/Newton，层内 Kirchhoff 通量与界面热阻通量；
%           每层可以独立使用 cp(T)/k(T)，密度固定，禁止表外外推。
% 模块三：记录背面温度、最终厚度剖面和包含界面半控制体储热的
%           两侧界面热流核对。背面默认绝热，可指定对流/辐射；
%           没有横向导热、烧蚀或相变。
if nargin<6, law=[]; end
grid=buildLayeredWallGrid(material); n=grid.nodeCount;
backBoundary=wallBackBoundary(material);
h=h(:); Taw=Taw(:); nf=numel(h);
validateattributes(h,{'numeric'},{'real','finite','nonnegative'});
validateattributes(Taw,{'numeric'},{'real','finite','positive'});
validateattributes(time,{'numeric'},{'real','finite','vector','nonempty','nonnegative'});
time=time(:)';
if numel(Taw)~=nf || time(1)~=0 || any(diff(time)<=0)
    error('RocketAeroHeat:LayeredWall:Input','热边界尺寸或从零开始的递增时间无效。');
end
validateattributes(material.emissivity,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
validateattributes(material.radiationTemperature,{'numeric'},{'real','finite','scalar','positive'});
if isscalar(initial), old=initial*ones(nf,n);
elseif isequal(size(initial),[nf n]), old=initial;
else, error('RocketAeroHeat:LayeredWall:Initial','初温须为标量或完整多层温度剖面。'); end
validateattributes(old,{'numeric'},{'real','finite','positive'});
layeredWallEnergyIncrement(grid,old,old); % 检查每层初温适用域。
rad=material.emissivity*5.670374419e-8; Tr=material.radiationTemperature;
surface=nan(nf,numel(time)); back=surface; surface(:,1)=old(:,1); back(:,1)=old(:,end);
interfaces=struct('available',false);
for step=2:numel(time)
    dt=time(step)-time(step-1); value=old;
    for iteration=1:60
        [res,d,lower,upper,q]=equations(value,old,dt);
        tolerance=1e-8+1e-10*max(1,abs(q)); err=sum(abs(res),2);
        if all(err<=tolerance), break; end
        change=triSolve(d,lower,upper,-res); fraction=ones(nf,1);
        for k=1:n
            rising=change(:,k)>0; falling=change(:,k)<0;
            fraction(rising)=min(fraction(rising),.99*(grid.temperatureMaximum(k)-value(rising,k))./change(rising,k));
            fraction(falling)=min(fraction(falling),.99*(grid.temperatureMinimum(k)-value(falling,k))./change(falling,k));
        end
        pending=err>tolerance; trial=value;
        for search=1:35
            candidate=value+fraction.*change; candidate(~pending,:)=trial(~pending,:);
            nextRes=equations(candidate,old,dt); nextErr=sum(abs(nextRes),2);
            accepted=pending & (nextErr<err | nextErr<=tolerance);
            trial(accepted,:)=candidate(accepted,:); pending(accepted)=false;
            if ~any(pending), break; end
            fraction(pending)=fraction(pending)/2;
        end
        if any(pending)
            error('RocketAeroHeat:LayeredWall:Convergence','多层导热不能在各层物性表范围内收敛，请检查材料范围与时间步长。');
        end
        value=trial;
    end
    [res,~,~,~,q]=equations(value,old,dt);
    if any(sum(abs(res),2)>1e-8+1e-10*max(1,abs(q)))
        error('RocketAeroHeat:LayeredWall:Convergence','多层导热残差未收敛。');
    end
    interfaces=interfaceAudit(old,value,dt);
    old=value; surface(:,step)=value(:,1); back(:,step)=value(:,end);
end
result=struct('method','逐层比热积分储热/Kirchhoff 通量 + 界面接触热阻 + 隐式 Euler/Newton', ...
    'scope','多层常密度材料，可逐层给 cp(T)/k(T)；绝热背面，无相变/烧蚀/横向导热', ...
    'surfaceTemperature',surface,'backTemperature',back,'finalDepthTemperature',old, ...
    'depth',grid.depth,'thermalNodes',n,'materialGrid',grid,'interfaceFinal',interfaces);
result.backBoundary=backBoundary;
if backBoundary.active
    result.method=[result.method '；' backBoundary.method];
    result.scope=strrep(result.scope,'绝热背面','指定内侧对流/辐射');
    result.scope=[result.scope '；' backBoundary.scope];
end

    function [res,d,lo,up,q]=equations(T,previous,dt)
        [flux,left,right]=layeredWallFlux(grid,T);
        res=layeredWallEnergyIncrement(grid,previous,T)/dt;
        d=zeros(nf,n);
        for layerIndex=1:numel(grid.layers)
            ids=grid.layerNodes{layerIndex}; p=wallMaterialProperties(grid.layers{layerIndex},T(:,ids));
            d(:,ids)=d(:,ids)+p.specificHeat.*grid.layerMass{layerIndex}/dt;
        end
        res(:,1:end-1)=res(:,1:end-1)-flux; res(:,2:end)=res(:,2:end)+flux;
        d(:,1:end-1)=d(:,1:end-1)-left; d(:,2:end)=d(:,2:end)+right;
        lo=left; up=-right;
        [hc,dh]=evaluateWallConvection(h,Taw,T(:,1),law);
        q=hc.*(Taw-T(:,1))-rad*(T(:,1).^4-Tr^4);
        res(:,1)=res(:,1)-q;
        d(:,1)=d(:,1)+hc-dh.*(Taw-T(:,1))+4*rad*T(:,1).^3;
        [backFlux,backSlope]=evaluateBackWallHeatFlux(T(:,end),backBoundary);
        res(:,end)=res(:,end)+backFlux;
        d(:,end)=d(:,end)+backSlope;
        q=abs(q)+abs(backFlux); % 残差尺度使用两侧交换量，避免近稳态相消。
        if any(~isfinite(res),'all') || any(~isfinite(d),'all') || any(~isfinite(q))
            error('RocketAeroHeat:LayeredWall:Numerical', ...
                '多层导热残差或系数溢出，请检查材料、接触热阻及热边界量级。');
        end
    end

    function audit=interfaceAudit(previous,T,dt)
        flux=layeredWallFlux(grid,T); ni=numel(grid.contactResistance);
        fromLeft=zeros(nf,ni); fromRight=fromLeft; contact=nan(nf,ni);
        for j=1:ni
            il=grid.interfaceNodes(j,1); ir=grid.interfaceNodes(j,2);
            eLeft=integrateWallMaterialProperty(grid.layers{j},'specificHeat',previous(:,il),T(:,il))*grid.layerMass{j}(end)/dt;
            eRight=integrateWallMaterialProperty(grid.layers{j+1},'specificHeat',previous(:,ir),T(:,ir))*grid.layerMass{j+1}(1)/dt;
            fromLeft(:,j)=-flux(:,il-1)-eLeft;
            fromRight(:,j)=-flux(:,ir)+eRight;
            if grid.contactResistance(j)>0
                contact(:,j)=(T(:,il)-T(:,ir))/grid.contactResistance(j);
            end
        end
        audit=struct('available',true,'leftTemperature',T(:,grid.interfaceNodes(:,1)), ...
            'rightTemperature',T(:,grid.interfaceNodes(:,2)), ...
            'heatFluxFromLeft',fromLeft,'heatFluxFromRight',fromRight, ...
            'contactHeatFlux',contact,'energyResidual',fromLeft-fromRight, ...
            'method','界面热流扣除各自半控制体的储热，正号表示由外向内');
    end
end

function x=triSolve(d,lo,up,rhs)
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
