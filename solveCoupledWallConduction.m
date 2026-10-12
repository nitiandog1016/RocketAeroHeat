function result = solveCoupledWallConduction(evaluate,initialHeating,time,material,initial)
%SOLVECOUPLEDWALLCONDUCTION 热边界重算与一维材料导热的隐式耦合。
% 模块一：同一时间步在表面温度场与换热边界之间作 Picard 迭代。
% evaluate(Tw) 或 evaluate(Tw,time) 返回逐面元 h、Taw、eligibleFace；可同时返回算法
% 编号。每次厚度求解都从上一时刻的完整剖面开始，不重复累计时间。
% 模块二：厚度有限体积 + 隐式 Euler，外表面对流/辐射；背面默认
% 绝热，可由 material.backBoundary 指定内侧对流/辐射。
% 所有有效面元均解一维导热，不依赖固定 Bi 阈值选择简化模型。
% 模块三：同时核对温度收敛、重算热流一致性与离散能量残差。
% 关联式切换导致固定步长无可收敛解时，撤销该次试探并二分时间步；
% 输出保留所有实际内部时间站，热量积分也使用这些站，不能平均补值。
% 未知壁温输入 NaN。迭代中触发适用域失效的面元保守停止求解，
% 不用零热流补齐、不自动再附；失效面元整条输出历程置 NaN。
% 外流由 evaluate 调用者在给定时间建立，壁温迭代内固定；
% 此处不求横向导热、烧蚀或非定常外流。

%% 1. 检查时间、材料及初始热边界，分配完整厚度状态
validateattributes(initial,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(time,{'numeric'},{'real','finite','vector','nonempty','nonnegative'});
t=double(time(:)');
if t(1)~=0 || any(diff(t)<=0)
    error('RocketAeroHeat:CoupledWall:Time','时间须从 0 开始且严格递增。');
end
if ~isa(evaluate,'function_handle')
    error('RocketAeroHeat:CoupledWall:Evaluator','需要热边界重算函数句柄。');
end
timeDependent=nargin(evaluate)==2;
backBoundary=wallBackBoundary(material);
layeredMaterial=isfield(material,'layers') && ~isempty(material.layers);
if layeredMaterial
    material=normalizeLayeredWallMaterial(material,initial);
    materialGrid=buildLayeredWallGrid(material);
end
materialProperties=wallMaterialProperties(material,initial);
material.specificHeat=materialProperties.specificHeat;
material.thermalConductivity=materialProperties.thermalConductivity;
variableMaterial=isfield(material,'propertyTable') && ~isempty(material.propertyTable);
for name={'density','specificHeat','thickness','thermalConductivity','radiationTemperature'}
    validateattributes(material.(name{1}),{'numeric'},{'real','finite','scalar','positive'});
end
validateattributes(material.emissivity,{'numeric'},{'real','finite','scalar','>=',0,'<=',1});
nf=numel(initialHeating.eligibleFace); nt=numel(t);
initialHeating=normalizeHeating(initialHeating);
if nf*nt>1e7, error('RocketAeroHeat:CoupledWall:Memory','输出矩阵过大。'); end
nodes=12; if isfield(material,'thermalNodes'), nodes=material.thermalNodes; end
validateattributes(nodes,{'numeric'},{'real','finite','scalar','integer','>=',3,'<=',100});
requestedTime=t; pending=t(2:end); depth=zeros(size(pending)); t=t(1); nt=1;
subdivisions=0;
active=validHeating(initialHeating);
initialActive=active;
profile=initial*ones(nf,nodes);
temperature=nan(nf,nt); hHistory=temperature; tawHistory=temperature;
backHistory=temperature; finalInterface=struct('available',false);
qConv=temperature; qRad=temperature;
methodID=zeros(nf,nt,'uint8'); sourceID=methodID;
reasonHistory=methodID; metadata=cell(1,nt);
firstFailureTime=nan(nf,1); failureCode=zeros(nf,1,'uint8');
iterations=zeros(1,nt); energyError=nan(1,nt);
rad=material.emissivity*5.670374419e-8;
dx=material.thickness/(nodes-1);
capacity=material.density*material.specificHeat*dx*ones(1,nodes);
capacity([1 end])=capacity([1 end])/2;
lastHeating=initialHeating;
record(1,initial*ones(nf,1),initialHeating);

%% 2. 每个隐式时间步重算非局部热边界；通过双重收敛检查才提交
while ~isempty(pending)
    step=numel(t)+1;
    if nf*step>1e7, error('RocketAeroHeat:CoupledWall:Memory','时间自适应后输出矩阵过大。'); end
    targetTime=pending(1); dt=targetTime-t(end); old=profile;
    savedActive=active; savedFailureTime=firstFailureTime; savedCode=failureCode;
    temperature(:,step)=NaN; hHistory(:,step)=NaN; tawHistory(:,step)=NaN;
    qConv(:,step)=NaN; qRad(:,step)=NaN;
    methodID(:,step)=0; sourceID(:,step)=0; energyError(step)=NaN;
    guess=old(:,1); previousChange=Inf; accepted=false;
    changeHistory=nan(1,30);
    for iteration=1:30
        field=nan(nf,1); field(active)=guess(active);
        heat=normalizeHeating(boundaryAt(field,targetTime)); dropInvalid(heat,targetTime);
        if ~any(active)
            record(step,field,heat); accepted=true; lastHeating=heat; break;
        end
        state=calcWallConduction1D(heat.heatTransferCoefficient(active), ...
            heat.recoveryTemperature(active),[0 dt],material,old(active,:),[]);
        next=state.finalDepthTemperature;
        change=max(abs(next(:,1)-guess(active)));
        changeHistory(iteration)=change;
        if change<1e-6
            % 用候选新温度重新评估整条边界层，不能仅验证冻结 h 的方程。
            trial=nan(nf,1); trial(active)=next(:,1);
            check=normalizeHeating(boundaryAt(trial,targetTime));
            before=active; dropInvalid(check,targetTime);
            if ~isequal(before,active), continue; end
            actualQ=check.heatTransferCoefficient(active).* ...
                (check.recoveryTemperature(active)-next(:,1));
            backFlux=evaluateBackWallHeatFlux(next(:,end),backBoundary);
            outerRadiation=rad*(next(:,1).^4-material.radiationTemperature^4);
            flux=actualQ-outerRadiation-backFlux;
            storage=sum((next-old(active,:)).*capacity,2)/dt;
            if variableMaterial
                specificEnergy=integrateWallMaterialProperty(material, ...
                    'specificHeat',old(active,:),next);
                storage=sum(specificEnergy.*(capacity/material.specificHeat),2)/dt;
            end
            if layeredMaterial
                storage=sum(layeredWallEnergyIncrement(materialGrid,old(active,:),next),2)/dt;
            end
            errorFlux=storage-flux;
            scale=max(1,abs(actualQ)+abs(outerRadiation)+abs(backFlux));
            if all(abs(errorFlux)<=1e-8+1e-9*scale)
                profile(active,:)=next;
                if layeredMaterial
                    finalInterface=state.interfaceFinal;
                    for auditKey={'leftTemperature','rightTemperature','heatFluxFromLeft', ...
                            'heatFluxFromRight','contactHeatFlux','energyResidual'}
                        auditAll=nan(nf,size(finalInterface.(auditKey{1}),2));
                        auditAll(active,:)=finalInterface.(auditKey{1});
                        finalInterface.(auditKey{1})=auditAll;
                    end
                end
                record(step,trial,check);
                energyError(step)=max(abs(errorFlux));
                lastHeating=check; accepted=true; break;
            end
        end
        % 出现振荡时仅松弛试探温度，不改变旧时刻储热或重复推进时间。
        fraction=1;
        if iteration>2 && change>=0.95*previousChange, fraction=0.5; end
        guess(active)=guess(active)+fraction*(next(:,1)-guess(active));
        previousChange=change;
        % 六次迭代没有实质下降，尽早细分；无需等待固定上限反复振荡。
        if iteration>=12 && change>1e-4 && change>0.8*changeHistory(iteration-6)
            break;
        end
    end
    iterations(step)=iteration;
    if ~accepted
        % 回滚失败试探产生的失效掩码，再在同一旧剖面上缩短时间步。
        active=savedActive; firstFailureTime=savedFailureTime; failureCode=savedCode;
        if depth(1)>=10 || dt<=1e-8*max(1,requestedTime(end))
            error('RocketAeroHeat:CoupledWall:Convergence', ...
                't=%.9g s 的耦合经时间细分仍未收敛，不能接受冻结热流结果。',targetTime);
        end
        midpoint=t(end)+dt/2;
        pending=[midpoint pending];
        depth=[depth(1)+1 depth(1)+1 depth(2:end)];
        subdivisions=subdivisions+1;
        continue;
    end
    t(end+1)=targetTime; %#ok<AGROW>
    pending(1)=[]; depth(1)=[];
end

%% 3. 输出与现有壁温/热量 GUI 一致的数据结构以及耦合审计量
temperature(~active,:)=NaN; hHistory(~active,:)=NaN;
tawHistory(~active,:)=NaN; qConv(~active,:)=NaN; qRad(~active,:)=NaN;
methodID(~active,:)=0; sourceID(~active,:)=0;
bi=hHistory*material.thickness/(2*material.thermalConductivity);
if variableMaterial
    currentProperties=wallMaterialProperties(material,temperature);
    bi=hHistory*material.thickness./(2*currentProperties.thermalConductivity);
end
result=struct('method', ...
    '沿程热边界重算 + Picard 耦合 + 一维有限体积隐式导热/辐射', ...
    'scope',['固定外流/准流线，每步重积分边界层并检查分离适用域；' ...
    '常物性材料、绝热背面，无横向导热或烧蚀；试探状态失效保守停止，' ...
    '不预测分离区热流/再附；所有有效面元均求一维厚度导热'], ...
    'time',t,'material',material,'initialTemperature',initial, ...
    'temperature',temperature,'convectiveHeatFlux',qConv, ...
    'radiativeHeatFlux',qRad,'netHeatFlux',qConv-qRad, ...
    'heatTransferCoefficientHistory',hHistory, ...
    'recoveryTemperatureHistory',tawHistory, ...
    'Biot',initialHeating.heatTransferCoefficient(:)*material.thickness/(2*material.thermalConductivity), ...
    'maximumEstimatedBiot',max(bi,[],2),'wallTemperatureFeedback',true, ...
    'validFace',active,'lumpedFace',false(nf,1),'oneDimensionalFace',active);
result.oneDimensional=struct('finalDepthTemperature',profile(active,:), ...
    'depth',linspace(0,material.thickness,nodes),'thermalNodes',nodes);
backHistory(~active,:)=NaN;
result.oneDimensional.backTemperature=backHistory(active,:);
if layeredMaterial
    result.Biot(:)=NaN; result.maximumEstimatedBiot(:)=NaN;
    result.oneDimensional.depth=materialGrid.depth;
    result.oneDimensional.materialGrid=materialGrid;
    result.oneDimensional.backTemperature=backHistory(active,:);
    if finalInterface.available
        for auditKey={'leftTemperature','rightTemperature','heatFluxFromLeft', ...
                'heatFluxFromRight','contactHeatFlux','energyResidual'}
            finalInterface.(auditKey{1})=finalInterface.(auditKey{1})(active,:);
        end
    end
    result.oneDimensional.interfaceFinal=finalInterface;
end
result.coupling=struct('iterationsByStep',iterations, ...
    'requestedTime',requestedTime,'timeSubdivisions',subdivisions, ...
    'maximumFluxResidualByStep',energyError,'firstFailureTime',firstFailureTime, ...
    'failureReasonCode',failureCode,'lostFace',initialActive & ~active, ...
    'faceMethodIDHistory',methodID,'turbulentSourceIDHistory',sourceID, ...
    'finalHeating',lastHeating);
result.coupling.coverageReasonCodeHistory=reasonHistory;
result.coupling.boundaryMetadata=metadata;
result.materialMethod=materialProperties.method;
if timeDependent
    result.method=['逐时刻外流热边界 + ' result.method];
    result.scope=['时间相关准定常外流/准流线，每个试探时间站由回调独立建立；' ...
        '每步重积分边界层与壁温耦合，不求非定常流场；常物性材料、绝热背面，' ...
        '无横向导热或烧蚀；适用域失效保守停止，整条历程留空'];
end
if variableMaterial
    result.method=[result.method ' + 材料 cp(T) 积分储热/Kirchhoff 变导热系数'];
    result.scope=strrep(result.scope,'常物性材料','单层常密度、用户温度相关 cp/k 材料');
end
if layeredMaterial
    result.method=[result.method ' + 逐层 cp(T) 储热/Kirchhoff 通量/界面接触热阻'];
    result.materialMethod='逐层独立材料物性与界面热阻；不使用平均物性求解';
    result.scope=strrep(result.scope,'常物性材料','多层常密度材料、各层可给 cp(T)/k(T)');
end
result.summary=struct('calculatedFaces',nnz(active),'lumpedFaces',0, ...
    'oneDimensionalFaces',nnz(active),'lostFaces',nnz(initialActive & ~active), ...
    'maximumFinalTemperature',NaN);
if any(active), result.summary.maximumFinalTemperature=max(temperature(active,end)); end
result=completeWallEnergyBalance(result);

    function heat=normalizeHeating(heat)
        for key={'heatTransferCoefficient','recoveryTemperature','eligibleFace'}
            if numel(heat.(key{1}))~=nf
                error('RocketAeroHeat:CoupledWall:Size','热边界输出面元数不匹配。');
            end
            heat.(key{1})=heat.(key{1})(:);
        end
    end
    function valid=validHeating(heat)
        valid=logical(heat.eligibleFace(:)) & ...
            isfinite(heat.heatTransferCoefficient(:)) & heat.heatTransferCoefficient(:)>=0 & ...
            isfinite(heat.recoveryTemperature(:)) & heat.recoveryTemperature(:)>0;
    end
    function dropInvalid(heat,atTime)
        lost=active & ~validHeating(heat);
        firstFailureTime(lost)=atTime;
        failureCode(lost)=uint8(255); % 未提供覆盖分类的通用回调不能冒充代码 0“有解”。
        if isfield(heat,'coverageReasonCode'), failureCode(lost)=heat.coverageReasonCode(lost); end
        active(lost)=false;
    end
    function record(index,wall,heat)
        reasonHistory(:,index)=uint8(255);
        if isfield(heat,'coverageReasonCode')
            reasonHistory(:,index)=heat.coverageReasonCode(:);
        end
        metadata{index}=[];
        if isfield(heat,'boundaryMetadata'), metadata{index}=heat.boundaryMetadata; end
        temperature(active,index)=wall(active);
        backHistory(:,index)=NaN; backHistory(active,index)=profile(active,end);
        hHistory(active,index)=heat.heatTransferCoefficient(active);
        tawHistory(active,index)=heat.recoveryTemperature(active);
        qConv(active,index)=hHistory(active,index).*(tawHistory(active,index)-wall(active));
        qRad(active,index)=rad*(wall(active).^4-material.radiationTemperature^4);
        if isfield(heat,'faceMethodID'), methodID(active,index)=heat.faceMethodID(active); end
        if isfield(heat,'turbulentBoundaryLayer') && isstruct(heat.turbulentBoundaryLayer)
            sourceID(active,index)=heat.turbulentBoundaryLayer.sourceID(active);
        end
    end
    function heat=boundaryAt(wall,atTime)
        % 自动子步用自己的时刻评估外流，不能复用原时间段右端工况。
        if timeDependent, heat=evaluate(wall,atTime);
        else, heat=evaluate(wall); end
    end
end
