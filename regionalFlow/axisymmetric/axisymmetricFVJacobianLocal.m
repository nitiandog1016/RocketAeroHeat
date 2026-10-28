function [J,report,pattern]=axisymmetricFVJacobianLocal(Q,mesh,gas,boundary,scale,options,pattern)
%AXISYMMETRICFVJACOBIANLOCAL 使用显式局部物理扰动的完整有色差分Jacobian。
% 模块零：仅改变扰动参考，原空间残差、稀疏支持、拆列与中央/单边回退逻辑保持。
% 原axisymmetricFVJacobian文件保持原字节；新增策略必须用独立工况/源码身份测试。
% 模块一：U=Q./scale，J=d[(原始rate)./scale]/d[U(:)]，沿MATLAB分量分块布局。
% 模块二：WLS边界迹、共享面双侧梯度/非正交修正和柱坐标源构造保守稀疏支持。
% 模块三：实际核对同色列不共享残差行；四/五守恒分量分别扰动，不能同时混色。
% 模块四：中央差分失败先缩扰动，再显式单边或拆分单列；绝不以零导数替代失败。
% 模块五：计入基准、正/负方向及失败的每一次残差调用；不重算dt或全局尺度。
% 参数/几何合同错误抛异常；耗尽有效扰动则J=[]且report.success=false供上层拒绝。
% 本模块不推进物理时间、不建立湍流/实验/真实箭体适用资格，新增物理解面积为0。

%% 一、固定参数、布局与差分选项；只对同一空间算子求导
if nargin<6||isempty(options),options=struct;end
if nargin<7,pattern=[];end
defaults=struct('differenceStep',2e-6,'differenceScheme','central','maxPerturbationRetries',8);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'AxisFV:JacobianOptions','Jacobian选项只允许differenceStep/differenceScheme/maxPerturbationRetries。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.differenceStep,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(options.maxPerturbationRetries,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
assert((ischar(options.differenceScheme)&&isrow(options.differenceScheme))|| ...
    (isstring(options.differenceScheme)&&isscalar(options.differenceScheme)&&~ismissing(options.differenceScheme)), ...
    'AxisFV:JacobianOptions','differenceScheme须为central或forward。');
options.differenceScheme=char(options.differenceScheme);
assert(ismember(options.differenceScheme,{'central','forward'}), ...
    'AxisFV:JacobianOptions','differenceScheme只接受central或forward。');
nq=size(Q,2);assert(ismember(nq,[4 5]),'AxisFV:JacobianState','需要四量NS或五量SA完整状态。');
validateattributes(scale,{'numeric'},{'real','finite','positive','numel',nq});scale=reshape(double(scale),1,nq);
assert(isstruct(mesh)&&isscalar(mesh)&&isfield(mesh,'kind')&& ...
    strcmp(mesh.kind,'axisymmetric-body-fitted-fv-v1'),'AxisFV:JacobianGeometry','需要正式贴体FV几何。');
nc=mesh.cellCount;
validateattributes(Q,{'numeric'},{'real','finite','size',[nc nq]});Q=double(Q);
gas=axisymmetricFVGas(gas);U=Q./scale;
[localReference,scaleAudit]=axisymmetricJacobianLocalReference(Q,gas,scale,options.differenceStep);
assert(all(isfinite(U),'all'),'AxisFV:JacobianScale','归一化守恒量必须有限。');

%% 二、保守支持：梯度依赖H，面依赖B''*(I+H)，残差依赖B*面依赖与源I+H
bf=mesh.boundaryFaces(:);nb=numel(bf);owner=mesh.faces.owner(:);
assert(size(mesh.gradientX,1)==nc&&isequal(size(mesh.gradientX),[nc nc+nb])&& ...
    isequal(size(mesh.gradientR),[nc nc+nb])&& ...
    isequal(size(mesh.cellFaceMultiplicity),[nc numel(owner)]), ...
    'AxisFV:JacobianGeometry','梯度样本与边界面顺序或面拓扑尺寸不一致。');
T=sparse((1:nb)',owner(bf),ones(nb,1),nb,nc);
Z=[speye(nc);T];
H=spones((abs(mesh.gradientX)+abs(mesh.gradientR))*Z);
B=spones(mesh.cellFaceMultiplicity);K=speye(nc)+H;
cellPattern=spones(B*(B'*K)+K);
binding=geometryBinding(mesh);
reused=~isempty(pattern);
if reused
    assert(isstruct(pattern)&&isscalar(pattern)&&all(isfield(pattern, ...
        {'kind','cellPattern','colors','geometryBinding'}))&& ...
        strcmp(pattern.kind,'axisymmetric-fv-colored-jacobian-pattern-v1'), ...
        'AxisFV:JacobianPattern','无效的Jacobian缓存结构。');
    assert(isequaln(pattern.geometryBinding,binding),'AxisFV:JacobianPatternGeometry', ...
        '复用的Jacobian缓存不属于当前几何/梯度/周期拓扑。');
    assert(issparse(pattern.cellPattern)&&isequal(pattern.cellPattern,cellPattern), ...
        'AxisFV:JacobianPatternSupport','缓存稀疏支持与当前保守依赖不一致。');
    colors=pattern.colors(:);
else
    % 列冲突由共享残差行判定，不能直接对cellPattern的一阶邻接着色。
    conflict=spones(cellPattern'*cellPattern);
    conflict=conflict-spdiags(diag(conflict),0,nc,nc);
    degree=full(sum(conflict,2));[~,order]=sortrows([-degree,(1:nc)'],[1 2]);
    colors=zeros(nc,1);
    for j=1:nc
        cellID=order(j);neighbours=find(conflict(:,cellID));used=colors(neighbours);used=used(used>0);
        candidate=1;while any(used==candidate),candidate=candidate+1;end
        colors(cellID)=candidate;
    end
    pattern=struct('kind','axisymmetric-fv-colored-jacobian-pattern-v1', ...
        'cellPattern',cellPattern,'colors',colors,'geometryBinding',binding, ...
        'layout','component-block-Q-colon','supportScope','B*(B''*(I+H))+(I+H), full four-by-four blocks');
end
validateattributes(colors,{'numeric'},{'real','finite','integer','positive','numel',nc});
colorCount=max(colors);
assert(colorCount<=nc&&isequal(unique(colors),(1:colorCount)'), ...
    'AxisFV:JacobianColors','颜色编号须连续且覆盖全部单元。');
for color=1:colorCount
    assert(all(full(sum(cellPattern(:,colors==color),2))<=1), ...
        'AxisFV:JacobianColorConflict','同色列共享残差行，不能解码有色差分。');
end
pattern.colors=colors;
if nq==5,pattern.supportScope='B*(B''*(I+H))+(I+H), full five-by-five blocks';end

%% 三、每个输入分量单独按颜色求导；所有调用和回退均计数
report=struct('success',false,'status','building','residualEvaluations',0, ...
    'baselineResidualEvaluations',0,'positiveResidualEvaluations',0,'negativeResidualEvaluations',0, ...
    'failedResidualEvaluations',0,'colorCount',colorCount,'centralBatches',0,'oneSidedBatches',0, ...
    'retries',0,'splitBatches',0,'backwardBatches',0,'rejectedPerturbations',0, ...
    'unrepresentablePerturbations',0,'arithmeticRejections',0,'columnsDifferentiated',0, ...
    'layout','component-block-Q-colon','differenceScheme',options.differenceScheme, ...
    'differenceStep',options.differenceStep,'maxPerturbationRetries',options.maxPerturbationRetries, ...
    'patternReused',reused,'cellPatternEntries',nnz(cellPattern), ...
    'failureIdentifiers',{{}},'failureCounts',zeros(0,1), ...
    'lastFailureIdentifier','','lastFailureMessage','', ...
    'scale',scale,'dtDifferentiated',false,'globalScaleDifferentiated',false, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false, ...
    'newSolvedSurfaceAreaM2',0, ...
    'method','Explicit local-conservative colored finite differences of unchanged spatial rate; legacy Jacobian untouched');
J=[];
[baseOK,base]=evaluate(Q,0);
if ~baseOK,report.status='failed-baseline-residual';return;end
[rowCell,colCell]=find(cellPattern);nzc=numel(rowCell);
offset=[0;cumsum(full(sum(cellPattern,1))')];
values=zeros(nzc,nq,nq);computed=false(nc,nq);columnSchemes=zeros(nc,nq,'int8');acceptedPlus=NaN(nc,nq);acceptedMinus=NaN(nc,nq);
for component=1:nq
    for color=1:colorCount
        cells=find(colors==color);
        if ~differentiate(cells,component)
            report.status='failed-perturbation';report.columnSchemes=columnSchemes;
            report.columnsDifferentiated=nnz(computed);return;
        end
    end
end
assert(all(computed,'all'),'AxisFV:JacobianIncomplete','未计算的导数列不得以零值交付。');
count=nq^2*nzc;rows=zeros(count,1);columns=rows;entries=rows;cursor=0;
for input=1:nq
    for output=1:nq
        range=cursor+(1:nzc);rows(range)=rowCell+(output-1)*nc;
        columns(range)=colCell+(input-1)*nc;entries(range)=values(:,output,input);cursor=cursor+nzc;
    end
end
assert(all(isfinite(entries)),'AxisFV:JacobianFinite','已计算的Jacobian必须有限。');
J=sparse(rows,columns,entries,nq*nc,nq*nc);
report.success=true;report.status='computed-complete-colored-jacobian';
report.columnsDifferentiated=nnz(computed);report.columnSchemes=columnSchemes;
report.localScaling=scaleAudit;report.acceptedPositiveNormalizedStep=acceptedPlus;report.acceptedNegativeNormalizedStep=acceptedMinus;
report.structuralEntries=count;report.storedNonzeroEntries=nnz(J);
assert(report.residualEvaluations==report.baselineResidualEvaluations+ ...
    report.positiveResidualEvaluations+report.negativeResidualEvaluations, ...
    'AxisFV:JacobianCount','残差调用计数不闭合。');

    function [ok,value]=evaluate(state,direction)
        % 在实际调用之前计数，因此被状态/温域/HLLC/特征边界拒绝的调用也不会遗漏。
        report.residualEvaluations=report.residualEvaluations+1;
        if direction==0,report.baselineResidualEvaluations=report.baselineResidualEvaluations+1;
        elseif direction>0,report.positiveResidualEvaluations=report.positiveResidualEvaluations+1;
        else,report.negativeResidualEvaluations=report.negativeResidualEvaluations+1;end
        ok=false;value=[];
        try
            if nq==5,rate=axisymmetricRANSResidual(state,mesh,gas,boundary);
            else,rate=axisymmetricNavierStokesResidual(state,mesh,gas,boundary);end
            value=rate./scale;
            assert(all(isfinite(value),'all'),'AxisFV:JacobianNormalizedRate','归一化残差非有限。');
            ok=true;
        catch exception
            report.failedResidualEvaluations=report.failedResidualEvaluations+1;
            report.rejectedPerturbations=report.rejectedPerturbations+double(direction~=0);
            failure(exception.identifier,exception.message);
            if direction==0||~retryable(exception.identifier),rethrow(exception);end
        end
    end

    function failure(identifier,message)
        report.lastFailureIdentifier=identifier;report.lastFailureMessage=message;
        index=find(strcmp(report.failureIdentifiers,identifier),1);
        if isempty(index)
            report.failureIdentifiers{end+1,1}=identifier;report.failureCounts(end+1,1)=1;
        else,report.failureCounts(index)=report.failureCounts(index)+1;end
    end

    function [ok,value,actual]=directionTrial(cells,component,h,direction)
        state=Q;state(cells,component)=(U(cells,component)+direction*h)*scale(component);
        actual=direction*(state(cells,component)-Q(cells,component))/scale(component);
        ok=false;value=[];
        if ~all(isfinite(actual)&actual>0)||~all(isfinite(state(cells,component)))
            report.unrepresentablePerturbations=report.unrepresentablePerturbations+1;
            report.rejectedPerturbations=report.rejectedPerturbations+1;
            failure('AxisFV:JacobianStep','扰动在当前浮点量级不可表示或非有限；未调用残差。');return;
        end
        [ok,value]=evaluate(state,direction);
    end

    function ok=store(cells,component,positive,hp,negative,hm,scheme)
        % 中央扰动的实际正负步长可能有舍入差，用非均匀三点权重而非假定完全对称。
        pieces=cell(numel(cells),1);ok=false;
        for local=1:numel(cells)
            cellID=cells(local);indices=(offset(cellID)+1):offset(cellID+1);rr=rowCell(indices);
            if scheme==2
                weight=hm(local)/(hp(local)+hm(local));
                pieces{local}=weight*(positive(rr,:)-base(rr,:))/hp(local)+ ...
                    (1-weight)*(base(rr,:)-negative(rr,:))/hm(local);
            elseif scheme==1,pieces{local}=(positive(rr,:)-base(rr,:))/hp(local);
            else,pieces{local}=(base(rr,:)-negative(rr,:))/hm(local);end
            if ~all(isfinite(pieces{local}),'all')
                report.arithmeticRejections=report.arithmeticRejections+1;
                failure('AxisFV:JacobianDifference','有限残差差分得到非有限导数；未填零。');return;
            end
        end
        for local=1:numel(cells)
            cellID=cells(local);indices=(offset(cellID)+1):offset(cellID+1);
            values(indices,:,component)=pieces{local};computed(cellID,component)=true;
            columnSchemes(cellID,component)=int8(scheme);
            if scheme==2||scheme==1,acceptedPlus(cellID,component)=hp(local);end
            if scheme==2||scheme==-1,acceptedMinus(cellID,component)=hm(local);end
        end
        if scheme==2,report.centralBatches=report.centralBatches+1;
        else,report.oneSidedBatches=report.oneSidedBatches+1;
            if scheme<0,report.backwardBatches=report.backwardBatches+1;end
        end
        ok=true;
    end

    function ok=differentiate(cells,component)
        ok=false;h=options.differenceStep*localReference(cells,component);
        fallback=[];
        for attempt=0:options.maxPerturbationRetries
            if attempt>0,h=h/2;report.retries=report.retries+1;end
            [plusOK,plus,hp]=directionTrial(cells,component,h,1);
            if strcmp(options.differenceScheme,'forward')&&plusOK
                if store(cells,component,plus,hp,[],[],1),ok=true;return;end
            end
            % 中央方案总尝试正负两向；forward的正向失败时才尝试后向。
            if strcmp(options.differenceScheme,'central')||~plusOK
                [minusOK,minus,hm]=directionTrial(cells,component,h,-1);
            else,minusOK=false;minus=[];hm=[];end
            if strcmp(options.differenceScheme,'central')
                if plusOK&&minusOK
                    if store(cells,component,plus,hp,minus,hm,2),ok=true;return;end
                end
                if isempty(fallback)
                    if plusOK,fallback=struct('sign',1,'value',plus,'step',hp);
                    elseif minusOK,fallback=struct('sign',-1,'value',minus,'step',hm);end
                end
            elseif minusOK
                if store(cells,component,[],[],minus,hm,-1),ok=true;return;end
            end
        end
        if ~isempty(fallback)
            if fallback.sign>0
                ok=store(cells,component,fallback.value,fallback.step,[],[],1);
            else,ok=store(cells,component,[],[],fallback.value,fallback.step,-1);end
            if ok,return;end
        end
        % 同色中不同单元可能分别只接受相反方向；拆列后再试，不能整批补零。
        if numel(cells)>1
            report.splitBatches=report.splitBatches+1;
            for local=1:numel(cells)
                if ~differentiate(cells(local),component),return;end
            end
            ok=true;
        end
    end
end

function binding=geometryBinding(mesh)
%GEOMETRYBINDING 保存影响依赖支持及几何意义的完整字段；缓存不能跨网格误用。
binding=struct('kind','axisymmetric-fv-jacobian-geometry-binding-v1', ...
    'X',mesh.X,'R',mesh.R,'cells',mesh.cells,'planarArea',mesh.planarArea, ...
    'planarCentroid',mesh.planarCentroid,'volumePerRadian',mesh.volumePerRadian, ...
    'volumeCentroid',mesh.volumeCentroid,'faces',mesh.faces,'boundaryFaces',mesh.boundaryFaces, ...
    'gradientX',mesh.gradientX,'gradientR',mesh.gradientR,'incidence',mesh.incidence, ...
    'cellFaceMultiplicity',mesh.cellFaceMultiplicity,'periodicI',mesh.periodicI);
end

function yes=retryable(identifier)
%RETRYABLE 仅状态/波支/数值扰动可退让；输入、气体和几何合同错误不得吞掉。
yes=ismember(identifier,{'AxisFV:State','AxisFV:Temperature','AxisFV:HLLC','AxisFV:HLLE', ...
    'AxisSA:State','AxisSA:Arithmetic','AxisFV:SAModelState', ...
    'AxisFV:AxisTraceState','AxisSA:AxisTraceState', ...
    'AxisFV:TransportTemperature','AxisFV:TransportState','AxisFV:WallThermalState', ...
    'AxisFV:PressureSource','AxisFV:OutflowMach', ...
    'RocketAeroHeat:NSDraft:CharacteristicSonic','RocketAeroHeat:NSDraft:CharacteristicRegime', ...
    'RocketAeroHeat:NSDraft:CharacteristicAmbiguous','RocketAeroHeat:NSDraft:CharacteristicState', ...
    'RocketAeroHeat:NSDraft:CharacteristicSound'});
end
