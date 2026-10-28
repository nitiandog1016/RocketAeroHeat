function [U,report,rate,info]=polishAxisymmetricBEState(U,U0,dt,mesh,gas,boundary,scale,target,options)
%POLISHAXISYMMETRICBESTATE 在真实浮点邻域作有预算的BE末期修正，不放宽方程目标。
% 模块一：只处理距原目标很近的状态；每个候选都完整重算NS/RANS，不使用残差替身。
% 模块二：以超出原逐分量目标的平方和作搜索目标，允许最大分量切换；最终仍要求原无穷范数目标。
% 模块三：每次只改一个归一化量的1/2/4个浮点间隔；以状态版本和提交索引保存可重放轨迹。
% 模块四：预算/无下降均显式退出；success=false不得作为接受BE或物理解。
if nargin<9,options=struct;end
defaults=struct('maximumEvaluations',512,'maximumSweeps',3,'maximumSearchCells',16,'triggerRatio',16);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'AxisPolish:Options','未知末期修正选项。');
for k=fieldnames(defaults)',if ~isfield(options,k{1}),options.(k{1})=defaults.(k{1});end,end
for k={'maximumEvaluations','maximumSweeps','maximumSearchCells'},validateattributes(options.(k{1}),{'numeric'},{'real','finite','integer','positive','scalar'});end
validateattributes(options.triggerRatio,{'numeric'},{'real','finite','scalar','>=',1});
validateattributes(dt,{'numeric'},{'real','finite','positive','scalar'});validateattributes(target,{'numeric'},{'real','finite','positive','scalar'});
assert(isequal(size(U),size(U0))&&size(U,1)==mesh.cellCount&&ismember(size(U,2),[4 5]),'AxisPolish:State','状态维度不匹配。');
scale=reshape(scale,1,[]);initialU=U;[rate,info]=evaluate(U);evaluations=1;G=U-U0-dt*(rate./scale);
[merit,norm]=measure(G);initialNorm=norm;initialMerit=merit;version=0;skipped=0;reason='maximum-sweeps';
log=struct('evaluation',{},'stateVersion',{},'sweep',{},'cell',{},'quantity',{},'stepFactor',{}, ...
    'beforeScaledValue',{},'candidateScaledValue',{},'norm',{},'merit',{},'errorIdentifier',{});
commits=struct('version',{},'evaluation',{},'cell',{},'quantity',{},'scaledValue',{},'norm',{},'merit',{});
success=norm<=target;
if norm>options.triggerRatio*target,reason='outside-near-target-trigger';
elseif success,reason='already-within-original-target';
else
    for sweep=1:options.maximumSweeps
        beforeMerit=merit;cellScore=max(abs(G),[],2);[~,order]=sort(cellScore,'descend');
        cells=order(1:min(options.maximumSearchCells,numel(order)));
        for ci=1:numel(cells)
            cellID=cells(ci);
            for quantity=1:size(U,2)
                bestMerit=merit;bestNorm=norm;bestValue=U(cellID,quantity);bestEvaluation=0;bestRate=[];bestInfo=[];bestG=[];
                for factor=[-1 1 -2 2 -4 4]
                    if evaluations>=options.maximumEvaluations,reason='evaluation-budget';break;end
                    candidate=U;candidate(cellID,quantity)=U(cellID,quantity)+factor*eps(U(cellID,quantity));
                    if candidate(cellID,quantity)*scale(quantity)==U(cellID,quantity)*scale(quantity),skipped=skipped+1;continue;end
                    identifier='';newNorm=NaN;newMerit=NaN;evaluations=evaluations+1;
                    try
                        [newRate,newInfo]=evaluate(candidate);newG=candidate-U0-dt*(newRate./scale);[newMerit,newNorm]=measure(newG);
                    catch exception
                        identifier=exception.identifier;
                        if ~ismember(identifier,{'AxisFV:State','AxisFV:Temperature','AxisSA:State','AxisFV:HLLC','AxisFV:HLLE', ...
                                'AxisFV:PressureSource','AxisFV:AxisTraceState','AxisSA:AxisTraceState','AxisFV:SAModelState', ...
                                'AxisFV:TransportTemperature','AxisFV:TransportState','AxisFV:WallThermalState','AxisSA:Arithmetic'}),rethrow(exception);end
                    end
                    log(end+1)=struct('evaluation',evaluations,'stateVersion',version,'sweep',sweep,'cell',cellID, ...
                        'quantity',quantity,'stepFactor',factor,'beforeScaledValue',U(cellID,quantity), ...
                        'candidateScaledValue',candidate(cellID,quantity),'norm',newNorm,'merit',newMerit,'errorIdentifier',identifier); %#ok<AGROW>
                    if isempty(identifier)&&newMerit<bestMerit
                        bestMerit=newMerit;bestNorm=newNorm;bestValue=candidate(cellID,quantity);bestEvaluation=evaluations;
                        bestRate=newRate;bestInfo=newInfo;bestG=newG;
                    end
                    if newNorm<=target,break;end
                end
                if bestEvaluation>0
                    U(cellID,quantity)=bestValue;rate=bestRate;info=bestInfo;G=bestG;merit=bestMerit;norm=bestNorm;version=version+1;
                    commits(end+1)=struct('version',version,'evaluation',bestEvaluation,'cell',cellID,'quantity',quantity, ...
                        'scaledValue',bestValue,'norm',norm,'merit',merit); %#ok<AGROW>
                end
                success=norm<=target;
                if success,reason='original-BE-target-met';break;end
                if evaluations>=options.maximumEvaluations,break;end
            end
            if success||evaluations>=options.maximumEvaluations,break;end
        end
        if success||evaluations>=options.maximumEvaluations,break;end
        if merit>=beforeMerit,reason='no-merit-descent';break;end
    end
end
report=struct('kind','bounded-representable-neighborhood-BE-polish-v1','success',success,'reason',reason, ...
    'options',options,'residualEvaluations',evaluations,'skippedUnchangedPhysicalProposals',skipped, ...
    'initialNorm',initialNorm,'finalNorm',norm,'originalTarget',target,'initialMerit',initialMerit,'finalMerit',merit, ...
    'commits',commits,'log',log,'stateChanged',~isequal(U,initialU), ...
    'infinityNormAcceptanceUnchanged',true,'conservationAcceptancePerformedByCaller',true, ...
    'meritDefinesSearchNotAcceptance',true,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    function [r,i]=evaluate(scaled)
        if size(scaled,2)==5,[r,i]=axisymmetricRANSResidual(scaled.*scale,mesh,gas,boundary);else,[r,i]=axisymmetricNavierStokesResidual(scaled.*scale,mesh,gas,boundary);end
    end
    function [m,n]=measure(g)
        n=max(abs(g),[],'all');excess=max(abs(g)/target-1,0);m=sum(excess.^2,'all');
    end
end
