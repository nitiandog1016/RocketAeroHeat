function operator=buildRocketSteadyGraphPreconditioner(shifted)
%BUILDROCKETSTEADYGRAPHPRECONDITIONER 显式候选：整体移除五量时间对角后重分解同线块。
% 模块一：C=B-I/tau=-Jgraph，仅作为右预条件近似；原完整NS和稳态方程不变。
% 模块二：五个守恒量一起处理，不单独移除质量/能量正则或改变物理边界。
% 模块三：重分解失败直接拒绝，不补对角、不裁状态，不把B旧因子冒作C因子。
assert(isfield(shifted,'wlsNormalFactor')&&shifted.fixedLinearPreconditioner&& ...
    isfinite(shifted.pseudoTimeStepSeconds)&&shifted.pseudoTimeStepSeconds>0, ...
    'SteadyGraph:Operator','需要已验证且带明确时间对角的五带近似算子。');
assert(~isfield(shifted,'steadyGraphPreconditionerReport'),'SteadyGraph:AlreadyUnshifted','不能重复移除时间对角。');
operator=shifted;nf=operator.sourceFaces;nl=operator.layers;n=nf*nl;inverseTime=1/operator.pseudoTimeStepSeconds;
for a=1:5,operator.diagonal(a,a,:)=operator.diagonal(a,a,:)-inverseTime;end
inverse=zeros(5,5,n);upper1=inverse;upper2=inverse;lower1=inverse;identity=repmat(eye(5),1,1,nf);maximumDefect=0;timer=tic;
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);D=operator.diagonal(:,:,ids);L=operator.lower(:,:,ids);U=operator.upper(:,:,ids);
    if layer>2
        L=L-pagemtimes(operator.secondLower(:,:,ids),upper1(:,:,ids-2*nf));
        D=D-pagemtimes(operator.secondLower(:,:,ids),upper2(:,:,ids-2*nf));
    end
    if layer>1,D=D-pagemtimes(L,upper1(:,:,ids-nf));U=U-pagemtimes(L,upper2(:,:,ids-nf));end
    solved=pagemldivide(D,cat(2,identity,U,operator.secondUpper(:,:,ids)));
    assert(all(isfinite(solved),'all'),'SteadyGraph:Factor','未正则化五带块不可分解；不自动补对角。');
    inverse(:,:,ids)=solved(:,1:5,:);upper1(:,:,ids)=solved(:,6:10,:);upper2(:,:,ids)=solved(:,11:15,:);lower1(:,:,ids)=L;
    maximumDefect=max(maximumDefect,max(abs(pagemtimes(D,inverse(:,:,ids))-identity),[],'all'));
end
assert(maximumDefect<1e-8,'SteadyGraph:FactorDefect','未正则化同线块未达到既有分解门槛。');
operator.wlsNormalFactor=struct('inverse',inverse,'upper1',upper1,'upper2',upper2,'effectiveLower1',lower1);
operator.kind='unshifted-five-conservation-equation-WLS-normal-and-colored-side-steady-graph-preconditioner';
operator.steadyGraphPreconditionerReport=struct('linearApproximation','C=B-I/tau=-Jgraph', ...
    'allFiveConservativeTimeDiagonalsRemoved',true,'removedInverseTimeStepPerSecond',inverseTime, ...
    'maximumFactorDefect',maximumDefect,'factorSeconds',toc(timer),'oldShiftedFactorReused',false, ...
    'physicalEquationChanged',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
