function [selected,report]=selectRocketMinimumStateDirection(cache,direction,rate,target)
%SELECTROCKETMINIMUMSTATEDIRECTION 在有界15维子空间内缩小当前方向的真实状态幅度。
% 按第4/5层与其余层、五个原始变量拆分；每列经完整当前NS模型重新作用。
% 原方向恰在该子空间内，不继承旧方向资格；有限子空间不是全局最优声明。
% 固定原0.1之内的残差球，最小化相对rho/u/a/T幅度，禁止状态裁剪或波系冻结。
Q=cache.Q;n=size(Q,1);assert(mod(n,11)==0,'MinimumState:Layers','本研究入口要求真实11层箭体。');
metric=buildRocketPrimitiveUpdateMetric(Q,cache.scale,cache.gas);raw=applyRocketPrimitiveUpdateMetric(metric,direction);
layer=floor((0:n-1)'/(n/11))+1;group=ones(n,1);group(layer==4)=2;group(layer==5)=3;
A=zeros(numel(Q),15);normalization=zeros(15,1);timer=tic;
for g=1:3
    for component=1:5
        column=(g-1)*5+component;z=zeros(size(Q));mask=group==g;z(mask,component)=raw(mask,component);
        normalization(column)=norm(z,'fro');assert(normalization(column)>0,'MinimumState:Basis','零状态列不能归一化。');
        value=invertRocketPrimitiveUpdateMetric(metric,z/normalization(column));image=applyCutCellSelectedDirection(cache,value);A(:,column)=image(:);
        fprintf('MINIMUM_STATE_MODEL column=%d/15\n',column);
    end
end
assert(isequal(size(rate),size(Q))&&all(isfinite(rate),'all'),'MinimumState:Rate','须传入当前原NS率。');
b=-rate./cache.scale;
[coefficients,selection]=solveRocketMinimumStateCoefficients(A,b(:),target);
selected=direction;
if selection.feasible
    z=zeros(size(Q));
    for g=1:3
        for component=1:5
            column=(g-1)*5+component;mask=group==g;z(mask,component)=raw(mask,component)*coefficients(column)/normalization(column);
        end
    end
    selected=invertRocketPrimitiveUpdateMetric(metric,z);
end
baselineModel=A*normalization;selectedModel=baselineModel;returnedCoefficients=normalization;
% 不可行时返回原方向；诊断像也必须属于实际返回方向，不能误报零候选的像。
if selection.feasible,selectedModel=A*coefficients;returnedCoefficients=coefficients;end
report=selection;report.elapsedSeconds=toc(timer);report.applied=selection.feasible;report.modelApplications=15;
report.baselineModelDefect=norm(baselineModel-b(:))/norm(b(:));report.normalization=normalization;report.coefficients=coefficients;
report.returnedCoefficients=returnedCoefficients;report.returnedModelDefect=norm(selectedModel-b(:))/norm(b(:));
report.originalStateMetricNorm=norm(raw,'fro');report.originalMaximumStateComponents=max(abs(raw),[],1);
report.selectedStateMetricNorm=norm(applyRocketPrimitiveUpdateMetric(metric,selected),'fro');
report.selectedMaximumStateComponents=max(abs(applyRocketPrimitiveUpdateMetric(metric,selected)),[],1);
report.selectedModel=reshape(selectedModel,size(Q));report.baselineModel=reshape(baselineModel,size(Q));
report.rhs=reshape(b,size(Q));report.globalMinimumClaimed=false;report.originalNSAcceptanceClaimed=false;
end
