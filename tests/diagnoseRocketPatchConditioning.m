function report=diagnoseRocketPatchConditioning(folder,fixture)
%DIAGNOSEROCKETPATCHCONDITIONING 固定真实Q24，以物理变量度量分解局部病态方向。
% 同一原NS Jacobian，比较精确Newton与达到原0.1的最小状态幅度方向。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);s=load(fixture,'patch');p=s.patch;Q=p.Q1;core=p.core;scale=p.scale;
[rate,~]=rocketFiniteStepPatchResidual(Q,p);rhs=rate./scale;n=numel(rhs);J=zeros(n);W=zeros(n);
metric=buildRocketPrimitiveUpdateMetric(Q(core,:),scale,p.gas);
for column=1:n
    [a,b]=ind2sub(size(rhs),column);h=1e-6*max(1,abs(Q(core(a),b)/scale(b)));plus=Q;minus=Q;
    plus(core(a),b)=plus(core(a),b)+h*scale(b);minus(core(a),b)=minus(core(a),b)-h*scale(b);
    rp=rocketFiniteStepPatchResidual(plus,p);rm=rocketFiniteStepPatchResidual(minus,p);col=(rp-rm)./scale/(2*h);J(:,column)=col(:);
    e=zeros(size(rhs));e(column)=1;z=applyRocketPrimitiveUpdateMetric(metric,e);W(:,column)=z(:);
end
A=J/W;[U,S,V]=svd(A,'econ');sv=diag(S);projection=U'*rhs(:);exact=J\(-rhs(:));
% 固定0.09留给原NS差分核对余量；原0.1门槛不变，非参数扫描。
target=.09;lo=0;hi=max(sv)^2;
for it=1:100
    lam=(lo+hi)/2;res=norm(lam./(sv.^2+lam).*projection)/norm(rhs(:));
    if res>target,hi=lam;else,lo=lam;end
end
z=-V*((sv./(sv.^2+lo)).*projection);small=W\z;
report=struct('singularValues',sv,'conditionNumber',sv(1)/sv(end), ...
    'exactRawStateMetricNorm',norm(W*exact),'inexactRawStateMetricNorm',norm(W*small), ...
    'exactLargestRawStateComponent',max(abs(W*exact)),'inexactLargestRawStateComponent',max(abs(W*small)), ...
    'inexactModelDefect',norm(J*small+rhs(:))/norm(rhs(:)),'physicalSolution',false);
directions=[exact small];trials=struct([]);
for kind=1:2
    d=zeros(size(Q));d(core,:)=reshape(directions(:,kind),size(rhs));h=1e-6/max(abs(d),[],'all');
    rp=rocketFiniteStepPatchResidual(Q+h*d.*scale,p);jv=(rp-rate)./scale/h;
    ns=norm(jv+rhs,'fro')/norm(rhs,'fro');
    for search=0:20
        alpha=2^-search;row=struct('kind',kind,'alpha',alpha,'originalNSDefect',ns,'rms',NaN,'maximum',NaN,'stateValid',false,'remainder',NaN);
        try
            rn=rocketFiniteStepPatchResidual(Q+alpha*d.*scale,p);row.stateValid=true;rr=rn./scale;row.rms=sqrt(mean(rr.^2,'all'));row.maximum=max(abs(rr),[],'all');
            e=rocketFiniteStepAgreement(rate,rn,jv,scale,alpha);row.remainder=e.relativeNonlinearRemainder;
        catch exception
            if ~startsWith(exception.identifier,'CutNS:')&&~startsWith(exception.identifier,'CutExactEuler:'),rethrow(exception);end
        end
        if isempty(trials),trials=row;else,trials(end+1)=row;end %#ok<AGROW>
    end
end
save(fullfile(folder,'conditioning.mat'),'report','J','W','A','rhs','exact','small','trials','-v7');
fid=fopen(fullfile(folder,'conditioning.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(rmfield(report,'singularValues'));
end
