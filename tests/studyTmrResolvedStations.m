function report=studyTmrResolvedStations(folder)
%STUDYTMRRESOLVEDSTATIONS 在公开风洞原节点上求解一个区分独立测站的网格。
% 模块一：绑定既有稳态并核对其壁量、完整面通量及单元源项。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
assert(isfolder(folder) && ~isfile(fullfile(folder,'initial.mat')),'试验目录已有初场或不存在。');
parentFile=fullfile(root,'outputs','v66_validation','tmr_wall_thermal','tmr_wall_thermal_evidence.mat');
assert(strcmp(fileSHA(parentFile),'d1b8d4479b23358611c653fd131e96759cf9556ed6555fbaea8c8a517694016f'));
d=load(parentFile,'mesh','implicit','gas','boundary');
gas=d.gas;boundary=d.boundary;parent=d.implicit;opts=parent.options;
assert(parent.steadyConverged && all(parent.history.residual(end-2:end)<=2e-6));
assert(opts.steadyResidualTolerance==2e-6 && opts.steadyRequiredSteps==3 && ...
    opts.conservationTolerance==5e-10 && boundary.jMin.temperature==311);
[parentRate,parentInfo]=axisymmetricRANSResidual(parent.Q,d.mesh,gas,boundary);
assert(isequal(parentInfo.faceFluxIntegralPerRadian,parent.finalInfo.faceFluxIntegralPerRadian) && ...
    isequal(parentInfo.cellSourceIntegral,parent.finalInfo.cellSourceIntegral) && ...
    isequal(parentInfo.wall.pressure,parent.finalInfo.wall.pressure) && ...
    isequal(parentInfo.wall.heatFluxIntoWall,parent.finalInfo.wall.heatFluxIntoWall), ...
    'TmrStations:Replay','父状态的原面通量、单元源项或壁量存在差异。');
assert(max(abs(parentRate)./parent.residualScale,[],'all')==parent.residualNorm,'父原残差不符。');

% 模块二：只用测站位置和原节点进行嵌套二分，各原法向分层保持不变。
grid=loadTmrMach7Grid(2);
baseI=1:4:grid.ni;jNodes=1:4:grid.nj;
assert(isequal(d.mesh.X,grid.X(baseI,jNodes)) && isequal(d.mesh.R,grid.R(baseI,jNodes)));
contractFile=fullfile(root,'tests','data','tmr_mach7','contract.json');
assert(strcmp(fileSHA(contractFile),'121b83d627f2ac573ea38dbcb155e177e321d2319c4f379474f4698d92935b19'));
c=jsondecode(fileread(contractFile));
stations={unique([c.measurements.pressure_rows.project_X_m]),unique([c.measurements.heating_rows.project_X_m])};
iNodes=baseI;
for pass=1:2
    split=false(numel(iNodes)-1,1);
    for k=1:numel(stations)
        idx=wallIntervals(grid.X(iNodes,1),stations{k});
        split=split | accumarray(idx,1,size(split))>1;
    end
    left=find(split);middle=(iNodes(left)+iNodes(left+1))/2;
    assert(all(middle==fix(middle)),'新增列必须存在于官方网格。');
    iNodes=sort([iNodes middle]);
end
assert(isequal(setdiff(iNodes,baseI),[67 71 75 87 111 115 119 122 123 127 131 143]));
mesh=buildAxisymmetricFVGeometry(grid.X(iNodes,jNodes),grid.R(iNodes,jNodes));
for k=1:numel(stations)
    idx=wallIntervals(mesh.X(:,1),stations{k});
    assert(numel(unique(idx))==numel(stations{k}),'仍有同类独立测站落在同一壁段。');
end
assert(mesh.cellCount==1300);
[i,j]=ind2sub(mesh.shape,(1:mesh.cellCount)');
coarseColumn=floor((iNodes(1:end-1)-1)/4)+1;
parentIDs=sub2ind(d.mesh.shape,reshape(coarseColumn(i),[],1),j);
Q=parent.Q(parentIDs,:);initialQ=Q;V=mesh.volumePerRadian;scale=parent.residualScale;
[rate,info]=axisymmetricRANSResidual(Q,mesh,gas,boundary);
initial=sum(Q.*V,1);parentIntegral=sum(parent.Q.*d.mesh.volumePerRadian,1);
save(fullfile(folder,'initial.mat'),'mesh','gas','boundary','initialQ','parentIDs','iNodes','jNodes', ...
    'initial','parentIntegral','scale','rate','info','opts','-v7');

% 模块三：逐个调用原生后向 Euler 步。每个已接受的完整状态先落盘，再判断稳态。
steps=0;rejections=0;evaluations=0;started=tic;cfl=opts.initialPseudoCFL;opts.maxSteps=1;
residual=max(abs(rate)./scale,[],'all');residuals=residual;net=zeros(1,5);src=zeros(1,5);
balance=0;steady=false;stop='step-budget';
while steps<75
    explicitStep=1/max(info.convectiveFrequency+info.diffusiveFrequency);
    opts.initialPseudoCFL=cfl;
    segment=advanceAxisymmetricSteady(Q,mesh,gas,boundary,opts);
    rejections=rejections+segment.rejectedAttempts;evaluations=evaluations+segment.residualEvaluations;
    if segment.acceptedSteps~=1
        save(fullfile(folder,'rejected.mat'),'segment','steps','-v7');
        stop=segment.stopReason;break
    end
    Q=segment.Q;steps=steps+1;
    [rate,info]=axisymmetricRANSResidual(Q,mesh,gas,boundary);
    residual=max(abs(rate)./scale,[],'all');
    assert(residual==segment.residualNorm && isequaln(info,segment.finalInfo),'接受状态未能逐值回读。');
    residuals(end+1)=residual; %#ok<AGROW>
    dt=segment.history.dt(end);
    net=net+dt*info.netBoundaryFlux;src=src+dt*info.sourceIntegral;
    balance=max(abs(sum(Q.*V,1)-initial+net-src)./(sum(V)*scale));
    cfl=min(opts.maximumPseudoCFL,opts.pseudoCFLGrowth*dt/explicitStep);
    history=segment.history;rejectLog=segment.rejectLog; %#ok<NASGU>
    file=fullfile(folder,sprintf('step_%03d.mat',steps));
    temp=fullfile(folder,'writing_step.mat');
    save(temp,'Q','residual','steps','dt','net','src','balance','cfl','history','rejectLog','-v7');
    assert(~isfile(file));movefile(temp,file);
    fprintf('station grid step %d residual %.9g balance %.4g elapsed %.2f s\n',steps,residual,balance,toc(started));
    assert(balance<=opts.conservationTolerance,'TmrStations:Balance','累计五量收支超出原容差。');
    steady=steps>=3 && all(residuals(end-2:end)<=opts.steadyResidualTolerance);
    if steady,stop='three-accepted-steady-endpoints';break;end
end
% 模块四：区分实际末场、停止原因与实验资格。物理时间及新增原箭体面积均为零。
implicit=struct('Q',Q,'finalInfo',info,'residualNorm',residual,'residualScale',scale, ...
    'history',struct('residual',residuals),'globalBalanceRelative',balance,'steadyConverged',steady, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false);
report=struct('cells',mesh.cellCount,'rings',mesh.shape(1),'acceptedSteps',steps,'rejectedAttempts',rejections, ...
    'solverResidualEvaluations',evaluations,'solverLoopSeconds',toc(started),'stopReason',stop, ...
    'initialResidual',residuals(1),'endResidual',residual,'globalBalanceRelative',balance,'steady',steady, ...
    'physicalTimeAdvanced',0,'givenTwK',311,'copiedParentStateIsInitialGuess',true,'originalRocketAreaAddedM2',0);
save(fullfile(folder,'terminal.mat'),'mesh','implicit','gas','boundary','report','-v7');
fid=fopen(fullfile(folder,'terminal.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
end

function idx=wallIntervals(x,stations)
% 原网格 X 向下游递减；末端点归入最后一个壁段，区间外保留拒绝。
idx=discretize(-stations(:),-x(:));assert(all(isfinite(idx)),'实验测站超出固壁范围。');
end

function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');md.update(fread(fid,inf,'*uint8'));
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
