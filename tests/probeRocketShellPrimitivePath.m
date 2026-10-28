function report=probeRocketShellPrimitivePath(folder)
%PROBEROCKETSHELLPRIMITIVEPATH 在同一真实Newton方向上核对新路径和完整NS残差下降。
% 模块一：原保守路径精确保持、零步身份及新路径链式切向导数校核。
% 模块二：真实三维第2/3方向的完整/半步调用原完整残差，越域候选明确单列。
% 模块三：试算结果与实际已接受半步分开，不推进物理时间或计生产面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'ShellTrialProbe:Overwrite','不得覆盖真实路径试算。');mkdir(folder);
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
s=load(fullfile(root,'outputs','v82_validation','exact_krylov01','steady_krylov_evidence.mat'),'result','gas','boundary');
scale=s.result.residualScale;gas=s.gas;boundary=s.boundary;rows=struct([]);candidates={};
for iteration=2:3
    t=s.result.trace{iteration};Q=t.initialQ;direction=t.linearUpdateU;
    old=rocketShellTrialState(Q,direction,scale,1,gas,'conservative');
    assert(isequaln(old,Q+direction.*scale)&&isequaln(rocketShellTrialState(Q,direction,scale,0,gas,'primitive-tangent'),Q));
    h=1e-6/max(abs(direction),[],'all');plus=rocketShellTrialState(Q,direction,scale,h,gas,'primitive-tangent');
    tangentError=max(abs((plus-Q)./scale/h-direction),[],'all')/max(abs(direction),[],'all');
    assert(tangentError<2e-5,'ShellTrialProbe:Tangent','新路径未保持同一一阶Newton方向。');
    for alpha=[1 .5]
    candidate=rocketShellTrialState(Q,direction,scale,alpha,gas,'primitive-tangent');timer=tic;identifier='';message='';
    try
        [rate,info]=cutCellNavierStokesResidual(candidate,d.domain,gas,boundary);
    catch ex
        identifier=ex.identifier;message=ex.message;
        if ~strcmp(identifier,'CutExactEuler:Temperature'),rethrow(ex);end
    end
    elapsed=toc(timer);merit=NaN;maximum=NaN;spatial=NaN;descent=false;temperatureRange=[NaN NaN];
    if isempty(identifier)
    z=rate./scale;merit=sqrt(mean(z.^2,'all'));maximum=max(abs(z),[],'all');
    f=d.domain.openFaces;w=d.domain.walls;outside=f.cellPlus==0;
    denominator=max(1,sum(abs(info.faceFlux(outside,:)).*f.area(outside),1)+sum(abs(info.wallFlux).*w.area,1));
    spatial=max(abs(info.balanceResidual)./denominator);assert(spatial<=1e-9);
    descent=merit<=(1-s.result.options.sufficientDecrease*alpha)*t.initialRmsMerit;
    temperatureRange=[min(info.primitive(:,5)),max(info.primitive(:,5))];
    candidates{end+1}=struct('sourceIteration',iteration,'alpha',alpha,'Q',candidate,'rate',rate,'info',info); %#ok<AGROW>
    else
    candidates{end+1}=struct('sourceIteration',iteration,'alpha',alpha,'Q',candidate,'errorIdentifier',identifier,'message',message); %#ok<AGROW>
    end
    row=struct('sourceIteration',iteration,'alpha',alpha,'tangentRelativeDifference',tangentError, ...
        'initialRmsMerit',t.initialRmsMerit,'originalAcceptedHalfStepRmsMerit',t.rmsMerit, ...
        'fullPrimitiveStepRmsMerit',merit,'fullPrimitiveStepMaximumResidual',maximum, ...
        'originalDescentRuleMet',descent,'spatialConservationRelative',spatial, ...
        'temperatureRangeK',temperatureRange,'errorIdentifier',identifier,'message',message,'elapsedSeconds',elapsed);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    fprintf('PRIMITIVE_PATH_FULL_NS iteration=%d alpha=%.3g rms=%.9g maximum=%.9g descent=%d error=%s\n',iteration,alpha,merit,maximum,descent,identifier);
    end
end
assert(rows(1).originalDescentRuleMet&&rows(2).originalDescentRuleMet&& ...
    strcmp(rows(3).errorIdentifier,'CutExactEuler:Temperature')&&rows(4).originalDescentRuleMet, ...
    'ShellTrialProbe:Scope','实际方向的下降或原温域拒绝未复现。');
report=struct('status','tested-real-primitive-path-with-full-NS-descent-and-explicit-interface-temperature-rejection','rows',rows, ...
    'originalEquationOrBoundaryChanged',false,'modelTemperatureRangeChanged',false,'physicalTimeAdvanced',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'primitive_path_evidence.mat'),'report','candidates','gas','boundary','-v7');
write(fullfile(folder,'primitive_path_summary.json'),jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
