function report=testRocketAxisymmetricCollarCase(folder)
%TESTROCKETAXISYMMETRICCOLLARCASE 核对新真实近壁网格的初态、完整五量方程和边界账。
% 模块一：同24层的两网格具有同一原墙、远场、气体、给定量和均匀Q，真实体积各自绑定。
% 模块二：实际调用完整NS+SA及墙热/应力/内部功，不把初态空间残差称为已推进解。
% 模块三：初态signed墙热和全五量账独立重算；策略/超域输入明确拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v70_validation','collar_case');end
item=createRocketAxisymmetricCollarCase;reference=createRocketAxisymmetricCollarCase(struct('gridPolicy','star-ray-reference'));
assert(isequaln(item,createRocketAxisymmetricCollarCase(item.options))&& ...
    isequaln(reference,createRocketAxisymmetricCollarCase(reference.options))&& ...
    isequal(item.Q,reference.Q)&&isequal(item.gas,reference.gas)&&isequal(item.boundary,reference.boundary)&& ...
    isequal(item.residualScale,reference.residualScale)&&isequaln(item.surfaceMesh,reference.surfaceMesh)&& ...
    isequal(item.domain.wallPolylineXR,reference.domain.wallPolylineXR)&& ...
    isequal(item.domain.outerBoundaryPolylineXR,reference.domain.outerBoundaryPolylineXR)&& ...
    ~isequal(item.mesh.X,reference.mesh.X)&&item.mesh.shape(2)==24&&reference.mesh.shape(2)==24);
assert(isequal(item.initialization.initialConservedIntegralPerRadian,sum(item.Q.*item.mesh.volumePerRadian,1))&& ...
    ~item.initialization.initialStateIsSolution&&~item.initialization.forebodyProfileImported&& ...
    ~item.initialization.externalCFDSolutionImported&&~item.inputContract.eligibleForSurfaceCoverage);
records=cell(2,1);rates=cell(2,1);infos=cell(2,1);cases={item,reference};
for k=1:2
    c=cases{k};[rate,info]=axisymmetricRANSResidual(c.Q,c.mesh,c.gas,c.boundary);
    assert(size(rate,2)==5&&all(isfinite(rate),'all')&&isequal(size(rate),size(c.Q))&& ...
        strcmp(info.viscousWorkVelocity.policy,'arithmetic')&&~isfield(info,'wallVelocityGradient')&& ...
        all(info.primitive(:,1)>0)&&all(info.sa.nuTilde>=0)&&~info.physicalSolution);
    flux=info.faceFluxIntegralPerRadian;src=info.cellSourceIntegral;volume=c.mesh.volumePerRadian;
    direct=(-c.mesh.incidence*flux+src)./volume;assert(isequal(direct,rate));
    defect=sum(rate.*volume,1)+info.netBoundaryFlux-info.sourceIntegral;
    formation=sum(abs(rate).*volume,1)+sum(abs(flux),1)+sum(abs(src),1);
    assert(all(abs(defect)<=2048*eps*formation),'RocketCollarCaseTest:Ledger','初态完整五量账不闭合。');
    w=info.wall;t=info.wallThermalGradient;[found,pos]=ismember(w.faceID,t.faceID);assert(all(found));
    q=-w.thermalConductivity.*t.normalGradientKm(pos);
    assert(all(w.heatFluxIntoWall<0)&&all(abs(q-w.heatFluxIntoWall)<=256*eps*(abs(q)+abs(w.heatFluxIntoWall)))&& ...
        all(w.temperatureInput==300)&&all(t.tangentialGradientInputKm==0)&& ...
        all(info.sa.faceMuT(w.faceID)==0)&&all(w.mechanicalPowerIntoFluidPerRadian==0));
    records{k}=struct('gridPolicy',c.inputContract.gridPolicy,'cellCount',c.mesh.cellCount, ...
        'initialScaledResidual',max(abs(rate)./c.residualScale,[],'all'),'fiveQuantityDefect',defect, ...
        'formationBudget',2048*eps*formation,'temperatureRangeK',[min(info.primitive(:,4)),max(info.primitive(:,4))], ...
        'signedWallHeatRangeWm2',[min(q),max(q)],'initialStateOnly',true);
    rates{k}=rate;infos{k}=info;
end
refusals=struct('label',{},'identifier',{});
reject('unknown-factory-option',@()createRocketAxisymmetricCollarCase(struct('guessSolution',true)));
reject('unknown-grid-policy',@()createRocketAxisymmetricCollarCase(struct('gridPolicy','smoothed')));
reject('reference-ignores-collar',@()createRocketAxisymmetricCollarCase(struct('gridPolicy','star-ray-reference','collar',struct('normalLayers',3))));
reject('nonzero-angle',@()createRocketAxisymmetricCollarCase(struct('flight',struct('alpha',20))));
reject('negative-SA-input',@()createRocketAxisymmetricCollarCase(struct('nuTildeRatio',-1)));
report=struct('status','passed-initial-complete-rocket-collar-and-same-layer-reference-spatial-equations', ...
    'matlabEngineProcessID',feature('getpid'),'cases',{records},'refusals',refusals, ...
    'physicalSolution',false,'flowAdvancePerformed',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'collar_case_evidence.mat'),'report','cases','rates','infos');
fid=fopen(fullfile(folder,'collar_case_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Collar/reference initial full five-quantity equations and signed wall heat passed; no accepted step yet.\n');
    function reject(label,action)
        try,action();catch e,assert(~isempty(e.identifier));refusals(end+1)=struct('label',label,'identifier',e.identifier);return;end
        error('RocketCollarCaseTest:Refusal','非法工况%s没有拒绝。',label);
    end
end
