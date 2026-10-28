function report=testAxisymmetricOpenBoundary(folder)
%TESTAXISYMMETRICOPENBOUNDARY 正反超声速/特征储库在轴对称四量上的合同。
% 本校核验证均匀迹及边界/源账，不代表非均匀真实边界层入口已闭合。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v62_validation','open_boundary');end
if ~exist(folder,'dir'),mkdir(folder);end
[X,R]=ndgrid([.1,.05,0],[.03,.04,.06]);mesh=buildAxisymmetricFVGeometry(X,R);
gas=struct('gamma',1.4,'R',287.05,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]);
wall=@(u)struct('type','wall','temperature',300,'velocity',[u 0]);
records=struct([]);infos=cell(1,4);
for k=1:4
    direction=1-2*mod(k,2);velocity=1100*direction;primitive=[1 velocity 0 300];
    Q=axisymmetricFVConservative(repmat(primitive,mesh.cellCount,1),gas);
    if k<=2
        inflow=struct('type','supersonicInflow','primitive',primitive);outflow=struct('type','supersonicOutflow');
        if velocity<0,bc=struct('iMin',inflow,'iMax',outflow);else,bc=struct('iMin',outflow,'iMax',inflow);end
        bc.jMin=wall(velocity);bc.jMax=wall(velocity);
    else
        far=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',[1 velocity 0 0 300]);
        bc=struct('iMin',far,'iMax',far,'jMin',far,'jMax',far);
    end
    [rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,bc);P=axisymmetricFVState(Q,gas);
    scale=[1,max(P(:,6)),max(P(:,6)),max(Q(:,4))];
    perStep=max(abs(rate)./scale./(info.convectiveFrequency+info.diffusiveFrequency),[],'all');
    assert(perStep<2e-12,'均匀全特征边界出现非舍入尺度改变。');
    assert(max(abs(info.conservationDefect)./max(1,abs(info.sourceIntegral)))<1e-12);
    assert(~info.physicalSolution&&~info.eligibleForSurfaceCoverage);
    item=struct('velocityX',velocity,'characteristic',k>2,'maximumLocalStepChange',perStep, ...
        'boundaryFlux',info.netBoundaryFlux,'axisSource',info.sourceIntegral); %#ok<AGROW>
    if isempty(records),records=item;else,records(k)=item;end
    infos{k}=info;
end
rejections=0;goodBC=bc;badBC=bc;badBC.iMax=struct('type','supersonicOutflow');
reject(@()axisymmetricNavierStokesResidual(Q,mesh,gas,badBC));
bad=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',[1 -1100 0 1 300]);
badBC=goodBC;badBC.iMin=bad;reject(@()axisymmetricNavierStokesResidual(Q,mesh,gas,badBC));
badBC=goodBC;badBC.iMin=struct('type','characteristic','mode','riemannFarfield');reject(@()axisymmetricNavierStokesResidual(Q,mesh,gas,badBC));
report=struct('status','passed-uniform-open-traces-only','records',records,'rejections',rejections, ...
    'actualBoundaryLayerInletValidated',false,'experimentalValidation',false,'newSolvedAreaM2',0);
save(fullfile(folder,'open_boundary_evidence.mat'),'mesh','gas','records','infos','report');
fid=fopen(fullfile(folder,'open_boundary_summary.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Axisymmetric open boundary passed: both axial directions,4-variable characteristic traces, %d refusals; no full inlet validation.\n',rejections);
    function reject(action)
        try,action();catch,rejections=rejections+1;return;end
        error('AxisFVTest:BoundaryRefusal','非法开放边界未被拒绝。');
    end
end
