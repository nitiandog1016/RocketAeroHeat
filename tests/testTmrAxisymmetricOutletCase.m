function report=testTmrAxisymmetricOutletCase(folder)
%TESTTMRAXISYMMETRICOUTLETCASE 冻结出口数值对照的来源/状态/实际空间账校核。
% 模块一：同网格同Q只改变明确的iMax四量及SA储库；给定量不当外流真值。
% 模块二：逐面核对冻结owner状态，拒绝非源点/非有限空间查询，不静默外推。
% 模块三：实际完整RANS残差及五量账，不以同Q变化称新稳态或实验通过。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v67_validation','outlet_contract');end
c=createTmrAxisymmetricOutletCase;b=createTmrAxisymmetricSensitivityCase('baseline');
assert(isequaln(c.mesh,b.mesh)&&isequaln(c.gas,b.gas)&&isequal(c.Q,b.Q)&&isequal(c.residualScale,b.residualScale));
for key={'iMin','jMin','jMax'},assert(isequaln(c.boundary.(key{1}),b.boundary.(key{1})));end
assert(strcmp(c.boundary.iMax.type,'riemannReservoir')&&isequaln(c.boundaryStudy.referenceBoundary,b.boundary.iMax));
p=c.boundaryStudy.profile;P=axisymmetricFVState(c.Q(:,1:4),c.gas);
assert(~p.dynamicZeroGradient&&~p.interpolated&&~p.extrapolated&&~p.experimentalValidation);
for k=1:numel(p.faceID)
    point=c.mesh.faces.centroid(p.faceID(k),:);owner=c.mesh.faces.owner(p.faceID(k));
    assert(isequal(c.boundary.iMax.primitive(point),P(owner,1:4))&& ...
        c.boundary.iMax.nuTilde(point)==c.Q(owner,5)/c.Q(owner,1));
end
refusals=0;
for point={p.centroidXR(1,:)+[1e-4 0],p.centroidXR(1,:)+[0 1e-4],[NaN 1],[1 2 3]}
    try,c.boundary.iMax.primitive(point{1});catch ex
        assert(strcmp(ex.identifier,'AxisFV:OutletProfilePoint'));refusals=refusals+1;continue;
    end
    error('AxisFV:OutletProfileTest','非法出口空间点未拒绝。');
end
[rate,info]=axisymmetricRANSResidual(c.Q,c.mesh,c.gas,c.boundary);
[oldRate,oldInfo]=axisymmetricRANSResidual(b.Q,b.mesh,b.gas,b.boundary);
assert(any(rate~=oldRate,'all')&&~info.physicalSolution&&~info.eligibleForSurfaceCoverage);
part=rate.*c.mesh.volumePerRadian;flux=info.faceFluxIntegralPerRadian(c.mesh.boundaryFaces,:);
source=info.cellSourceIntegral;defect=sum(part,1)+sum(flux,1)-sum(source,1);
scale=sum(abs(part),1)+sum(abs(flux),1)+sum(abs(source),1);
budget=2e-12*scale+10000*eps*sum(abs(info.faceFluxIntegralPerRadian),1);
assert(all(abs(defect)<=budget),'AxisFV:OutletProfileLedger','新出口数值算子的五量边界源账未闭合。');
report=struct('status','passed-explicit-frozen-outlet-reservoir-spatial-contract', ...
    'source',c.parentBinding,'outletFaces',numel(p.faceID),'refusals',refusals, ...
    'sameQOriginalResidual',max(abs(rate)./c.residualScale,[],'all'), ...
    'legacyOriginalResidual',max(abs(oldRate)./c.residualScale,[],'all'), ...
    'fiveQuantityDefect',defect,'fiveQuantityBudget',budget,'steadySolvePerformed',false, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'outlet_contract_evidence.mat'),'report','c','rate','info','oldRate','oldInfo');
fid=fopen(fullfile(folder,'outlet_contract_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Frozen outlet profile contract: %d real faces, %d refusals, five-quantity ledger; no steady/physical claim.\n',numel(p.faceID),refusals);
end
