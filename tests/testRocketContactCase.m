function report=testRocketContactCase(folder)
%TESTROCKETCONTACTCASE 原HLLC受控工况同Q/同网格/同物理边界完整比较。
% 核对黏性、墙、源和实际开放面保持；内部实际质量通量供同一SA第五量。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v78_validation','factory_test');end
assert(~exist(folder,'dir'),'RocketContactTest:Overwrite','不能覆盖实际测试。');rows=struct([]);
for kinds={'collar72','nested144'}
    kind=kinds{1};c=createRocketAxisymmetricContactCase(struct('kind',kind));assert(isequaln(c,createRocketAxisymmetricContactCase(c.options)));
    old=c.boundary;old.interiorConvectiveFlux='hlle-einfeldt';[r0,i0]=axisymmetricRANSResidual(c.Q,c.mesh,c.gas,old);
    [r1,i1]=axisymmetricRANSResidual(c.Q,c.mesh,c.gas,c.boundary);
    assert(isequal(i0.faceViscousFluxPerRadian,i1.faceViscousFluxPerRadian)&&isequal(i0.cellSourceIntegral,i1.cellSourceIntegral)&& ...
        isequaln(i0.wall,i1.wall)&&isequal(i0.faceConvectiveFluxPerRadian(c.mesh.boundaryFaces,:),i1.faceConvectiveFluxPerRadian(c.mesh.boundaryFaces,:))&& ...
        isequal(i1.sa.faceMassFluxPerRadian,i1.faceConvectiveFluxPerRadian(:,1))&& ...
        strcmp(i1.interiorConvectiveFlux.policy,'hllc')&&~i1.interiorConvectiveFlux.stateClipped&&~i1.interiorConvectiveFlux.adaptiveFluxFallback);
    assert(~isequal(r0,r1)&&max(abs(r0)./c.residualScale,[],'all')<=2e-6&& ...
        ~c.initialization.seed.initialGuessIsBoundaryCondition&&~isfield(c.initialization,'spatialRootHistory'));
    row=struct('kind',kind,'cellCount',c.mesh.cellCount,'sameConservativeInitialQUsed',true,'physicalTermsAndOpenBoundariesPreserved',true);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
report=struct('status','passed-same-real-Q-grid-physical-terms-original-HLLC-and-mass-consistent-SA', ...
    'cases',rows,'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'contact_factory_test_evidence.mat'),'report');
fid=fopen(fullfile(folder,'contact_factory_test_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
