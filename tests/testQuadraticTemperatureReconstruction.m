function report=testQuadraticTemperatureReconstruction(folder)
%TESTQUADRATICTEMPERATURERECONSTRUCTION 真实FV均值二次再现与输入/秩合同测试。
% 模块用途：只验证新标量空间重构，不代表完整NS温度闭合、时间稳定或物理覆盖。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,folder=fullfile(root,'outputs','v60_validation','library_reproduction');end
if ~isfolder(folder),mkdir(folder);end
maxNumCompThreads(1);addpath(fullfile(root,'internalFlow'));setupCutCellFlow;
val=validationFolder(root,'v59');prior=validationFolder(root,'v58');names={'actual','sphere_grid8','sphere_grid12','sphere_grid16'};
g=[8 -5 3];H=[30 4 -2;4 -20 3;-2 3 10];k=.5*1.4/(1.4-1)/.72;
point=struct('type','temperaturePoint','value',@(x,id)300+x*g'+.5*sum((x*H).*x,2));bc=allBoundaries(point);
evidence=cell(4,1);checks=cell(4,1);
for j=1:4
    if j==1
        a=load(fullfile(prior,'cut_cell_transport_draft','actual_aggregation_results.mat'),'group');
        d=load(fullfile(prior,'cut_cell_draft','cut_cell_geometry_results.mat'),'actual');source=d.actual;
    else,a=load(fullfile(val,'geometry',[names{j},'.mat']),'source','group');source=a.source;end
    b=load(fullfile(val,'curved_accuracy_oracle',[names{j},'_fv_volume_moments.mat']),'moments');c=a.group;
    m=buildCutCellVolumeMoments(source,c);
    field=struct('values',b.moments.quadraticTemperatureVolumeAverage(:),'valueKind','volume-average-temperature','units','K');
    started=tic;stencil=buildQuadraticTemperatureStencil(c,m,bc);out=quadraticTemperatureFourierFlux(field,c,k,bc,stencil);
    exactGradient=g+c.cells.centroid*H;gradientError=max(abs(out.reconstruction.gradient-exactGradient),[],'all');
    hessianError=max(abs(out.reconstruction.hessian-repmat(reshape(H,1,3,3),numel(c.cells.volume),1,1)),[],'all');
    wallExact=-k*sum((g+c.walls.centroid*H).*c.walls.normal,2);wallError=max(abs(out.wallHeatFluxIntoWall-wallExact));
    laplaceError=max(abs(out.thermalDiffusionRate-k*trace(H)));
    area=c.walls.area;rms=sqrt(sum(area.*(out.wallHeatFluxIntoWall-wallExact).^2)/sum(area));
    expectedPower=sum(area.*wallExact);powerError=sum(area.*out.wallHeatFluxIntoWall)-expectedPower;
    assert(gradientError<1e-7&&hessianError<1e-6&&wallError<1e-6,'CutQuadraticTest:Reproduction','二次FV场没有再现。');
    assert(laplaceError<1e-5,'CutQuadraticTest:Divergence','真实面积/体积上的二次Fourier散度不一致。');
    assert(abs(out.thermalBalance)<1e-8,'CutQuadraticTest:Conservation','共享热通量全域账不一致。');
    assert(all(stencil.rank==9)&&all(stencil.scaledRelativeSingularValue>1e-10));
    assert(sum(cellfun(@(s)numel(s.boundarySampleIDs),stencil.local))==numel(stencil.samples.owner),'CutQuadraticTest:Dust','不能遗漏小正墙片。');
    uniform=struct('type','temperaturePoint','value',300);ubc=allBoundaries(uniform);
    uniformField=field;uniformField.values(:)=300;u=quadraticTemperatureFourierFlux(uniformField,c,k,ubc,stencil);
    assert(all(u.wallHeatFluxIntoWall==0)&&all(u.thermalDiffusionRate==0));
    if j==2
        % 零法向梯度给定的是热流边界，不会自动求得墙温。
        nbc=ubc;nbc.solidWall=struct('type','zeroNormalGradient');
        ns=buildQuadraticTemperatureStencil(c,m,nbc);nu=quadraticTemperatureFourierFlux(uniformField,c,k,nbc,ns);
        assert(all(nu.wallHeatFluxIntoWall==0)&&all(isnan(nu.wallTemperatureInput)));
        assert(all(strcmp(nu.wallHeatValueKind,'prescribed-zero-normal-gradient')));
        assert(all(strcmp(nu.wallTemperatureValueKind,'unknown-wall-temperature')));
    end
    check=struct('caseName',names{j},'fluidGroups',numel(c.cells.volume),'gradientMaximumError',gradientError, ...
        'hessianMaximumError',hessianError,'wallHeatMaximumError',wallError,'wallHeatAreaWeightedRMS',rms, ...
        'thermalDivergenceMaximumError',laplaceError,'wallHeatPowerError',powerError, ...
        'thermalBalance',out.thermalBalance,'maximumNeighborRing',max(stencil.ringCount), ...
        'minimumScaledRelativeSingularValue',min(stencil.scaledRelativeSingularValue),'elapsedSeconds',toc(started));
    checks{j}=check;evidence{j}=struct('caseName',names{j},'meanField',field,'stencil',stencil,'output',out,'check',check);
    fprintf('QUADRATIC_FV_REPRODUCTION %s gradient=%.9g Hessian=%.9g wallq=%.9g laplace=%.9g rings=%d\n', ...
        names{j},gradientError,hessianError,wallError,laplaceError,check.maximumNeighborRing);
end
% 输入身份、模板失效与负Cov均明确拒绝，不以成功拟合掩盖错误数据。
badField=field;badField.valueKind='EOS-primitive-from-cell-average-Q';
reject(@()evaluateQuadraticTemperature(badField,c,bc,stencil),'CutQuadratic:MeanIdentity');
stale=c;stale.walls.area(1)=stale.walls.area(1)*1.001;
reject(@()evaluateQuadraticTemperature(field,stale,bc,stencil),'CutQuadratic:GeometryIdentity');
reject(@()buildQuadraticTemperatureStencil(stale,m,bc),'CutQuadratic:MomentGeometryIdentity');
missingIdentity=rmfield(m,'geometryIdentity');
reject(@()buildQuadraticTemperatureStencil(c,missingIdentity,bc),'CutQuadratic:MomentGeometryIdentity');
badGroup=c;badGroup.walls.centroid(1,1)=badGroup.walls.centroid(1,1)+.001;
reject(@()buildCutCellVolumeMoments(source,badGroup),'CutMom:GroupMap');
bad=m;bad.covariance(1,1,1)=-100;
reject(@()buildQuadraticTemperatureStencil(c,bad,bc),'CutQuadratic:Moments');
% 复数Hermitian正定矩也不是物理几何Cov，不能只用共轭对称/eig判据放过。
complexMoments=m;complexMoments.covariance(1,1,2)=complexMoments.covariance(1,1,2)+1i*1e-12;
complexMoments.covariance(1,2,1)=complexMoments.covariance(1,2,1)-1i*1e-12;
reject(@()buildQuadraticTemperatureStencil(c,complexMoments,bc),'CutQuadratic:Moments');
[isolated,im]=oneCell;
reject(@()buildQuadraticTemperatureStencil(isolated,im,bc),'CutQuadratic:InsufficientStencil');
reject(@()quadraticTemperatureFourierFlux(field,c,1e308,bc,stencil),'CutQuadratic:NonfiniteFlux');
complexField=field;complexField.values=complexField.values+1i;
reject(@()evaluateQuadraticTemperature(complexField,c,bc,stencil),'CutQuadratic:MeanValues');
report=struct('status','tested-experimental-library-quadratic-FV-reproduction','checks',{checks},'refusalGroups',10, ...
    'uniformPhysicalZeroRetained',true,'trueFVTemperatureInputRequired',true, ...
    'zeroNormalGradientDoesNotSolveWallTemperature',true,'finiteFluxRequired',true, ...
    'momentsBoundToCompleteCVGeometry',true, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'fullNavierStokesStabilityProven',false,'newSolvedSurfaceAreaM2',0);
save(fullfile(folder,'quadratic_temperature_evidence.mat'),'report','evidence','-v7');
fid=fopen(fullfile(folder,'quadratic_temperature_summary.json'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('QUADRATIC_FV_RECONSTRUCTION_PASS four fixed geometries, constant zero, explicit rank/input/cache refusals. Physical area0.\n');
end
function p=validationFolder(root,version)
p=fullfile(root,'outputs',[version,'_validation']);if ~isfolder(p),p=fullfile(root,'validation',version);end
end
function b=allBoundaries(p),b=struct('xmin',p,'xmax',p,'ymin',p,'ymax',p,'zmin',p,'zmax',p,'solidWall',p);end
function reject(action,id)
try,action();catch err,assert(strcmp(err.identifier,id),err.message);return;end
error('CutQuadraticTest:MissingRefusal','预期未拒绝：%s',id);
end
function [c,m]=oneCell
% 一个全格只有六个面中心样本，交叉曲率不可识别，必须拒绝而非伪造9维解。
c.cells=struct('volume',1,'centroid',[0 0 0]);
c.openFaces=struct('area',ones(6,1),'centroid',[-.5 0 0;.5 0 0;0 -.5 0;0 .5 0;0 0 -.5;0 0 .5], ...
    'normal',[1 0 0;1 0 0;0 1 0;0 1 0;0 0 1;0 0 1], ...
    'cellMinus',[0;1;0;1;0;1],'cellPlus',[1;0;1;0;1;0],'axis',[1;1;2;2;3;3]);
c.walls=struct('area',zeros(0,1),'ownerCell',zeros(0,1),'centroid',zeros(0,3),'normal',zeros(0,3),'originalFaceID',zeros(0,1));
m=struct('volume',1,'centroid',[0 0 0],'covariance',reshape(eye(3)/12,1,3,3));
m.geometryIdentity=struct('volume',1,'centroid',c.cells.centroid,'openFaces',c.openFaces,'walls',c.walls);
end
