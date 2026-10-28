function report=studyCubicTemperatureReconstruction(folder)
%STUDYCUBICTEMPERATURERECONSTRUCTION 真三次FV温度的墙热精度对照。
% 模块用途：用独立三阶体积矩提供温度均值，用面积积分提供墙热参考，
% 检查二次重构在非二次场中的真实误差。MATLAB自行构造二阶几何矩。
% 旧线性NS只在常rho/u人工场中作为空间Fourier基线；不推进或改生产场。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,folder=fullfile(root,'outputs','v60_validation','library_cubic');end
if ~isfolder(folder),mkdir(folder);end
maxNumCompThreads(1);addpath(fullfile(root,'internalFlow'));setupCutCellFlow;
v58=validationFolder(root,'v58');v59=validationFolder(root,'v59');v60=validationFolder(root,'v60');
names={'actual','sphere_grid8','sphere_grid12','sphere_grid16'};checks=cell(4,1);evidence=cell(4,1);
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);cv=gas.R/(gas.gamma-1);
for j=1:4
    if j==1
        a=load(fullfile(v58,'cut_cell_transport_draft/actual_aggregation_results.mat'),'group');
        d=load(fullfile(v58,'cut_cell_draft/cut_cell_geometry_results.mat'),'actual');source=d.actual;spacing=.0625;
    else
        a=load(fullfile(v59,'geometry',[names{j},'.mat']),'source','group','identity');
        source=a.source;spacing=a.identity.spacing;
    end
    c=a.group;d=load(fullfile(v60,'reconstruction_review',[names{j},'_cubic_oracle.mat']),'cubicOracle');o=d.cubicOracle;
    assert(max(abs(o.volume(:)-c.cells.volume)./c.cells.volume)<1e-10 && max(abs(o.sourceCentroid-c.cells.centroid),[],'all')<1e-11);
    [m,details]=buildCutCellVolumeMoments(source,c); %#ok<ASGLU>
    % 独立oracle只提供测试精确均值/面积积分，不是运行几何或重构的依赖。
    reconstructedQuadraticMean=o.T0+c.cells.centroid*o.g(:)+.5*sum((c.cells.centroid*o.H).*c.cells.centroid,2);
    for id=1:numel(m.volume)
        reconstructedQuadraticMean(id)=reconstructedQuadraticMean(id)+.5*sum(o.H.*reshape(m.covariance(id,:,:),3,3),'all');
    end
    assert(max(abs(reconstructedQuadraticMean-o.quadraticTemperatureVolumeAverage(:)))<1e-10);
    l=o.l(:);g=o.g(:);H=o.H;T0=o.T0;alpha=o.a;k=o.conductivity;
    pointFunction=@(x,id)T0+x*g+.5*sum((x*H).*x,2)+alpha*(x*l).^3;
    bc=allBoundaries(struct('type','temperaturePoint','value',pointFunction));
    field=struct('values',o.temperatureVolumeAverage(:),'valueKind','volume-average-temperature','units','K');
    stencil=buildQuadraticTemperatureStencil(c,m,bc);out=quadraticTemperatureFourierFlux(field,c,k,bc,stencil);
    oldBC=allBoundaries(struct('type','wall','temperature',pointFunction,'velocity',[0 0 0]));
    Q=[ones(numel(m.volume),1),zeros(numel(m.volume),3),cv*field.values];
    [~,old]=cutCellNavierStokesResidual(Q,c,gas,oldBC);
    qref=o.wallHeatFluxAreaAverage(:);A=c.walls.area;V=c.cells.volume;
    newError=out.wallHeatFluxIntoWall-qref;oldError=old.wallHeatFluxIntoWall-qref;
    newRMS=rmsWeighted(newError,A);oldRMS=rmsWeighted(oldError,A);refRMS=rmsWeighted(qref,A);
    exactDivergence=k*(trace(H)+6*alpha*sum(l.^2)*(c.cells.centroid*l));
    divergenceRMS=rmsWeighted(out.thermalDiffusionRate-exactDivergence,V);
    cellPointError=out.reconstruction.gradient-o.gradientAtSourceCentroid;
    cellMeanError=out.reconstruction.gradient-o.gradientVolumeAverage;
    oldGradient=reshape(old.gradient(4,:,:),3,[])';
    centroidQref=-k*sum((o.g+c.walls.centroid*H+3*alpha*(c.walls.centroid*l).^2*l').*c.walls.normal,2);
    referenceCentroidDifference=max(abs(qref-centroidQref));
    assert(all(isfinite(newError))&&isfinite(newRMS)&&all(isfinite(out.thermalDiffusionRate)));
    assert(abs(out.thermalBalance)<1e-8,'CubicStudy:Conservation','共享热通量的全域热账不平。');
    assert(~out.physicalSolution&&~out.eligibleForSurfaceCoverage&&out.newSolvedSurfaceAreaM2==0);
    check=struct('caseName',names{j},'spacing',spacing,'fluidGroups',numel(V),'wallPatches',numel(A), ...
        'quadraticWallHeatAreaWeightedRMS',newRMS,'linearWallHeatAreaWeightedRMS',oldRMS, ...
        'wallHeatReferenceRMS',refRMS,'quadraticWallHeatRelativeRMS',newRMS/refRMS, ...
        'linearWallHeatRelativeRMS',oldRMS/refRMS,'quadraticToLinearRMSRatio',newRMS/oldRMS, ...
        'quadraticWallHeatMaximumError',max(abs(newError)),'linearWallHeatMaximumError',max(abs(oldError)), ...
        'quadraticGradientPointVolumeWeightedRMS',rmsWeighted(sqrt(sum(cellPointError.^2,2)),V), ...
        'quadraticGradientMeanVolumeWeightedRMS',rmsWeighted(sqrt(sum(cellMeanError.^2,2)),V), ...
        'linearGradientPointVolumeWeightedRMS',rmsWeighted(sqrt(sum((oldGradient-o.gradientAtSourceCentroid).^2,2)),V), ...
        'thermalDivergenceVolumeWeightedRMS',divergenceRMS,'thermalBalance',out.thermalBalance, ...
        'wallHeatPowerError',sum(A.*newError),'linearWallHeatPowerError',sum(A.*oldError), ...
        'areaMeanVersusCentroidReferenceMaximumDifference',referenceCentroidDifference, ...
        'maximumNeighborRing',max(stencil.ringCount),'minimumScaledRelativeSingularValue',min(stencil.scaledRelativeSingularValue));
    checks{j}=check;evidence{j}=struct('caseName',names{j},'meanField',field,'moments',m, ...
        'output',out,'stencil',stencil,'linearOutput',old,'exactWallHeatAreaAverage',qref, ...
        'exactThermalDivergenceVolumeAverage',exactDivergence,'check',check);
    fprintf('CUBIC_FV_RESULT %s linearRMS=%.9g quadraticRMS=%.9g relative=%.9g ratio=%.9g divergence=%.9g\n', ...
        names{j},oldRMS,newRMS,newRMS/refRMS,newRMS/oldRMS,divergenceRMS);
end
h=cellfun(@(s)s.spacing,checks(2:4));newErrors=cellfun(@(s)s.quadraticWallHeatAreaWeightedRMS,checks(2:4));oldErrors=cellfun(@(s)s.linearWallHeatAreaWeightedRMS,checks(2:4));
newOrder=log(newErrors(1:2)./newErrors(2:3))./log(h(1:2)./h(2:3));oldOrder=log(oldErrors(1:2)./oldErrors(2:3))./log(h(1:2)./h(2:3));
report=struct('status','tested-cubic-true-FV-scalar-spatial-study','checks',{checks}, ...
    'quadraticFixedGeometryObservedOrders',newOrder,'linearFixedGeometryObservedOrders',oldOrder, ...
    'referenceIsWallAreaAverage',true,'inputIsTrueTemperatureVolumeAverage',true, ...
    'geometryMomentsComputedByMATLAB',true,'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'fullNavierStokesStabilityProven',false,'newSolvedSurfaceAreaM2',0, ...
    'scope','三次制造场的标量Fourier空间误差；不是完整NS闭合、稳定性、飞行热流或实验精度');
save(fullfile(folder,'cubic_temperature_evidence.mat'),'report','evidence','-v7');
fid=fopen(fullfile(folder,'cubic_temperature_summary.json'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('CUBIC_FV_STUDY_PASS quadraticOrders=%s linearOrders=%s physicalArea=0\n',mat2str(newOrder,8),mat2str(oldOrder,8));
end
function p=validationFolder(root,version)
p=fullfile(root,'outputs',[version,'_validation']);
switch version
    case 'v58',required=fullfile('cut_cell_transport_draft','actual_aggregation_results.mat');
    case 'v59',required=fullfile('geometry','sphere_grid8.mat');
    case 'v60',required=fullfile('reconstruction_review','actual_cubic_oracle.mat');
end
if ~isfile(fullfile(p,required)),p=fullfile(root,'validation',version);end
end
function b=allBoundaries(p),b=struct('xmin',p,'xmax',p,'ymin',p,'ymax',p,'zmin',p,'zmax',p,'solidWall',p);end
function value=rmsWeighted(x,w),value=sqrt(sum(w.*x.^2)/sum(w));end
