function report=studyCurvedWallQuadraticAccuracy(folder)
%STUDYCURVEDWALLQUADRATICACCURACY 真实FV二次能量场的曲面Fourier精度研究。
% 模块用途：以独立边界矩oracle取得精确体积平均；三档固定曲面仅改变CV。
% 不把质心primitive样本当体积平均，不把本制造场当无源NS稳态/飞行实验。
% 原始残差的对流会响应制造压力梯度，故热扩散散度从实际黏性通量独立积分。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,folder=fullfile(root,'outputs','v59_validation','accuracy');end
if ~isfolder(folder),mkdir(folder);end
maxNumCompThreads(1);addpath(fullfile(root,'internalFlow'),'-begin');setupCutCellFlow;
base=fullfile(root,'outputs','v59_validation');
if ~isfile(fullfile(base,'curved_accuracy_oracle','actual_fv_volume_moments.mat')),base=fullfile(root,'validation','v59');end
oracle=fullfile(base,'curved_accuracy_oracle');
prior=fullfile(root,'outputs','v58_validation');
if ~isfile(fullfile(prior,'cut_cell_transport_draft','actual_aggregation_results.mat')),prior=fullfile(root,'validation','v58');end
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
cv=gas.R/(gas.gamma-1);k=gas.mu*gas.gamma*cv/gas.Pr;
g=[8 -5 3];H=[30 4 -2;4 -20 3;-2 3 10];T0=300;
wall=struct('type','wall','temperature',@(x,id)T0+x*g'+.5*sum((x*H).*x,2),'velocity',[0 0 0]);bc=struct;
for name={'xmin','xmax','ymin','ymax','zmin','zmax','solidWall'},bc.(name{1})=wall;end
names={'actual_rocket','sphere_grid8','sphere_grid12','sphere_grid16'};evidence=cell(4,1);statistics=cell(4,1);
for j=1:numel(names)
    if j==1
        a=load(fullfile(prior,'cut_cell_transport_draft','actual_aggregation_results.mat'),'group');
        b=load(fullfile(oracle,'actual_fv_volume_moments.mat'),'moments');spacing=.0625;
    else
        name=names{j};a=load(fullfile(base,'geometry',[name,'.mat']),'group','identity');spacing=a.identity.spacing;
        b=load(fullfile(oracle,[name,'_fv_volume_moments.mat']),'moments');
    end
    c=a.group;m=b.moments;n=numel(c.cells.volume);A=c.walls.area;
    % 独立SciPy MAT可把一维量保存成行向量；只统一数组方向，不改变任何几何/物理值。
    for key={'volume','sourceVolume','quadraticTemperatureVolumeAverage','quadraticTemperatureAtSourceCentroid','quadraticTotalEnergyVolumeAverage'}
        m.(key{1})=m.(key{1})(:);
    end
    assert(isequal(size(m.volume),size(c.cells.volume)) && max(abs(m.volume-c.cells.volume)./c.cells.volume)<1e-10);
    assert(max(abs(m.centroid-c.cells.centroid),[],'all')<1e-11);
    assert(isequal(m.H,H)&&isequal(m.g,g)&&all(m.xref==0)&&m.T0==T0,'CurvedAccuracy:Oracle','制造场身份不一致。');
    % 常rho/u使保守态的体积平均可由二阶矩精确给出；不存在primitive非线性换算歧义。
    Q=[ones(n,1),zeros(n,3),cv*m.quadraticTemperatureVolumeAverage];
    [rate,info]=cutCellNavierStokesResidual(Q,c,gas,bc);
    exactGradient=g+c.walls.centroid*H;exactHeat=-k*sum(exactGradient.*c.walls.normal,2);
    error=info.wallHeatFluxIntoWall-exactHeat;referenceRMS=sqrt(sum(A.*exactHeat.^2)/sum(A));
    heatRMS=sqrt(sum(A.*error.^2)/sum(A));heatL1=sum(A.*abs(error))/sum(A);
    cellGradient=reshape(info.gradient(4,:,:),3,[])';gradientError=cellGradient-m.quadraticGradientVolumeAverage;
    gradientRMS=sqrt(sum(c.cells.volume.*sum(gradientError.^2,2))/sum(c.cells.volume));
    diffusion=thermalDiffusion(info,c);expectedDivergence=k*trace(H);
    divergenceRMS=sqrt(sum(c.cells.volume.*(diffusion-expectedDivergence).^2)/sum(c.cells.volume));
    totalHeat=sum(A.*info.wallHeatFluxIntoWall);exactTotalHeat=sum(A.*exactHeat);
    budget=sum(c.cells.volume.*diffusion)+info.netBoundaryFlux(5);
    % 外边界全部是静止墙且初态u=0：净边界能量只有给定温度下的Fourier通量。
    assert(abs(budget)<1e-8*max(1,abs(info.netBoundaryFlux(5))),'CurvedAccuracy:Budget','扩散积分与实际墙边界账不一致。');
    % 另存旧“质心点值”输入，仅量化其与真实FV输入的差异，不能替换主结果。
    centroidQ=[ones(n,1),zeros(n,3),cv*m.quadraticTemperatureAtSourceCentroid];
    [~,centroidInfo]=cutCellNavierStokesResidual(centroidQ,c,gas,bc);
    centroidError=centroidInfo.wallHeatFluxIntoWall-exactHeat;
    stat=struct('caseName',names{j},'spacing',spacing,'fluidGroups',n,'wallPatches',numel(A), ...
        'wallHeatAreaWeightedRMS',heatRMS,'wallHeatAreaWeightedL1',heatL1, ...
        'wallHeatMaximumError',max(abs(error)),'wallHeatReferenceRMS',referenceRMS, ...
        'wallHeatRMSRelative',heatRMS/referenceRMS,'cellGradientVolumeWeightedRMS',gradientRMS, ...
        'thermalDivergenceVolumeWeightedRMS',divergenceRMS,'expectedThermalDivergence',expectedDivergence, ...
        'wallHeatPower',totalHeat,'exactWallHeatPower',exactTotalHeat, ...
        'wallHeatPowerError',totalHeat-exactTotalHeat,'globalThermalBudget',budget, ...
        'fvTemperatureMeanVersusCentroidMaximumDifference',max(abs(m.quadraticTemperatureVolumeAverage-m.quadraticTemperatureAtSourceCentroid)), ...
        'centroidInputWallHeatRMS',sqrt(sum(A.*centroidError.^2)/sum(A)));
    statistics{j}=stat;evidence{j}=struct('caseName',names{j},'Q',Q,'rate',rate,'info',info, ...
        'exactWallHeat',exactHeat,'thermalDiffusionRate',diffusion,'statistics',stat);
    fprintf('QUADRATIC_CURVED_RESULT %s h=%.6g heatRMS=%.9g rel=%.6g max=%.9g heatPowerError=%.9g\n', ...
        names{j},spacing,heatRMS,stat.wallHeatRMSRelative,stat.wallHeatMaximumError,stat.wallHeatPowerError);
end
h=cellfun(@(s)s.spacing,statistics(2:4));errors=cellfun(@(s)s.wallHeatAreaWeightedRMS,statistics(2:4));
order=log(errors(1:2)./errors(2:3))./log(h(1:2)./h(2:3));
report=struct('status','completed-manufactured-curved-accuracy-study','statistics',{statistics}, ...
    'fixedSurfaceWallHeatRMSOrders',order,'finestWallHeatRMSBelowCoarsest',errors(3)<errors(1), ...
    'quadraticInputsAreTrueConservativeVolumeAverages',true,'centroidPrimitiveSamplesUsedAsVolumeAverages',false, ...
    'geometryFixedAcrossGrids',true,'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'newSolvedSurfaceAreaM2',0,'testGas',gas,'T0',T0,'g',g,'H',H, ...
    'scope','固定320面曲面及实际箭体粗网格；精确二次FV能量均值、空间Fourier壁热/梯度/散度误差，不是无源NS稳态或实验验证');
save(fullfile(folder,'curved_quadratic_evidence.mat'),'report','evidence','-v7');
fid=fopen(fullfile(folder,'curved_quadratic_summary.json'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('CURVED_QUADRATIC_STUDY_COMPLETED rmsOrders=%s finestBelowCoarsest=%d newSolvedArea=0\n',mat2str(order,7),report.finestWallHeatRMSBelowCoarsest);
end
function rate=thermalDiffusion(info,c)
% 模块用途：从实际面通量独立组装纯热扩散，不用有压力梯度的总能量残差冒充导热方程。
n=numel(c.cells.volume);f=c.openFaces;w=c.walls;integral=zeros(n,1);
for k=find(f.area>0)'
    flux=f.area(k)*info.viscousFaceFlux(k,5);m=f.cellMinus(k);p=f.cellPlus(k);
    if m>0,integral(m)=integral(m)+flux;end
    if p>0,integral(p)=integral(p)-flux;end
end
for k=find(w.ownerCell>0)',id=w.ownerCell(k);integral(id)=integral(id)-w.area(k)*info.wallHeatFluxIntoWall(k);end
rate=integral./c.cells.volume;
end
