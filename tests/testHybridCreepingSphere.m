function testHybridCreepingSphere(folder)
%TESTHYBRIDCREEPINGSPHERE 运行有双向穿流的球外混合算例并保存逐壁面、界面和网格对照。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(isfolder(folder)&&~isfolder(fullfile(folder,'fine')),'RAH:Sphere:Output','须使用新的运行目录。');
c=struct('radiusM',.001,'centerM',[0 0 0],'velocityMps',[-2e-5 0 0], ...
    'densityKgM3',950,'viscosityPaS',5,'conductivityWmK',.15,'specificHeatJkgK',1500, ...
    'ambientTemperatureK',300,'ambientPressurePa',101325);
[mesh,meta]=buildCanonicalBody('sphere',struct('radius',c.radiusM,'resolution',64));
names={'coarse','medium','fine','moved','conduction','isothermal'};
sizes=[12 16;24 32;48 64;48 64;12 16;12 16];cases=struct;
for k=1:numel(names)
    % 模块一：同一原面上改变体网格或人工切面，零速和等温另作明确的方程控制。
    name=names{k};out=fullfile(folder,name);mkdir(out);
    runCondition=c;tw=302;cut=2;
    if strcmp(name,'moved'),cut=3;end
    if strcmp(name,'conduction'),runCondition.velocityMps(:)=0;end
    if strcmp(name,'isothermal'),tw=300;end
    opt=struct('innerCells',sizes(k,1),'outerCells',sizes(k,1),'angularCells',sizes(k,2), ...
        'cutRatio',cut,'bathRatio',8,'saveFolder',out);
    cfg=createCreepingSphereHybridCase(mesh,c,302,opt);
    cfg.hybridSurface.condition=runCondition;cfg.hybridSurface.condition.thermalBathRadiusM=8*c.radiusM;
    cfg.hybridSurface.wallTemperatureK(:)=tw;
    [model,~,report]=runAeroHeatSingleCase(mesh,cfg);
    assert(report.fullSurfaceComputed && report.engineeringFirstAreaM2==0 && report.regionalAddedAreaM2==report.totalArea);
    assert(numel(report.attempts)==1 && report.attempts.accepted && report.faceCount==size(mesh.F,1));
    assert(isequal(model.hybridWall.faceID,(0:report.faceCount-1)') && isequal(model.F,mesh.F));
    assert(all(model.hybridWall.nearWallCellTemperatureK>=min(tw,c.ambientTemperatureK) & ...
        model.hybridWall.nearWallCellTemperatureK<=max(tw,c.ambientTemperatureK)));
    source=load(fullfile(out,'terminal.mat'),'flow','audit','thermal','wallQ','ok');
    assert(source.ok);audit=source.audit;
    speed=norm(runCondition.velocityMps);a=c.radiusM;
    centers=(mesh.V(mesh.F(:,1),:)+mesh.V(mesh.F(:,2),:)+mesh.V(mesh.F(:,3),:))/3;
    mu=-centers(:,1)./vecnorm(centers,2,2);
    % 模块二：仅在返回结果后使用公开整域球解核对速度、压力扰动、阻力与两个半球质量。
    s=linspace(1,cut,51);y=deval(source.flow,s);
    exactF=s.^2-1.5*s+.5./s;exactFp=2*s-1.5-.5./s.^2;
    velocityError=max([abs((y(1,:)-exactF)./s.^2),abs((y(2,:)-exactFp)./(2*s))]);
    dpScale=c.viscosityPaS*speed/a;
    deltaP=model.hybridWall.pressurePa-c.ambientPressurePa;
    pressureError=max(abs(deltaP+1.5*dpScale*mu))/max(dpScale,realmin);
    rr=-(y(4,1)-2*y(2,1)+4*y(1,1))/2+2*(y(2,1)-2*y(1,1));
    rt=-y(3,1)/2+y(2,1)-y(1,1);
    dragRatio=NaN;if speed>0,dragRatio=(rr-2*rt)*2/9;end
    exactMass=pi*c.densityKgM3*speed*a^2*(cut^2-1.5*cut+.5/cut);
    massError=max(abs([audit.interfaceMassInKgS,audit.interfaceMassOutKgS]-exactMass))/max(exactMass,realmin);
    assert(velocityError<1e-6 && pressureError<1e-6 && (speed==0||abs(dragRatio-1)<1e-6) && massError<1e-6);
    signedCutHeat=audit.interfaceAdvectedHeatNetOutW+audit.interfaceConductedHeatNetOutW;
    assert(abs(signedCutHeat+audit.curvedWallPowerW)<1e-9*max(abs(audit.curvedWallPowerW),realmin));
    if speed>0 && tw~=300,assert(audit.interfaceAdvectedHeatAbsoluteW>0);end
    if strcmp(name,'conduction')
        expected=-4*pi*c.conductivityWmK*(tw-c.ambientTemperatureK)*a*8/7;
        assert(abs(audit.curvedWallPowerW/expected-1)<1e-9);
    elseif strcmp(name,'isothermal')
        assert(all(model.hybridWall.heatFluxWm2==0) && report.wholeSurfaceHeatPowerW==0);
    end
    assert(tw==300 || all(model.hybridWall.heatFluxWm2<0));
    tableFlux=table(source.thermal.muEdges(1:end-1),source.thermal.muEdges(2:end),source.wallQ, ...
        'VariableNames',{'muLow','muHigh','heatIntoWallWm2'});
    writetable(tableFlux,fullfile(out,'wall_bands.csv'));
    cases.(name)=struct('faces',report.faceCount,'triangleAreaM2',report.totalArea,'exactSphereAreaM2',4*pi*a^2, ...
        'wallTrianglePowerW',report.wholeSurfaceHeatPowerW,'velocityError',velocityError, ...
        'pressurePerturbationError',pressureError,'dragRatio',dragRatio,'hemisphereMassError',massError,'interface',audit);
    if strcmp(name,'fine')
        writetable(model.hybridWall,fullfile(out,'wall.csv'));
        page=exportAerothermalWeb(model,fullfile(out,'workbench'));
        raw=strtrim(fileread(fullfile(fileparts(page),'case-data.js')));prefix='window.ROCKET_CASE=';
        assert(startsWith(raw,prefix) && endsWith(raw,';'));d=jsondecode(raw(numel(prefix)+1:end-1));
        for pair={'pressure','heatFlux'}
            f=d.fields(strcmp({d.fields.key},pair{1}));
            if strcmp(pair{1},'pressure'),values=model.hybridWall.pressurePa;else,values=model.hybridWall.heatFluxWm2;end
            assert(all(f.status==1) && max(abs(1000*f.values-values))<1e-8);
            if strcmp(pair{1},'pressure'),method='Stokes';else,method='Fourier';end
            assert(all(contains(string(f.methodNames),method)));
        end
        assert(contains(d.resultCondition.caption,'300 K') && contains(d.caseNote,'2e-05 m/s'));
        for key={'edgeTemperature','recoveryTemperature','wallSolved'}
            f=d.fields(strcmp({d.fields.key},key{1}));
            assert(isscalar(f) && all(f.status==0) && all(isnan(f.values)));
        end
        f=d.fields(strcmp({d.fields.key},'wallInput'));assert(all(f.status==2 & f.values==302));
        f=d.fields(strcmp({d.fields.key},'nearWallCellTemperature'));assert(all(f.status==1 & abs(f.values-model.hybridWall.nearWallCellTemperatureK)<1e-10));
    end
end
% 模块三：不闭合面、非等温壁和超出低速范围的输入显式拒绝。
bad=mesh;bad.F(1,:)=[];cfg=createCreepingSphereHybridCase(bad,c,302);
expect(@()runAeroHeatSingleCase(bad,cfg),'RAH:Sphere:Geometry');
broken=c;broken.velocityMps=100*c.velocityMps;expect(@()createCreepingSphereHybridCase(mesh,broken,302),'RAH:Sphere:Regime');
tem=repmat(302,size(mesh.F,1),1);tem(1)=301;expect(@()createCreepingSphereHybridCase(mesh,c,tem),'RAH:Sphere:Wall');
cfg=createCreepingSphereHybridCase(mesh,c,302);cfg.hybridSurface.condition.velocityMps=broken.velocityMps;
expect(@()runAeroHeatSingleCase(mesh,cfg),'RAH:Sphere:Regime');
cfg.hybridSurface.condition.velocityMps=c.velocityMps;cfg.hybridSurface.wallTemperatureK=tem;
expect(@()runAeroHeatSingleCase(mesh,cfg),'RAH:Sphere:Wall');
summary=struct('cases',cases,'sourceGeometry',meta.kind,'condition',c,'bathRatio',8,'wallTemperatureK',302, ...
    'experimentalValidation',false,'rocketCoverageGainM2',0,'browserInteractionTested',false);
fid=fopen(fullfile(folder,'result.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));
end

function expect(run,identifier)
try
    run();
catch err
    assert(strcmp(err.identifier,identifier),'意外的拒绝原因：%s',err.message);return;
end
error('RAH:Sphere:MissingRefusal','应拒绝该输入：%s',identifier);
end
