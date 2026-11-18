function testHybridSharedFluid(folder)
%TESTHYBRIDSHAREDFLUID 在一个低速库埃特流体域上分别求两个实体壁的热流。
% 平行板充分发展，法向无穿流。截断端传递切向速度、温度与区域质量，
% 对照返回的切向牵引、法向应力和总能流。中间切面只参与守恒核对。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'internalFlow'));
if ~isfolder(folder),mkdir(folder);end
assert(~isfile(fullfile(folder,'initial.mat')),'RAH:Hybrid:Output','需要新的求解目录。');
c=struct('height',.06,'cut',.03,'speed',3,'rhoMean',1,'R',1,'gamma',1.4,'mu',.5,'Pr',.72);
quad=[0 0 0;1 0 0;1 0 1;0 0 1];
mesh=struct('V',[quad;quad+[0 c.height 0]],'F',[1 2 3;1 3 4;5 7 6;5 8 7]);
tw=[300;300;301;301];k=c.mu*c.gamma*c.R/((c.gamma-1)*c.Pr);
cfg=struct('analysisMode','hybridSurface','hybridSurface',struct('condition',c,'wallTemperatureK',tw, ...
    'region',["lower";"lower";"upper";"upper"],'engineeringSolver',@lowerGap,'localSolver',@upperGap));
[model,~,report]=runAeroHeatSingleCase(mesh,cfg);
assert(report.fullSurfaceComputed && report.totalArea==2 && report.engineeringFirstAreaM2==1 && report.regionalAddedAreaM2==1);
% 整隙解析对照只在求解结束后计算，从实体壁及整体气体质量独立复算。
etaT=@(z)tw(1)+(tw(end)-tw(1))*z+c.mu*c.speed^2/(2*k).*z.*(1-z);
pExact=c.rhoMean*c.R/integral(@(z)1./etaT(z),0,1);
wallPower=c.mu*c.speed^2/c.height;
conduction=k*(tw(end)-tw(1))/c.height;
qExact=[conduction+wallPower/2;conduction+wallPower/2;-conduction+wallPower/2;-conduction+wallPower/2];
w=model.hybridWall;
pError=max(abs(w.pressurePa/pExact-1));qError=max(abs(w.heatFluxWm2./qExact-1));
assert(pError<1e-4 && qError<1e-5 && abs(report.wholeSurfaceHeatPowerW/wallPower-1)<1e-5);
page=exportAerothermalWeb(model,fullfile(folder,'workbench'));
src=strtrim(fileread(fullfile(fileparts(page),'case-data.js')));prefix='window.ROCKET_CASE=';
assert(startsWith(src,prefix)&&endsWith(src,';')); d=jsondecode(src(numel(prefix)+1:end-1));
for key={'heatFlux','pressure'}
    field=d.fields(strcmp({d.fields.key},key{1}));
    if strcmp(key{1},'pressure'),values=w.pressurePa;else,values=w.heatFluxWm2;end
    assert(all(field.status==1) && all(abs(1000*field.values-values)<1e-9));
end
writetable(w,fullfile(folder,'wall.csv'));
a=load(fullfile(folder,'states.mat'),'ns','interfaceCheck');
summary=struct('scope','connected fully developed periodic low-speed Couette control', ...
    'physicalWallAreaM2',report.totalArea,'engineeringAreaM2',report.engineeringFirstAreaM2, ...
    'localCFDAreaM2',report.regionalAddedAreaM2,'jointComputedAreaM2',report.jointComputedAreaM2, ...
    'fullSurfaceComputed',report.fullSurfaceComputed,'internalCutCountedInWallArea',false, ...
    'pressureMaxRelativeError',pError,'heatMaxRelativeError',qError,'totalHeatW',report.wholeSurfaceHeatPowerW, ...
    'wallMechanicalWorkW',wallPower,'acceptedSteps',a.ns.acceptedSteps,'lastEightResidual',a.ns.history.residualNorm(end-7:end), ...
    'globalBalanceRelative',a.ns.globalBalanceRelative,'interfaceCheck',a.interfaceCheck, ...
    'experimentalValidation',false,'rocketCoverageGainM2',0);
fid=fopen(fullfile(folder,'result.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));

    function wall=lowerGap(~,condition,wallT)
        % 模块一：用两个半间隙的剪切阻力、导热及黏性功计算界面状态。
        L=[condition.cut,condition.height-condition.cut]; lo=wallT(1);hi=wallT(end);
        tau=condition.speed/sum(L./condition.mu); ui=tau*L(1)/condition.mu;
        diss=condition.mu*[ui,condition.speed-ui].^2/(2*k);
        ti=((lo+diss(1))/L(1)+(hi+diss(2))/L(2))/sum(1./L);
        lowInvT=integral(@(x)1./(lo+(ti-lo)*x+diss(1)*x.*(1-x)),0,1);
        highInvT=integral(@(x)1./(ti+(hi-ti)*x+diss(2)*x.*(1-x)),0,1);
        pressure=condition.rhoMean*condition.R*condition.height/(L*[lowInvT;highInvT]);
        lowHeat=k*(ti-lo+diss(1))/L(1);cutHeat=k*(ti-lo-diss(1))/L(1);
        wall=table([0;1],[pressure;pressure],[lowHeat;lowHeat],repmat("平行剪切层质量与能量积分",2,1), ...
            'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
        wall.Properties.UserData=struct('height',L(2),'speed',ui,'temperature',ti,'pressure',pressure, ...
            'meanDensity',pressure*highInvT/condition.R,'shearPa',tau,'heatWm2',cutHeat, ...
            'energyWm2',cutHeat+ui*tau,'targetRegion',"upper");
    end

    function local=upperGap(g,condition,wallT,regionIDs,targets,engineering)
        % 模块二：从工程计算取得中间切面条件，局部域使用原 NS/Fourier 算子。
        cut=engineering.Properties.UserData;
        assert(all(g.region(regionIDs+1)==cut.targetRegion) && isequal(regionIDs,targets));
        n=8;gas=struct('gamma',condition.gamma,'R',condition.R,'mu',condition.mu,'Pr',condition.Pr,'temperatureRange',[200 500]);
        periodic=struct('type','periodic');bound=struct('xmin',periodic,'xmax',periodic,'zmin',periodic,'zmax',periodic);
        bound.ymin=struct('type','wall','temperature',cut.temperature,'velocity',[cut.speed 0 0]);
        bound.ymax=struct('type','wall','temperature',wallT(targets(1)+1),'velocity',[condition.speed 0 0]);
        grid=struct('spacing',[1 cut.height/n 1]);
        initial=zeros(1,n,1,5);initial(:,:,:,1)=cut.meanDensity;
        initial(:,:,:,5)=cut.meanDensity*gas.R*cut.temperature/(gas.gamma-1);
        opts=struct('steadyResidualTolerance',2e-6,'steadyRequiredSteps',8,'minimumSteadySteps',8, ...
            'maxSteps',4000,'convectiveCFL',.45,'diffusiveCFL',.45);
        save(fullfile(folder,'initial.mat'),'initial','gas','grid','bound','cut','opts');
        opts.acceptedStepObserver=@keepFrame;
        ns=advanceCompressibleNavierStokes3D(initial,grid,false(1,n,1),gas,bound,opts);
        cutFace=ns.wallFaces(arrayfun(@(f)f.normal(2)<0,ns.wallFaces));assert(isscalar(cutFace));
        surface=ns.wallFaces(arrayfun(@(f)f.normal(2)>0,ns.wallFaces));assert(isscalar(surface));
        interfaceCheck=struct('pressureRelative',abs(cutFace.pressure/cut.pressure-1), ...
            'shearRelative',abs(cutFace.traction(1)/cut.shearPa-1),'heatRelative',abs(cutFace.heatFluxIntoWall/cut.heatWm2-1), ...
            'totalEnergyRelative',abs(cutFace.totalEnergyFluxOut/cut.energyWm2-1), ...
            'normalVelocityMps',dot(cutFace.normal,cutFace.velocityInput), ...
            'normalMassFlux',ns.finalSpatialInfo.faceFlux{2}(1,1));
        ok=ns.steadyConverged && ns.globalBalanceRelative<1e-10 && interfaceCheck.pressureRelative<1e-4 && ...
            interfaceCheck.shearRelative<1e-5 && interfaceCheck.heatRelative<1e-5 && interfaceCheck.totalEnergyRelative<1e-5 && ...
            interfaceCheck.normalVelocityMps==0 && interfaceCheck.normalMassFlux==0;
        save(fullfile(folder,'states.mat'),'ns','interfaceCheck');
        local=struct('accepted',ok,'reason','shared fluid: steady, interface and balance check','wall',table());
        if ok
            local.wall=table(targets,repmat(surface.pressure,numel(targets),1),repmat(surface.heatFluxIntoWall,numel(targets),1), ...
                repmat("局部有限体积 Navier–Stokes 壁通量",numel(targets),1),'VariableNames',engineering.Properties.VariableNames);
        end
    end

    function keepFrame(step,time,Q)
        if step==1 || mod(step,50)==0,save(fullfile(folder,sprintf('step%04d.mat',step)),'step','time','Q');end
    end
end
