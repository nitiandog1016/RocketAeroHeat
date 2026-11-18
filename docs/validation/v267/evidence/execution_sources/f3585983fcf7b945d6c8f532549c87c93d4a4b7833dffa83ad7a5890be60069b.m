function [report,mesh,cfg]=demoHybridWallCalculation(folder)
%DEMOHYBRIDWALLCALCULATION 用两个独立平行间隙检验工程解与局部 CFD 联合求解。
% 两个方向周期，每个方格的两块恒温壁构成全部湿壁。控制算例使用已验证
% 的人工恒物性气体，目的是用解析解逐面对照；此处没有实验输入。
if nargin<1,folder=fullfile(tempdir,'rah_hybrid_gap');end
if ~isfolder(folder),mkdir(folder);end
assert(~isfile(fullfile(folder,'states.mat')),'RAH:Hybrid:Output','请选择新的输出目录。');
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'internalFlow'));
c=struct('heightM',0.06,'meanDensity',1,'upperSpeedMps',3,'gamma',1.4,'R',1,'mu',0.5,'Pr',0.72);
% 模块一：建四个单位面积的实际矩形壁，三角面上的编号与工作台一致。
mesh=struct('V',zeros(16,3),'F',zeros(8,3));
for w=1:4
    x=2*floor((w-1)/2);y=mod(w-1,2)*c.heightM;
    start=(w-1)*4;
    mesh.V(start+(1:4),:)=[x y 0;x+1 y 0;x+1 y 1;x y 1];
    tri=[1 2 3;1 3 4];if mod(w,2)==0,tri=tri(:,[1 3 2]);end
    mesh.F(2*w-1:2*w,:)=start+tri;
end
tw=[300;300;330;330;300;300;300;300];
regions=[repmat("stationary-gap",4,1);repmat("moving-gap",4,1)];
cfg=struct('analysisMode','hybridSurface','hybridSurface',struct('condition',c,'wallTemperatureK',tw, ...
    'region',regions,'engineeringSolver',@stationaryGap,'localSolver',@movingGap));
[current,~,report]=runAeroHeatSingleCase(mesh,cfg);
report.viewer=exportAerothermalWeb(current,fullfile(folder,'workbench'));
writetable(current.hybridWall,fullfile(folder,'wall.csv'));

    function rows=stationaryGap(gasMesh,condition,wallT)
        % 模块二：静止均质间隙解析导热；由气体质量与温度积分推算压力。
        faceID=find(gasMesh.region=="stationary-gap")-1;
        lo=wallT(faceID(1)+1);hi=wallT(faceID(end)+1);dT=hi-lo;
        if dT==0,meanT=lo;else,meanT=dT/log1p(dT/lo);end
        pressurePa=repmat(condition.meanDensity*condition.R*meanT,4,1);
        conductivity=condition.mu*condition.gamma*condition.R/((condition.gamma-1)*condition.Pr);
        heatFluxWm2=conductivity*dT/condition.heightM*[1;1;-1;-1];
        method=repmat("稳态平壁导热与封闭气体质量积分",4,1);
        rows=flipud(table(faceID,pressurePa,heatFluxWm2,method));
    end

    function regionResult=movingGap(gasMesh,condition,wallT,regionIDs,targetIDs,~)
        % 模块三：缺项触发原黏性内核，壁面值由最后接受态派生。
        assert(all(gasMesh.region(regionIDs+1)=="moving-gap") && isequal(sort(regionIDs),sort(targetIDs)));
        gas=struct('gamma',condition.gamma,'R',condition.R,'mu',condition.mu,'Pr',condition.Pr,'temperatureRange',[200 500]);
        ny=8;periodic=struct('type','periodic');
        bc=struct('xmin',periodic,'xmax',periodic,'zmin',periodic,'zmax',periodic);
        bc.ymin=struct('type','wall','temperature',wallT(regionIDs(1)+1),'velocity',[0 0 0]);
        bc.ymax=struct('type','wall','temperature',wallT(regionIDs(end)+1),'velocity',[condition.upperSpeedMps 0 0]);
        initial=zeros(1,ny,1,5);initial(:,:,:,1)=condition.meanDensity;
        initial(:,:,:,5)=condition.meanDensity*condition.R*bc.ymin.temperature/(condition.gamma-1);
        options=struct('steadyResidualTolerance',2e-6,'steadyRequiredSteps',8,'minimumSteadySteps',8, ...
            'maxSteps',4000,'convectiveCFL',0.45,'diffusiveCFL',0.45);
        grid=struct('spacing',[1 condition.heightM/ny 1]);
        save(fullfile(folder,'initial.mat'),'initial','grid','bc','gas','options');
        options.acceptedStepObserver=@(step,time,Q)saveGapFrame(folder,step,time,Q);
        result=advanceCompressibleNavierStokes3D(initial,grid,false(1,ny,1),gas,bc,options);
        options=rmfield(options,'acceptedStepObserver');
        save(fullfile(folder,'states.mat'),'result','grid','bc','gas','options');
        accepted=result.steadyConverged && result.globalBalanceRelative<1e-10;
        regionResult=struct('accepted',accepted,'reason',result.stopReason,'wall',table());
        if ~accepted,return;end
        faceID=regionIDs; pressurePa=nan(4,1); heatFluxWm2=nan(4,1);
        y=mean(reshape(gasMesh.V(gasMesh.F(faceID+1,:),2),[],3),2);
        for side=[-1 1]
            w=result.wallFaces(arrayfun(@(item)item.normal(2)==side,result.wallFaces));
            assert(isscalar(w)&&w.area==1);
            if side<0,onWall=y==0;else,onWall=abs(y-condition.heightM)<1e-12;end
            assert(nnz(onWall)==2);
            pressurePa(onWall)=w.pressure; heatFluxWm2(onWall)=w.heatFluxIntoWall;
        end
        method=repmat("局部有限体积 Navier–Stokes 稳态壁通量",4,1);
        regionResult.wall=flipud(table(faceID,pressurePa,heatFluxWm2,method));
    end
end

function saveGapFrame(folder,step,time,Q)
% 保留首步及每 50 个完整接受步的全场，正常终态另存。
if step==1 || mod(step,50)==0
    save(fullfile(folder,sprintf('step%04d.mat',step)),'step','time','Q');
end
end
