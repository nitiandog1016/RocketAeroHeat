function report=testTmrAxisymmetricStartup(folder)
%TESTTMRAXISYMMETRICSTARTUP 在公开原网格内真实推进20步Mach7区域NS。
% 模块一：只用WIND公开计算边界的M/p0/T0构造内部均匀初场，无外部CFD解。
% 模块二：常gamma、冻结入口Sutherland值的常mu层流原型，不冒充SA/SST。
% 模块三：保持等温无滑移壁和显式特征储库，真实RK2更新并核对四量账。
% 20步仅为初始化/组装/开放边界短推进校核，未达到穿越时间或稳态；
% 初始墙比气流温暖，signed q可能为负，不能取绝对值冒充气动加热。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'),fullfile(root,'tests'));
if nargin<1,folder=fullfile(root,'outputs','v62_validation','tmr_startup');end
if ~exist(folder,'dir'),mkdir(folder);end
benchmark=loadTmrMach7Benchmark;source=benchmark.contract.computational_truncated_inflow;
grid=loadTmrMach7Grid(4);mesh=buildAxisymmetricFVGeometry(grid.X,grid.R);
gamma=1.4;Rgas=287.05;factor=1+(gamma-1)*source.Mach^2/2;
temperature=source.total_temperature_K_WIND_input/factor;
pressure=source.total_pressure_Pa_WIND_input/factor^(gamma/(gamma-1));
velocity=-source.Mach*sqrt(gamma*Rgas*temperature);rho=pressure/(Rgas*temperature);
mu=1.716e-5*(temperature/273.15)^1.5*(273.15+110.4)/(temperature+110.4);
gas=struct('gamma',gamma,'R',Rgas,'mu',mu,'Pr',.72,'temperatureRange',[40 2000]);
primitive=[rho velocity 0 temperature];Q=axisymmetricFVConservative(repmat(primitive,mesh.cellCount,1),gas);
far=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',[rho velocity 0 0 temperature]);
boundary=struct('iMin',far,'iMax',far,'jMin',struct('type','wall','temperature',311,'velocity',[0 0]),'jMax',far);
options=struct('maxSteps',20,'convectiveCFL',.25,'diffusiveCFL',.25,'steadyResidualTolerance',0);
timer=tic;result=advanceAxisymmetricNavierStokes(Q,mesh,gas,boundary,options);seconds=toc(timer);
assert(result.acceptedSteps==20&&~result.steadyConverged&&~result.completedRequestedTime&& ...
    strcmp(result.stopReason,'maximum-steps-reached')&&~isequal(result.Q,Q));
assert(result.globalBalanceRelative<1e-10&&all(isfinite(result.Q),'all'));
P=axisymmetricFVState(result.Q,gas);wall=result.finalInfo.wall;
assert(all(wall.temperatureInput==311)&&all(isfinite(wall.heatFluxIntoWall))&&all(P(:,1)>0));
crossing=(max(grid.X(:,1))-min(grid.X(:,1)))/abs(velocity);
assert(result.time<.01*crossing,'短推进测试不得改成未审查的长段计算。');
report=struct('status','passed-unconverged-public-grid-Mach7-internal-startup-only', ...
    'Mach',source.Mach,'prescribedUniformPrimitive',primitive,'gas',gas,'sourceGrid',grid.source, ...
    'sourceBoundaryIdentity','Published WIND total M/p0/T0; characteristic reservoir uses that prescribed uniform state', ...
    'transportIdentity','Constant mu frozen at prescribed inlet Sutherland value, constantPr.72/caloricgamma1.4; not turbulent TMR air validation', ...
    'sourceBoundaryDifference','All open faces use explicit zero-swirl characteristic Riemann reservoirs; no claimed exact WIND/NMF turbulent boundary reproduction', ...
    'cellCount',mesh.cellCount,'acceptedSteps',result.acceptedSteps,'seconds',seconds,'time',result.time, ...
    'nominalCrossingTime',crossing,'timeOverCrossing',result.time/crossing, ...
    'residualNorm',result.residualNorm,'minimumTemperature',min(P(:,4)),'maximumTemperature',max(P(:,4)), ...
    'wallHeatFluxSignedRange',[min(wall.heatFluxIntoWall),max(wall.heatFluxIntoWall)], ...
    'globalBalanceRelative',result.globalBalanceRelative,'steadyConverged',false, ...
    'experimentalValidation',false,'physicalSolution',false,'eligibleForSurfaceCoverage',false,'newSolvedAreaM2',0);
save(fullfile(folder,'tmr_startup_evidence.mat'),'result','mesh','gas','boundary','options','report','primitive');
fid=fopen(fullfile(folder,'tmr_startup_summary.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('TMR actual-grid Mach%.4g startup: %d real steps, %.9g s / crossing %.9g s, balance%.3g; NOT steady/experimental.\n', ...
    source.Mach,result.acceptedSteps,result.time,crossing,result.globalBalanceRelative);
end
