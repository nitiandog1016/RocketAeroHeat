function report=prepareActualRocketFullDomain(out)
%PREPAREACTUALROCKETFULLDOMAIN 原19752面/精确30km来流准备完整混合域及正守恒初场。
% 模块一：原基线V/F、来流向量及给定壁温逐项读取，拒绝用旧近似导入工况。
% 模块二：近壁厚度按参考状态平板尺度估计，仅用于网格/初猜，不用于壁热预测。
% 模块三：18层保留首6层高度；仍须独立网格与后续物理验证，不改变生产默认。
assert(~isfolder(out));mkdir(out);
source='outputs/v56_validation/state_validity/state_validity_actual_rocket.mat';
saved=load(source,'result','fs','cfg');s=struct('V',saved.result.V,'F',double(saved.result.F));
d=load('outputs/v113_validation/actual_input.mat','V','F');
assert(isequal(s.V,d.V) && isequal(s.F,d.F) && size(s.F,1)==19752);
fs=saved.fs;Tw=saved.cfg.flight.Twall;
assert(fs.Mach==5 && fs.alpha==20 && fs.beta==0 && fs.altitude==30000 && Tw==300);
L=max(s.V(:,1))-min(s.V(:,1));Pr=.72;
Taw=fs.temperature*(1+sqrt(Pr)*(fs.gamma-1)/2*fs.Mach^2);
Tref=fs.temperature+.5*(Tw-fs.temperature)+.22*(Taw-fs.temperature);
mu=1.458e-6*Tref^1.5/(Tref+110.4);rho=fs.pressure/(fs.R*Tref);
thickness=5*sqrt(mu*L/(rho*fs.speed));
heights=[0 cumsum(2e-5*1.4.^(0:17))];assert(heights(end)>1.5*thickness);
opt=struct('layerHeights',heights,'downstreamLengths',2,'upstreamLengths',1,'lateralLengths',1, ...
    'outerMinSize',.02,'farSize',.35,'pythonExecutable','D:/anacondapython/python.exe','maximumSeconds',120);
folder=fullfile(out,'mesh');volume=buildSU2HybridVolume(s,folder,opt);
flight=struct('Mach',fs.Mach,'temperature',fs.temperature,'pressure',fs.pressure, ...
    'wallTemperature',Tw,'alpha',fs.alpha,'beta',fs.beta,'iterations',600,'outputEvery',600,'namedCheckpoints',true);
engine=fullfile(out,'engine');manifest=writeSU2RocketCase(volume,engine,flight);
assert(norm(manifest.expectedVelocity-fs.flowVelocityVector(:)')<1e-10);
initial=writeSU2HybridInitialState(volume,engine,thickness);
manifest=jsondecode(fileread(fullfile(engine,'case_manifest.json')));
projection=buildSU2WallProjection(volume,s);save(fullfile(engine,'projection.mat'),'projection','manifest','-v7');
save(fullfile(out,'surface.mat'),'s','fs','flight','-v7');
report=struct('baselineSource',source,'originalFaces',size(s.F,1),'bodyArea',volume.bodyArea, ...
    'bodyLengthM',L,'ReL',fs.reynoldsPerMeter*L,'referenceTemperatureK',Tref, ...
    'laminarThicknessScaleM',thickness,'collarThicknessM',heights(end),'layers',18, ...
    'firstSixHeightsRetained',true,'flight',flight,'initialization',initial, ...
    'meshReport',volume.meshConstruction,'topology',volume.topology, ...
    'flowAssumption','laminar numerical starting study; transition/RANS and thermal accuracy remain unqualified', ...
    'newQualifiedArea',0,'physicsQualified',false,'productionDefaultChanged',false);
fid=fopen(fullfile(out,'preparation.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Exact real 20-degree full domain prepared; 18 layers and conservative initial guess. No solver started.\n');
end
