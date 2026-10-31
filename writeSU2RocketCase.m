function info=writeSU2RocketCase(volume,folder,flight)
%WRITESU2ROCKETCASE 导出完整三维网格与明确物性的SU2配置，不启动求解。
% 模块一：保留棱柱/六面体/四面体及真实墙面、完整远场，输入须已通过几何审核。
% 模块二：项目来流[-cosαcosβ,-sinβ,+sinαcosβ]转换到SU2的攻角/侧滑。
% 模块三：给定Tw与流体温度分开；几何预处理不授予热流资格。
if nargin<3,flight=struct;end
defaults=struct('Mach',4,'temperature',300,'pressure',5000,'wallTemperature',300, ...
    'alpha',0,'beta',0,'iterations',1,'outputEvery',1000,'namedCheckpoints',false,'thermalBoundary',struct([]));
assert(all(ismember(fieldnames(flight),fieldnames(defaults))));
for key=fieldnames(defaults)',k=key{1};if ~isfield(flight,k),flight.(k)=defaults.(k);end,end
for key={'Mach','temperature','pressure','wallTemperature','iterations'}
    validateattributes(flight.(key{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for key={'alpha','beta'},validateattributes(flight.(key{1}),{'numeric'},{'scalar','real','finite','>',-90,'<',90});end
assert(flight.Mach>1 && flight.iterations==fix(flight.iterations));
validateattributes(flight.outputEvery,{'numeric'},{'scalar','integer','positive'});
validateattributes(flight.namedCheckpoints,{'logical'},{'scalar'});
% 具名检查点会保留多个时刻，须由调用方显式启用并限定输出预算；默认不扩大存储。
overwrite='YES';if flight.namedCheckpoints,overwrite='NO';end
assert(strcmp(volume.kind,'complete-rocket-three-dimensional-su2-volume') && ~volume.physicalSolution);
% 热边界在写目录前验证。空间Tw必须先有实际调用器，禁止静默降级为统一Tw。
spec=flight.thermalBoundary;
if isempty(spec),spec=struct('mode','isothermal','temperature',flight.wallTemperature);end
thermal=resolveSU2ThermalBoundary(struct('V',volume.sourceVertices,'F',volume.sourceFaces),spec);
spatial=strcmp(thermal.mode,'axialTemperature');
if spatial
    runtimeFile=fullfile(fileparts(mfilename('fullpath')),'tools','runtime','su2','8.6.0-rah-axial-wall-v1','runtime.json');
    assert(isfile(runtimeFile),'RocketAeroHeat:SU2:SpatialThermalAdapter', ...
        '尚未构建原生空间壁温运行时；不能退回官方包的统一壁温。');
    runtime=jsondecode(fileread(runtimeFile));assert(strcmp(runtime.flavor,thermal.requiredRuntimeFlavor));
    thermalLine=sprintf('MARKER_ISOTHERMAL= (wall, %.17g)',spec.temperature(1));
elseif strcmp(thermal.mode,'adiabatic'),thermalLine='MARKER_HEATFLUX= (wall, 0)';
else,thermalLine=sprintf('MARKER_ISOTHERMAL= (wall, %.17g)',spec.temperature);end
assert(~exist(folder,'dir'),'RocketAeroHeat:SU2:Existing','不能覆写已有计算目录。');mkdir(folder);
fid=fopen(fullfile(folder,'mesh.su2'),'w');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'NDIME= 3\nNELEM= %d\n',numel(volume.cells));
for k=1:numel(volume.cells)
    cellNodes=volume.cells{k};
    % SU2棱柱端面次序与本项目正四面体分解相反，显式转换而非依赖自动调序。
    if volume.cellTypes(k)==13,cellNodes=cellNodes([2 1 3 5 4 6]);end
    fprintf(fid,'%d',volume.cellTypes(k));fprintf(fid,' %d',cellNodes-1);fprintf(fid,'\n');
end
fprintf(fid,'NPOIN= %d\n',size(volume.vertices,1));
fprintf(fid,'%.17g %.17g %.17g %d\n',[volume.vertices,(0:size(volume.vertices,1)-1)']');
fprintf(fid,'NMARK= 2\n');markers={volume.wall,volume.farfield};names={'wall','farfield'};
for j=1:2
    faces=markers{j};fprintf(fid,'MARKER_TAG= %s\nMARKER_ELEMS= %d\n',names{j},numel(faces));
    for k=1:numel(faces)
        % SU2原生边界法向指向流体域内部；几何体积审计中的外向面保持不变。
        f=faces{k};f=f([1 end:-1:2]);
        if numel(f)==3,type=5;elseif numel(f)==4,type=9;else,error('无效边界面。');end
        fprintf(fid,'%d',type);fprintf(fid,' %d',f-1);fprintf(fid,'\n');
    end
end
clear guard;
R=287.05287;gamma=1.4;muRef=1.458e-6*273.15^1.5/(273.15+110.4);
refArea=pi*max(hypot(volume.sourceVertices(:,2),volume.sourceVertices(:,3)))^2;
cfg={'SOLVER= NAVIER_STOKES';'KIND_TURB_MODEL= NONE';'MATH_PROBLEM= DIRECT'; ...
    'RESTART_SOL= NO';'MESH_FILENAME= mesh.su2';'MESH_FORMAT= SU2'; ...
    'SYSTEM_MEASUREMENTS= SI';'REF_DIMENSIONALIZATION= DIMENSIONAL';'FLUID_MODEL= IDEAL_GAS'; ...
    sprintf('GAMMA_VALUE= %.17g',gamma);sprintf('GAS_CONSTANT= %.17g',R); ...
    'INIT_OPTION= TD_CONDITIONS';sprintf('MACH_NUMBER= %.17g',flight.Mach); ...
    sprintf('AOA= %.17g',180-flight.alpha);sprintf('SIDESLIP_ANGLE= %.17g',-flight.beta); ...
    sprintf('FREESTREAM_TEMPERATURE= %.17g',flight.temperature);sprintf('FREESTREAM_PRESSURE= %.17g',flight.pressure); ...
    'VISCOSITY_MODEL= SUTHERLAND';sprintf('MU_REF= %.17g',muRef);'MU_T_REF= 273.15';'SUTHERLAND_CONSTANT= 110.4'; ...
    'CONDUCTIVITY_MODEL= CONSTANT_PRANDTL';'PRANDTL_LAM= 0.72'; ...
    thermalLine;'MARKER_FAR= (farfield)'; ...
    'MARKER_PLOTTING= (wall)';'MARKER_MONITORING= (wall)';sprintf('REF_AREA= %.17g',refArea); ...
    'NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES';'CONV_NUM_METHOD_FLOW= AUSM'; ...
    'MUSCL_FLOW= YES';'SLOPE_LIMITER_FLOW= VENKATAKRISHNAN_WANG';'VENKAT_LIMITER_COEFF= 0.05'; ...
    'TIME_DISCRE_FLOW= EULER_IMPLICIT';'LINEAR_SOLVER= BCGSTAB';'LINEAR_SOLVER_PREC= LU_SGS'; ...
    'LINEAR_SOLVER_ERROR= 1E-5';'LINEAR_SOLVER_ITER= 10';'MGLEVEL= 0'; ...
    'CFL_NUMBER= 0.1';'CFL_ADAPT= YES';'CFL_ADAPT_PARAM= (0.5, 1.1, 0.05, 5.0)'; ...
    sprintf('ITER= %d',flight.iterations);'CONV_FIELD= (RMS_DENSITY, RMS_ENERGY)'; ...
    'CONV_RESIDUAL_MINVAL= -8';'CONV_STARTITER= 10';'TABULAR_FORMAT= CSV'; ...
    'HISTORY_OUTPUT= (ITER, RMS_RES, AERO_COEFF, HEAT)';'CONV_FILENAME= history'; ...
    'SURFACE_FILENAME= surface_flow';'RESTART_FILENAME= restart_flow';'VOLUME_FILENAME= volume_flow'; ...
    'OUTPUT_FILES= (SURFACE_PARAVIEW_ASCII, RESTART_ASCII)'; ...
    sprintf('OUTPUT_WRT_FREQ= %d',min(flight.iterations,flight.outputEvery));['WRT_RESTART_OVERWRITE= ' overwrite];['WRT_SURFACE_OVERWRITE= ' overwrite]; ...
    'SCREEN_WRT_FREQ_INNER= 50';'VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)'};
if spatial
    cfg{end+1}='MARKER_PYTHON_CUSTOM= (wall)';
    % 统一标记的默认温度只供上游初始化；实际每个节点必须被表值覆盖并取得回执。
    fid=fopen(fullfile(folder,'wall_temperature.csv'),'w','n','UTF-8');assert(fid>=0);
    guard=onCleanup(@()fclose(fid));fprintf(fid,'RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K\n');
    fprintf(fid,'%.17g,%.17g\n',[spec.x(:),spec.temperature(:)]');clear guard;
end
fid=fopen(fullfile(folder,'case.cfg'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',cfg{:});clear guard;
direction=[-cosd(flight.alpha)*cosd(flight.beta),-sind(flight.beta),sind(flight.alpha)*cosd(flight.beta)];
info=struct('kind',volume.kind,'flight',flight,'thermalBoundary',thermal,'expectedVelocity',flight.Mach*sqrt(gamma*R*flight.temperature)*direction, ...
    'nodes',size(volume.vertices,1),'cells',numel(volume.cells),'hexes',sum(volume.cellTypes==12), ...
    'wedges',sum(volume.cellTypes==13),'tetrahedra',sum(volume.cellTypes==10),'wallFaces',numel(volume.wall),'wallArea',volume.bodyArea, ...
    'bodyVolume',volume.bodyVolume,'fluidVolume',volume.fluidVolume,'minimumSubTetVolume',volume.minimumSubTetVolume, ...
    'regions',{unique(volume.wallRegion)},'referenceArea',refArea,'originalSurfaceFaceCount',size(volume.sourceFaces,1), ...
    'coordinateConvention','X-axis; nose +X; zero-incidence flow -X; positive alpha flow +Z; SI', ...
    'thermalPropertyMode','caloricallyPerfect; gamma1.4; Sutherland; Pr0.72', ...
    'nativeOrientation','SU2v8.6.0 prism inward triangular ends; boundary normals into fluid domain', ...
    'physicalSolution',false,'productionQualified',false,'fieldProjectionImplemented',true);
if spatial
    info.requiredRuntimeFlavor=thermal.requiredRuntimeFlavor;
    info.thermalProfileSHA256=sha(fullfile(folder,'wall_temperature.csv'));
end
fid=fopen(fullfile(folder,'case_manifest.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(info));
end
function h=sha(file)
% 文件身份与实际温度表一起传给监督器，避免边界被改写后继续使用旧清单。
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
md.update(fread(fid,inf,'*uint8'));h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
