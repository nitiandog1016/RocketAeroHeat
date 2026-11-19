function [current,report,page]=runCM011TunnelCase(folder,grid,savedRun,energyObservation,sectorRun,sectorObservation)
%RUNCM011TUNNELCASE 在普通单工况入口内计算低速等温 CM-011 缩比试件，导出离线三维结果。
% 字段和方法逐面返回。第三个参数可显式回读本项目同一工况的已核对原生终态。
% 第四个参数可附同一终场的边界能量观测；来源和收支核对后在热流视图显示。
% 第五个参数可附低速同工况的三维扇区目录，在同一组原面上比较两种离散。
% 第六个参数是该三维保存态的可选热收支观测，单独绑定到三维结果。
root=fileparts(mfilename('fullpath'));outputs=fullfile(root,'outputs');
if nargin<1 || isempty(folder)
    folder=fullfile(outputs,['cm011_' char(java.util.UUID.randomUUID())]);
end
if nargin<2,grid='coarse';end
if nargin<3,savedRun='';end
if nargin<4,energyObservation='';end
if nargin<5,sectorRun='';end
if nargin<6,sectorObservation='';end
assert(isempty(sectorObservation) || ~isempty(sectorRun),'RAH:CM011:SectorObservation','请指定三维观测对应的运行目录。');
assert(isempty(energyObservation) || ~isempty(savedRun),'RAH:CM011:Observation', ...
    '读入能量观测时请同时指定该观测对应的原生终场目录。');
grid=validatestring(grid,{'coarse','fine','outer'});
assert(isempty(sectorRun) || (~isempty(savedRun) && strcmp(grid,'fine')), ...
    'RAH:CM011:Sector','三维对照需要绑定的 12D 细网终场。');
folder=canonicalPath(folder);
assert(startsWith(folder,[char(java.io.File(outputs).getCanonicalPath()) filesep],'IgnoreCase',ispc), ...
    'RAH:CM011:Output','请选择项目 outputs 下的新结果目录。');
assert(java.io.File(folder).mkdirs(),'RAH:CM011:Output','须使用新的结果目录。');
prepare=fullfile(root,'scripts','cm011Tunnel.py');reader=fullfile(root,'scripts','readCM011Tunnel.py');
solver=fullfile(root,'scripts','runEmbeddedSU2.py');
fromArchive=~isempty(savedRun);
if fromArchive
    native=canonicalPath(savedRun);
else
    native=fullfile(folder,'flow');
    pythonCommand({prepare,native,grid,'--prepare-only'},'prepare');
end
body=jsondecode(fileread(fullfile(native,'body.json')));
assert(strcmp(body.grid,grid),'RAH:CM011:Grid','指定网格与保存试件不符。');
mesh=struct('V',body.V,'F',body.F);
condition=struct('caseLabel','CM-011 低速等温缩比试件 · 轴对称光滑主体', ...
    'velocityMps',[-0.1 0 0],'ambientPressurePa',101325,'ambientTemperatureK',300, ...
    'wallTemperatureK',300.1,'diameterM',body.diameterM,'gravityModel','忽略浮力', ...
    'geometryScope','按盾面外径缩尺 1:1000；无实物偏置、倾角、附件和支撑。');
labels=string(body.regionNames);region=labels(body.regionCode+1);
setup=struct('condition',condition,'region',region(:),'wallTemperatureK',repmat(condition.wallTemperatureK,size(mesh.F,1),1), ...
    'engineeringSolver',@noWallEstimate,'localSolver',@localWall);
rows=[];stats=[];readFolder=fullfile(folder,'audit');
[current,~,report]=runAeroHeatSingleCase(mesh,struct('analysisMode','hybridSurface','hybridSurface',setup));
assert(report.fullSurfaceComputed && isequal(current.hybridWall.faceID,rows.faceID) && ...
    abs(report.wholeSurfaceHeatPowerW/stats.triangularPowerW-1)<1e-10,'RAH:CM011:Mapping','原面编号或功率不一致。');
% 壁节点流体温度、第一层流体温度和给定壁温分别来自不同的数据列。
current.hybridWall.wallFluidTemperatureK=rows.wallFluidTemperatureK;
current.hybridWall.wallFluidTemperatureMethod=repmat("SU2 能量方程的壁节点流体温度",height(rows),1);
current.hybridWall.pressureGaugePa=current.hybridWall.pressurePa-condition.ambientPressurePa;
report.native=stats;report.reusedNativeRun=fromArchive;
if isfield(stats,'nativeEnergyBalance'),current.hybridEnergyBalance=stats.nativeEnergyBalance;end
% 只有这一可选分支会回读三维帧。主结果保留二维终场的身份和原面分母。
if ~isempty(sectorRun)
    sectorFolder=fullfile(folder,'sector_audit');
    sectorArgs={fullfile(root,'scripts','readCM011TunnelSector.py'), ...
        canonicalPath(sectorRun),native,sectorFolder};
    if ~isempty(sectorObservation)
        sectorArgs=[sectorArgs,{'--energy-observation',canonicalPath(sectorObservation)}];
    end
    pythonCommand(sectorArgs,'sector_audit');
    sector=jsondecode(fileread(fullfile(sectorFolder,'summary.json')));
    sectorRows=readtable(fullfile(sectorFolder,'wall.csv'));
    assert(sector.originalFaces==height(rows) && isequal(rows.faceID,sectorRows.faceID), ...
        'RAH:CM011:SectorMapping','三维对照与主场原面编号不一致。');
    current.cm011Sector=struct('wall',sectorRows,'summary',sector);
    report.sector=sector;
end
stlwrite(triangulation(mesh.F,mesh.V),fullfile(folder,'CM011-mm-scale-SI.stl'));
page=exportAerothermalWeb(current,fullfile(folder,'workbench'));
save(fullfile(folder,'result.mat'),'current','report','condition');
fid=fopen(fullfile(folder,'summary.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
summary=struct('condition',condition,'coverage',report.coverageMessage,'reusedNativeRun',fromArchive,'native',stats);
if isfield(report,'sector'),summary.sector=report.sector;end
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));
fprintf('%s\n%s\n',report.coverageMessage,page);

    function empty=noWallEstimate(~,~,~)
        % 当前低速区统一由黏性求解器提供壁面场。
        empty=table(zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
            'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
    end
    function result=localWall(g,c,tw,~,target,~)
        assert(isequal(g.V,mesh.V) && isequal(g.F,mesh.F) && isequal(c,condition) && all(tw==c.wallTemperatureK));
        % 四个壁区共享同一个连通外域，每个原面只返回一次。
        if isempty(rows)
            if ~fromArchive,pythonCommand({solver,'--request',fullfile(native,'engine_request.json')},'solve');end
            auditArgs={};
            if ~isempty(energyObservation)
                auditArgs={'--native-balance',canonicalPath(energyObservation)};
            end
            pythonCommand([{reader,native,readFolder},auditArgs],'audit');
            stats=jsondecode(fileread(fullfile(readFolder,'summary.json')));
            if ~fromArchive && ~stats.acceptedNumericalWall && max(stats.last20MaximumLogRMS)<-4
                next=fullfile(folder,'continued');
                pythonCommand({prepare,next,grid,'--restart',native},'continue');
                native=next;readFolder=fullfile(folder,'continued_audit');
                pythonCommand({reader,native,readFolder},'continued_audit');
                stats=jsondecode(fileread(fullfile(readFolder,'summary.json')));
            end
            assert(stats.acceptedNumericalWall,'RAH:CM011:State','当前终场未达到数值门槛，详见 %s。',readFolder);
            rows=readtable(fullfile(readFolder,'wall.csv'));
            assert(isequal(rows.faceID,(0:size(g.F,1)-1)'),'RAH:CM011:IDs','热壁原面编号不完整。');
        end
        wall=rows(target+1,{'faceID','pressurePa','heatFluxWm2','nearWallCellTemperatureK'});
        wall.method=repmat("轴对称不可压 SU2 层流与能量方程",height(wall),1);
        wall.pressureMethod=repmat("SU2 数值相对压力，加试验环境参考压力",height(wall),1);
        wall.heatFluxMethod=repmat("SU2 原生等温壁 Fourier 热流，按母线节点映射",height(wall),1);
        wall.nearWallTemperatureMethod=repmat("SU2 顶点中心有限体积，首层流体节点温度",height(wall),1);
        result=struct('accepted',true,'reason',"四分量残差、壁面热流窗口与原生终场回读通过；实验精度待验", ...
            'wall',wall,'diagnostics',struct('sourceVTKSHA256',{stats.sourceVTKSHA256}));
    end
    function pythonCommand(args,name)
        % 独立参数启动项目 Python，路径空格或括号保持原样。失败日志留在本次目录。
        items=[{'python'},args];argv=javaArray('java.lang.String',numel(items));
        for k=1:numel(items),argv(k)=java.lang.String(items{k});end
        cmd=java.lang.ProcessBuilder(argv);cmd.directory(java.io.File(root));cmd.redirectErrorStream(true);
        logPath=fullfile(folder,[name '.log']);cmd.redirectOutput(java.io.File(logPath));
        proc=cmd.start();code=proc.waitFor();
        assert(code==0,'RAH:CM011:Process','%s 阶段退出 %d，见 %s。',name,code,logPath);
    end
end

function result=canonicalPath(value)
% 模块用途：相对路径随 MATLAB 当前目录解析，包括批处理启动后曾经 cd 的情况。
file=java.io.File(char(value));
if ~file.isAbsolute(),file=java.io.File(pwd,char(value));end
result=char(file.getCanonicalPath());
end
