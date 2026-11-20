function [current,report,page]=runThermistorBeadCase(folder,mode,windowOpen)
%RUNTHERMISTORBEADCASE 读取 B 型热敏珠的冻结结果，或从冷初场计算同一工况。
% [mesh,report]=runThermistorBeadCase() 只读随项目发布的原面表。
% [mesh,report,page]=runThermistorBeadCase([], 'solve') 运行内嵌求解器并导出离线页。
% 显式指定 outputs 下的新目录，或请求第三个输出，都会在那里导出离线页。
root=fileparts(mfilename('fullpath'));outputs=fullfile(root,'outputs');
if nargin<1,folder='';end
if nargin<2,mode='archive';end
if nargin<3,windowOpen=[];end
mode=validatestring(mode,{'archive','solve'});fresh=strcmp(mode,'solve');
writePage=fresh || ~isempty(folder) || nargout>=3;page='';
if writePage
    if isempty(folder),folder=fullfile(outputs,['thermistor_' char(java.util.UUID.randomUUID())]);end
    targetDir=java.io.File(char(folder));
    if ~targetDir.isAbsolute(),targetDir=java.io.File(pwd,char(folder));end
    folder=char(targetDir.getCanonicalPath());
    assert(startsWith(folder,[char(java.io.File(outputs).getCanonicalPath()) filesep],'IgnoreCase',ispc) && ...
        java.io.File(folder).mkdirs(),'RAH:Bead:Output','请选择项目 outputs 下尚未使用的目录。');
end

% 模块一：新计算经过同网格三维动量、强温差和反力核对；浏览时校验归档后再解包。
caseDir=fullfile(root,'cases','thermistor_b');
if fresh
    native=fullfile(folder,'native');
    pulse=[];
    if ~isempty(windowOpen),pulse=@checkWindow;end
    runCM011Python({fullfile(root,'scripts','runThermistorCase.py'),native},fullfile(folder,'solve.log'),pulse);
    geo=fullfile(native,'strong');result=fullfile(native,'result');
else
    zip=fullfile(caseDir,'surface.zip');
    assert(strcmp(cm011FileSHA256(zip),'11ba9c504f8df41a1be7c317d1047e77ab76800d3eae7b569d0b34333ffbb563'), ...
        'RAH:Bead:Archive','B 型珠体的归档指纹不符。');
    geo=tempname;mkdir(geo);clean=onCleanup(@()removeExtracted(geo)); %#ok<NASGU>
    unzip(zip,geo);result=geo;
end
spec=jsondecode(fileread(fullfile(geo,'case.json')));
stats=jsondecode(fileread(fullfile(result,'summary.json')));
rowsFile=fullfile(result,'wall.csv');bodyFile=fullfile(geo,'body.json');
assert(strcmp(spec.grid,'outer') && strcmp(stats.grid,'outer') && ...
    strcmp(cm011FileSHA256(fullfile(geo,'case.json')),'46d981fce8d8b3f456f0a36cfde7ccd84cd47cb721d2d71a31ceea87fded038c') && ...
    strcmp(cm011FileSHA256(bodyFile),stats.bodySHA256) && strcmp(cm011FileSHA256(rowsFile),stats.wallTableSHA256), ...
    'RAH:Bead:Binding','工况、原面表或三角面与数值回执不符。');
body=jsondecode(fileread(bodyFile));mesh=struct('V',body.V,'F',body.F);
rows=readtable(rowsFile);n=size(mesh.F,1);
assert(stats.fullSurfaceComputed && stats.originalFaces==n && n==5216 && ...
    isequal(rows.faceID,(0:n-1)') && stats.boundaryRelativeDifference<.005 && stats.heatChange20<.002 && ...
    all(structfun(@(v)v < -8,stats.last20MaxLogRMS)) && all(structfun(@(v)v < -8,stats.continued20MaxLogRMS)), ...
    'RAH:Bead:Receipt','原面编号、热量账或收敛窗未通过。');
assert(all(rows.pressurePa>0 & rows.heatFluxWm2<0) && all(rows.nearWallCellTemperatureK>spec.ambientK & ...
    rows.nearWallCellTemperatureK<spec.wallK) && max(abs(rows.pressurePa-spec.pressurePa-rows.pressureGaugePa))<2e-10, ...
    'RAH:Bead:Fields','压力、热流或首层温度读数不符。');
f=mesh.F;v=mesh.V;area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(all(abs(rows.areaM2-area)<1e-10*area),'RAH:Bead:Area','逐面面积与原几何不符。');

% 模块二：同一名义等温珠体的全部原面交给区域求解入口，保留每个变量的算法来源。
condition=struct('caseLabel','B 型热敏珠 · 名义长椭球 · 沿长轴来流','velocityMps',[-spec.freestreamMps 0 0], ...
    'ambientPressurePa',spec.pressurePa,'ambientTemperatureK',spec.ambientK,'wallTemperatureK',spec.wallK, ...
    'alphaDeg',0,'betaDeg',0,'diameterM',spec.diameterM,'geometryScope','光滑等温珠体，未建引线');
setup=struct('condition',condition,'region',repmat("热敏珠完整外表面",n,1), ...
    'wallTemperatureK',repmat(spec.wallK,n,1),'engineeringSolver',@emptyWall,'localSolver',@localWall);
[current,~,report]=runAeroHeatSingleCase(mesh,struct('analysisMode','hybridSurface','hybridSurface',setup));
assert(report.fullSurfaceComputed && isequal(current.hybridWall.faceID,rows.faceID) && ...
    abs(report.wholeSurfaceHeatPowerW/stats.powerW-1)<1e-11,'RAH:Bead:Power','完整原面功率与回执不符。');
current.hybridWall.pressureGaugePa=rows.pressureGaugePa;
identity=struct('kind','SU2-8.6.0-thermistor-B','sha256',stats.sourceFrameSHA256, ...
    'reactionSHA256',stats.sourceReactionSHA256,'wallSHA256',stats.wallTableSHA256,'bodySHA256',stats.bodySHA256);
current.evidenceIdentity.sourceFrame=identity;current.evidenceIdentity.sourceTraceable=true;
current.evidenceIdentity.productionQualified=false;
current.evidenceIdentity.scope='名义等温珠体的数值解。实验只给整珠平均换热系数。';

% 模块三：按论文的名义面积比较总功率；逐面验证面积、给定温度与首层温度独立保存。
comparison=jsondecode(fileread(fullfile(root,'docs','validation','v282','evidence','comparison.json')));
meanH=-report.wholeSurfaceHeatPowerW/((spec.wallK-spec.ambientK)*spec.areaSourceM2);
assert(abs(meanH/stats.nominalExperimentAreaCoefficientWm2K-1)<1e-11,'RAH:Bead:Mean','名义面积换热积分不符。');
current.thermistorBead=struct('experiment',comparison.experiment,'meanH',meanH, ...
    'relativeDifference',meanH/comparison.experiment.hWm2K-1,'verifiedJointAreaM2',0);
report.native=stats;report.reusedArchive=~fresh;report.comparison=current.thermistorBead;
if writePage,page=exportAerothermalWeb(current,fullfile(folder,'workbench'));end

    function checkWindow()
        % 关闭主窗口后通知当前原生监督器收回本次计算，其他工况继续独立运行。
        drawnow limitrate;
        if ~windowOpen()
            fid=fopen(fullfile(folder,'thermistor_cancel.request'),'w');
            if fid>=0,fclose(fid);end
            error('RAH:Bead:Closed','热敏珠窗口已关闭。');
        end
    end
    function empty=emptyWall(~,~,~)
        % 对这个热敏珠，采用同一连通外域的数值解覆盖全部表面。
        empty=rows([],{'faceID','pressurePa','heatFluxWm2'});empty.method=strings(0,1);
    end
    function solution=localWall(~,~,~,~,target,~)
        wall=rows(target+1,{'faceID','pressurePa','heatFluxWm2','nearWallCellTemperatureK'});
        count=height(wall);
        wall.method=repmat("SU2 8.6.0 三维层流，30° 对称扇区",count,1);
        wall.pressureMethod=repmat("原生壁节点相对压力按面内权重平均，加环境静压",count,1);
        wall.heatFluxMethod=repmat("强等温对流扩散方程的壁反力，按原三角面积守恒分配",count,1);
        wall.nearWallTemperatureMethod=repmat("首层流体节点的温差标量换算为 K，再按面内权重平均",count,1);
        solution=struct('accepted',true,'reason','六个分量的末窗残差、温界与边界热量账通过', ...
            'wall',wall,'diagnostics',struct('sourceFrameSHA256',stats.sourceFrameSHA256));
    end
end

function removeExtracted(folder)
% 模块四：释放此调用独有的四个临时文件。新计算的现场和用户导出保存在项目中。
for name={'case.json','body.json','wall.csv','summary.json'}
    path=fullfile(folder,name{1});if isfile(path),delete(path);end
end
rmdir(folder);
end
