function [mesh,stations,s]=loadShortSupportSurface(g,job,resolution)
% 模块一：核对同一几何、温界、收敛端点和五区交叠文件，从存档或新算目录取数。
if nargin<3,resolution='coarse';end
root=fileparts(mfilename('fullpath'));cases=fullfile(root,'cases','short_sphere');
transferName='supportedTransfer65.npz';
if strcmp(resolution,'angular'),transferName='supportedAngularTransfer65.npz';end
if nargin<2 || isempty(job)
    if strcmp(resolution,'angular')
        pair=readShortSupportAngularPair(cases);archive=fullfile(cases,'supportedAngularSurface65.zip');
        assert(strcmp(cm011FileSHA256(archive),pair.fineArchiveSHA256));
    elseif strcmp(resolution,'fine')
        pair=readShortSupportGridPair(cases);archive=fullfile(cases,'supportedFine65.zip');
        assert(strcmp(cm011FileSHA256(archive),pair.fineArchiveSHA256));
    else
        archive=fullfile(cases,'supportedSurface65.zip');
        assert(strcmp(cm011FileSHA256(archive),'2e5b065e3b7bd7f15dc46a22bdc6f5a041bf57ebb9a158c31cb46426953c772b'));
    end
else
    ack=jsondecode(fileread(fullfile(job,'completed.json')));archive=fullfile(job,'surface','surface.zip');
    assert(strcmp(cm011FileSHA256(archive),ack.archiveSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(ack.source,'thermal.json')),ack.thermalSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(ack.source,'early','initial_flow.csv')),ack.frameSHA256));
end
temp=tempname;mkdir(temp);clean=onCleanup(@()discard(temp)); %#ok<NASGU>
unzip(archive,temp);s=jsondecode(fileread(fullfile(temp,'summary.json')));
assert(strcmp(cm011FileSHA256(fullfile(cases,'test65.json')),s.caseSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(cases,'supported65.zip')),s.bodySHA256) && ...
    strcmp(cm011FileSHA256(fullfile(cases,transferName)),s.transferSHA256));
assert(s.thermal.passed && abs(s.thermal.early.physicalGapFraction)<.005 && ...
    abs(s.thermal.late.physicalGapFraction)<.005 && s.thermal.thermalChange<.002 && ...
    strcmp(s.thermal.early.sourceFrame,s.frameSHA256));
for key={'wall','stations'}
    name=[key{1} '.csv'];assert(strcmp(cm011FileSHA256(fullfile(temp,name)),s.files.(matlab.lang.makeValidName(name))));
end
w=readtable(fullfile(temp,'wall.csv'));stations=readtable(fullfile(temp,'stations.csv'));
unzip(fullfile(cases,'supported65.zip'),temp);
geom=jsondecode(fileread(fullfile(temp,'geometry.json')));
assert(isequal(geom.F,g.F) && isequal(geom.regionID,g.regionID) && max(abs(geom.V-g.V),[],'all')<1e-14, ...
    'RAH:Short:FaceOrder','当前几何原面顺序与冻结解不一致。');
n=size(g.F,1);c=s.condition;
assert(isequal(w.faceID,(1:n)') && isequal(w.regionID,g.regionID) && max(abs(w.areaM2-g.areaM2))<1e-17);
assert(max(abs(w.absolutePa-w.gaugePa-c.airPa))<1e-8 && abs(w.areaM2'*w.fluxInWm2+s.thermal.early.wallW)<1e-10);
assert(isequal(stations.originalFaceID,g.stations.faceID) && s.jointExperimentalAreaM2==0);
% 模块二：用统一的原面适配器交付压力、热流、给定壁温及对应首层流体温度。
names=["外露银球";"+Y 支管外壁";"−Y 支管外壁";"+Y 湿端盖";"−Y 湿端盖"];
condition=struct('caseLabel','Short Test 65 · 双支管银球 · 稳态','velocityMps',[-c.speedMps 0 0], ...
    'ambientPressurePa',c.airPa,'ambientTemperatureK',c.airK,'reynolds',c.re,'alphaDeg',0,'betaDeg',0);
setup=struct('condition',condition,'region',names(g.regionID),'wallTemperatureK',w.wallInputK, ...
    'engineeringSolver',@empty,'localSolver',@wall);
[mesh,~,coverage]=runHybridSurfaceCase(g,struct('hybridSurface',setup));
assert(coverage.fullSurfaceComputed && abs(coverage.wholeSurfaceHeatPowerW+s.thermal.early.wallW)<1e-10);
mesh.hybridWall.pressureGaugePa=w.gaugePa;mesh.hybridWall.pressureGaugeMethod=mesh.hybridWall.pressureMethod;
mesh.shortSupportSurface=true;
id=mesh.evidenceIdentity;id.meshSHA256=s.thermal.meshSHA256;
id.sourceFrame=struct('kind','SU2 8.6.0 · 双支管银球有界稳态 · 统一远场质量通量','sha256',s.frameSHA256, ...
    'caseSHA256',s.caseSHA256,'reactionSHA256',s.releaseSHA256,'wallSHA256',s.files.wall_csv);
id.sourceTraceable=true;id.productionQualified=false;
id.scope='名义一球径双支管、湿端盖及给定壁温。稳态层流研究解；空间误差和逐面实验精度待验证。';
mesh.evidenceIdentity=id;s.coverage=coverage;
    function t=empty(~,~,~)
        t=table(zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1),'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
    end
    function out=wall(~,~,~,~,ids,~)
        j=ids+1;count=numel(j);
        t=table(ids,w.absolutePa(j),w.fluxInWm2(j),repmat("SU2 8.6.0 三维不可压层流",count,1), ...
            'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
        t.pressureMethod=repmat("SU2 壁节点表压按对偶面积和五区交叠加权，绝对压加环境压力",count,1);
        t.heatFluxMethod=repmat("同态释放壁反力与节点余量通量修正，按五区交叠守恒分配，正值入壁",count,1);
        t.nearWallCellTemperatureK=w.firstK(j);
        t.nearWallTemperatureMethod=repmat("对应壁节点的首层温差标量，由对偶面积和分区交叠换算到 K",count,1);
        out=struct('accepted',true,'reason','有界稳态、同态热量守恒及热稳定均通过','wall',t);
    end
end
function discard(folder)
% 临时目录仅展开已核实的几何包及结果包，逐项释放。
for name={'summary.json','geometry.json','parents.csv','stations.csv','supported_m.stl','wall.csv'}
    p=fullfile(folder,name{1});if isfile(p),delete(p);end
end
rmdir(folder);
end
