function [result,stations,report]=runShortSphereWorkbench(mode,pulse)
% 模块一：从冻结工况冷态计算或读取已登记基流，核对源帧及同态回执。
if nargin<1,mode='archive';end
if nargin<2,pulse=[];end
mode=validatestring(mode,{'archive','solve'});root=fileparts(mfilename('fullpath'));
cases=fullfile(root,'cases','short_sphere');caseFile=fullfile(cases,'test65.json');source=jsondecode(fileread(caseFile));
key=jsondecode(fileread(fullfile(cases,'surface65.json')));job='';
assert(strcmp(cm011FileSHA256(caseFile),key.caseSHA256),'RAH:Short:Source','银球工况指纹不符。');
if strcmp(mode,'solve')
    job=fullfile(root,'outputs',['short65_' char(java.util.UUID.randomUUID())]);
    runCM011Python({fullfile(root,'scripts','shortSphere65.py'),job},[job '.log'],pulse,true);
    ack=jsondecode(fileread(fullfile(job,'completed.json')));folder=fullfile(job,'surface');expected=ack.summary;
    assert(strcmp(ack.case,key.caseSHA256) && strcmp(cm011FileSHA256(fullfile(job,'endpoints.json')),ack.endpoints));
    endpoint=jsondecode(fileread(fullfile(job,'endpoints.json')));
else
    archive=fullfile(cases,'surface65.zip');assert(strcmp(cm011FileSHA256(archive),key.archiveSHA256));
    folder=tempname;mkdir(folder);cleanup=onCleanup(@()discard(folder)); %#ok<NASGU>
    unzip(archive,folder);expected=key.summarySHA256;
end
path=fullfile(folder,'summary.json');assert(strcmp(cm011FileSHA256(path),expected));
s=jsondecode(fileread(path));assert(strcmp(s.caseSHA256,key.caseSHA256));
if strcmp(mode,'solve')
    stage=fullfile(job,endpoint.stage);
    assert(strcmp(endpoint.frame,s.frameSHA256) && strcmp(cm011FileSHA256(fullfile(stage,'restart.csv')),s.frameSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(stage,'reaction','restart.csv')),s.reactionSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(job,'prepared','mesh.su2')),s.meshSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(folder,'surface.zip')),ack.archive),'RAH:Short:State','终帧、反力或归档与本次任务不符。');
end
for name={'wall.csv','body.json','stations.csv'}
    file=name{1};assert(strcmp(cm011FileSHA256(fullfile(folder,file)),s.files.(matlab.lang.makeValidName(file))));
end
g=s.balance;assert(abs(g.balance)<.005 && g.extra20Change<.002 && all(structfun(@(v)v<-8,g.window)) && all(structfun(@(v)v<-8,g.extraWindow)));
% 模块二：以本试件自己的每一个三角面交付数值，保留给定壁温、首层温度及热流方向。
body=jsondecode(fileread(fullfile(folder,'body.json')));w=readtable(fullfile(folder,'wall.csv'));
stations=readtable(fullfile(folder,'stations.csv'));n=size(body.F,1);c=s.condition;
assert(n==s.faces && isequal(w.faceID,(1:n)') && isequal(stations.angleDeg,source.angleDeg) && isequal(stations.reportedNu,source.localNuSurface));
assert(max(abs(stations.deltaNu-(stations.computedNu-stations.reportedNu)))<1e-12);
f=body.F;v=body.V;area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(max(abs(w.areaM2-area))<1e-18 && abs(area'*w.fluxInWm2+s.wholeOutW)<1e-11 && max(abs(w.absolutePa-w.gaugePa-c.airPa))<1e-8);
region=repmat("迎风半球",n,1);region(w.angleDeg>90)="背风半球";
condition=struct('caseLabel','Short Test 65 · 光滑银球稳态基流','velocityMps',[-c.speedMps 0 0], ...
    'ambientPressurePa',c.airPa,'ambientTemperatureK',c.airK,'reynolds',c.re,'alphaDeg',0,'betaDeg',0);
setup=struct('condition',condition,'region',region,'wallTemperatureK',w.wallInputK,'engineeringSolver',@empty,'localSolver',@wall);
[result,~,coverage]=runHybridSurfaceCase(body,struct('hybridSurface',setup));
assert(coverage.fullSurfaceComputed && abs(coverage.wholeSurfaceHeatPowerW+s.wholeOutW)<1e-11);
result.hybridWall.pressureGaugePa=w.gaugePa;
result.hybridWall.pressureGaugeMethod=repmat("SU2 节点相对压按壁面线性权重求平均",n,1);
result.shortSphere=struct('angleDeg',w.angleDeg,'doi',source.doi,'test',source.test);
id=result.evidenceIdentity;id.meshSHA256=s.meshSHA256;
id.sourceFrame=struct('kind','SU2 8.6.0 全周稳态层流基流','sha256',s.frameSHA256,'reactionSHA256',s.reactionSHA256, ...
    'caseSHA256',s.caseSHA256,'wallSHA256',s.files.wall_csv);
id.sourceTraceable=true;id.productionQualified=false;id.scope='光滑银球稳态基流；支管、来流湍动和非定常尾迹未计入。局部 Nu 按试验定义对照。';
result.evidenceIdentity=id;report=s;report.coverage=coverage;report.action=mode;report.jobDirectory=job;
    function t=empty(~,~,~)
        t=table(zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1),'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
    end
    function out=wall(~,~,~,~,ids,~)
        j=ids+1;count=numel(j);
        t=table(ids,w.absolutePa(j),w.fluxInWm2(j),repmat("SU2 8.6.0 三维不可压层流及强壁温被动标量",count,1), ...
            'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
        t.pressureMethod=repmat("节点表压换算 Pa 后按三角面平均并加环境压力",count,1);
        t.heatFluxMethod=repmat("固定终帧释放强标量壁行，节点反力按面积守恒积分，正值入壁",count,1);
        t.nearWallCellTemperatureK=w.firstK(j);t.nearWallTemperatureMethod=repmat("同位壁节点对应的实际首层标量，按面积权重及 90 K 温标换算",count,1);
        out=struct('accepted',true,'reason','稳态末窗、追加步热稳定及同态热量账通过','wall',t);
    end
end
function discard(folder)
% 模块三：只释放本次展开的四个已登记小文件。
for name={'wall.csv','stations.csv','body.json','summary.json'}
    path=fullfile(folder,name{1});if isfile(path),delete(path);end
end
rmdir(folder);
end
