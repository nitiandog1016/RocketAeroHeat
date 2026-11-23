function [result,stations,report]=runShortSphereWorkbench(mode,pulse)
% 模块一：从冻结工况冷态计算或读取已登记基流，核对源帧及同态回执。
if nargin<1,mode='archive';end
if nargin<2,pulse=[];end
mode=validatestring(mode,{'archive','solve','wakeFine','wakeCoarse','wakeStart','wakeSolve'});root=fileparts(mfilename('fullpath'));
cases=fullfile(root,'cases','short_sphere');caseFile=fullfile(cases,'test65.json');source=jsondecode(fileread(caseFile));
wake=startsWith(mode,'wake');branch='fine';if strcmp(mode,'wakeCoarse'),branch='coarse';end
stem='surface65';if any(strcmp(mode,{'wakeFine','wakeCoarse'})),stem=['wake_' branch];end
key=jsondecode(fileread(fullfile(cases,[stem '.json'])));job='';
assert(strcmp(cm011FileSHA256(caseFile),key.caseSHA256),'RAH:Short:Source','银球工况指纹不符。');
if strcmp(mode,'solve')
    job=fullfile(root,'outputs',['short65_' char(java.util.UUID.randomUUID())]);
    runCM011Python({fullfile(root,'scripts','shortSphere65.py'),job},[job '.log'],pulse,true);
    ack=jsondecode(fileread(fullfile(job,'completed.json')));folder=fullfile(job,'surface');expected=ack.summary;
    assert(strcmp(ack.case,key.caseSHA256) && strcmp(cm011FileSHA256(fullfile(job,'endpoints.json')),ack.endpoints));
    endpoint=jsondecode(fileread(fullfile(job,'endpoints.json')));
elseif any(strcmp(mode,{'wakeStart','wakeSolve'}))
    % 从真实基流出发，分别核对粗步长和细步长在同一时刻的终点。
    horizon=8;if strcmp(mode,'wakeStart'),horizon=.25;end
    job=fullfile(root,'outputs',['short65_wake_' char(java.util.UUID.randomUUID())]);
    runCM011Python({fullfile(root,'scripts','runShortSphereWakeCase.py'),job,'--through',num2str(horizon)},[job '.log'],pulse,true);
    name=sprintf('result_%03d',round(horizon*8));ack=jsondecode(fileread(fullfile(job,[name '.json'])));
    folder=fullfile(job,name,branch);expected=ack.summary.(branch);
    assert(ack.timeDU==horizon && strcmp(ack.case,key.caseSHA256) && strcmp(cm011FileSHA256(fullfile(job,name,'comparison.json')),ack.comparison));
else
    archive=fullfile(cases,[stem '.zip']);assert(strcmp(cm011FileSHA256(archive),key.archiveSHA256));
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
g=s.balance;assert(abs(g.balance)<.005 && all(structfun(@(v)v<-8,g.window)));
if wake
    assert(g.assemblyChange<.002 && abs(g.airStorageW-g.independentStorageW)<1e-9 && ...
        abs(s.elapsedPhysicalTimeS-s.physicalStep*s.dtSeconds)<1e-13 && strcmp(s.baseSHA256,'794193ab0027408bea43686d9145e1e2e0b6d133f76ac37bab102641eca079d9'));
    if ~isempty(job)
        step=fullfile(job,branch,sprintf('step%03d',s.physicalStep));
        assert(strcmp(cm011FileSHA256(fullfile(step,sprintf('restart_%05d.csv',s.physicalStep))),s.frameSHA256) && ...
            strcmp(cm011FileSHA256(fullfile(step,sprintf('seed_%05d.csv',s.physicalStep-1))),s.previousFrameSHA256) && ...
            strcmp(cm011FileSHA256(fullfile(step,'readout','release_native',sprintf('restart_%05d.csv',s.physicalStep))),s.reactionSHA256));
    end
else
    assert(g.extra20Change<.002 && all(structfun(@(v)v<-8,g.extraWindow)));
end
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
if wake
    condition.caseLabel=sprintf('Short Test 65 · 尾迹瞬时场 · tU/D = %g，ΔtU/D = %g',s.timeDU,s.dtDU);
    condition.elapsedPhysicalTimeS=s.elapsedPhysicalTimeS;condition.dtSeconds=s.dtSeconds;condition.initialStateSHA256=s.seedSHA256;
end
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
if wake
    result.shortSphere.timeDU=s.timeDU;result.shortSphere.dtDU=s.dtDU;
    id.sourceFrame.kind='SU2 8.6.0 三维 BDF1 尾迹瞬时场';id.sourceFrame.seedSHA256=s.seedSHA256;
    id.sourceFrame.previousFrameSHA256=s.previousFrameSHA256;id.sourceFrame.physicalStep=s.physicalStep;
    id.scope='给定小幅尾迹扰动后的瞬时解；七站比较使用实验参考值，支管与实测入口湍动未建模。';
end
result.evidenceIdentity=id;report=s;report.coverage=coverage;report.action=mode;report.jobDirectory=job;
% 模块三：已发布的尾迹剖面只在原网格、工况和同一物理终帧上显示。
if wake && s.timeDU==8
    path=fullfile(cases,'wakeResolution.json');
    assert(strcmp(cm011FileSHA256(path),'af201413a60e6388d05d3ed255820d1b89152502d80c5e2dd9633d3c60c26485'), ...
        'RAH:Short:Resolution','银球尾迹空间诊断文件校验失败。');
    diagnostic=jsondecode(fileread(path));
    selected=diagnostic.(branch);
    if strcmp(selected.frameSHA256,s.frameSHA256) && strcmp(diagnostic.meshSHA256,s.meshSHA256) && ...
            strcmp(diagnostic.caseSHA256,s.caseSHA256) && selected.timeS==s.elapsedPhysicalTimeS
        report.resolution=selected;
    end
end
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
        reason='稳态末窗、追加步热稳定及同态热量账通过';
        if wake
            t.method(:)="三维 BDF1 不可压层流及逐节点强温界标量";
            t.heatFluxMethod(:)="同一物理步固定端点装配壁反力，按节点面积守恒积分；释放壁行独立复核";
            reason='真实时间层残差、同态壁反力与含 BDF1 储能的收支通过';
        end
        out=struct('accepted',true,'reason',reason,'wall',t);
    end
end
function discard(folder)
% 模块三：只释放本次展开的四个已登记小文件。
for name={'wall.csv','stations.csv','body.json','summary.json'}
    path=fullfile(folder,name{1});if isfile(path),delete(path);end
end
rmdir(folder);
end
