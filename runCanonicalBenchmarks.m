function report=runCanonicalBenchmarks(outputFolder,options)
%RUNCANONICALBENCHMARKS 经典模型的完整工程链路和公开理论表格校核。
% 用法：runCanonicalBenchmarks('output/canonical')。
% 模块用途：建立可重复工况，按解析区域审计每个场量，导出三维结果。
% 覆盖率仅说明计算状态，不代表误差验证；不自动放宽物理适用范围。
if nargin<1,outputFolder=fullfile(fileparts(mfilename('fullpath')),'output','canonical');end
if nargin<2,options=struct;end
resolutions=getOption('resolutions',[24 48]);alphas=getOption('alphas',[0 20]);
kinds=getOption('kinds',{'sphere','cone','coneCylinder','hemisphereCylinder','bluntedCone'});
exportWeb=getOption('exportWeb',true);
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end
report=struct('version',1,'evidenceType','theoretical table check and engineering coverage audit', ...
    'experimentalValidation',false,'shockTable',verifyPublishedShockTable(),'cases',struct([]));
assert(report.shockTable.passed,'RocketAeroHeat:Benchmark:ShockTable','公开正激波表格校核失败');
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.externalFlow.normalModel=getOption('normalModel','validatedSections');
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.beta=0;cfg.flight.Twall=300;
for j=1:numel(kinds)
    for n=resolutions
        [geometry,meta]=buildCanonicalBody(kinds{j},struct('resolution',n));
        evalc('geometry=calcFaceGeometry(geometry);geometry=fixNormalDirection(geometry,cfg);');
        evalc('geometry=buildSurfaceCoordinates(geometry,cfg);');
        geometry=analyzeLocalGeometry(geometry,cfg);
        for alpha=alphas
            caseCfg=cfg;caseCfg.flight.alpha=alpha;
            % 已知球头几何半径为输入真值；不把半径当作外流解。
            if isfinite(meta.noseRadius),caseCfg.stagnation.physicalNoseRadius=meta.noseRadius;end
            timer=tic;mesh=calcFlightStationState(geometry,caseCfg.flight,caseCfg);
            data=buildAerothermalViewerData(mesh);
            entry=struct('name',sprintf('%s_n%d_a%g',kinds{j},n,alpha), ...
                'geometry',rmfield(meta,'faceRegion'),'flight',caseCfg.flight, ...
                'faces',size(mesh.F,1),'seconds',toc(timer),'area',sum(mesh.area), ...
                'fields',struct,'regions',struct([]),'heatingGaps',data.heatingGaps, ...
                'axisymmetricMarch',mesh.axisymmetricMarch, ...
                'externalNormalModel',cfg.externalFlow.normalModel, ...
                'normalReconstruction',struct('status',mesh.boundaryLayerEdge.normalReconstruction.status, ...
                'validatedAreaFraction',sum(mesh.area(mesh.boundaryLayerEdge.normalReconstruction.appliedFace))/sum(mesh.area), ...
                'maximumCorrectionAngleDeg',max(mesh.boundaryLayerEdge.normalReconstruction.correctionAngleDeg)), ...
                'noseRadiusEstimate',mesh.stagnation.noseRadiusEstimate, ...
                'streamlineSummary',mesh.surfaceStreamlines.summary);
            for k=1:numel(data.fields)
                field=data.fields(k);entry.fields.(field.key)=summarize(field,true(size(mesh.F,1),1));
            end
            labels=unique(meta.faceRegion,'stable');
            for k=1:numel(labels)
                mask=strcmp(meta.faceRegion,labels{k});
                row=struct('name',labels{k},'faces',nnz(mask),'area',sum(mesh.area(mask)),'fields',struct);
                for f=1:numel(data.fields)
                    field=data.fields(f);row.fields.(field.key)=summarize(field,mask);
                end
                if k==1,entry.regions=row;else,entry.regions(end+1)=row;end
            end
            if isempty(report.cases),report.cases=entry;else,report.cases(end+1)=entry;end
            if exportWeb && n==max(resolutions)
                exportAerothermalWeb(mesh,fullfile(outputFolder,entry.name));
            end
            fprintf('%s: %d faces, computed heat area %.2f%%\n',entry.name,entry.faces,100*entry.fields.heatFlux.computedFraction);
            % 每个完成的算例即保存，长任务中断时仍有可信的阶段证据。
            writeJson(fullfile(outputFolder,'report.json'),report);
        end
    end
end
save(fullfile(outputFolder,'report.mat'),'report');
writeMarkdown(fullfile(outputFolder,'report.md'),report);

    function value=getOption(key,fallback)
        value=fallback;if isfield(options,key),value=options.(key);end
    end
    function s=summarize(field,mask)
        a=mesh.area(mask);v=field.values(mask);status=field.status(mask);total=sum(a);
        valid=status==1 & isfinite(v);
        s=struct('unit',field.unit,'computedFraction',sum(a(status==1))/total, ...
            'prescribedFraction',sum(a(status==2))/total,'missingFraction',sum(a(status==0))/total, ...
            'computedAreaWeightedMean',NaN,'computedMinimum',NaN,'computedMaximum',NaN, ...
            'methods',{field.methodNames(unique(field.methodIndex(mask)+1))});
        if any(valid)
            s.computedAreaWeightedMean=sum(a(valid).*v(valid))/sum(a(valid));
            s.computedMinimum=min(v(valid));s.computedMaximum=max(v(valid));
        end
    end
end
function writeJson(path,data)
% 模块：保存机器可读明细，未解数值序列化为 null。
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(data,'PrettyPrint',true,'ConvertInfAndNaN',true));
end
function writeMarkdown(path,report)
% 模块：生成科研校核摘要，按区域保留缺口；不混淆验证种类。
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# 经典模型校核与区域覆盖\n\n');
fprintf(fid,'NACA 1135 表 II（印刷页 25）四行共 20 个正激波数值在半个印刷末位精度内通过：%d。此项是理论表格校核，不是实验。\n\n',report.shockTable.passed);
fprintf(fid,'工况：M=4，海拔 20 km，输入壁温 300 K，零侧滑；长度单位 m。所有覆盖率按三角面面积加权。给定压力不计入计算覆盖。\n\n');
fprintf(fid,'| 算例 | 面数 | 压力计算/给定 | 外缘压力计算 | 热流计算 |\n|---|---:|---:|---:|---:|\n');
for i=1:numel(report.cases)
    e=report.cases(i);f=e.fields;
    fprintf(fid,'| %s | %d | %.2f%% / %.2f%% | %.2f%% | %.2f%% |\n',e.name,e.faces,100*f.pressure.computedFraction,100*f.pressure.prescribedFraction,100*f.edgePressure.computedFraction,100*f.heatFlux.computedFraction);
end
fprintf(fid,'\n不同分辨率结果用于暴露离散敏感性；本报告不据此声称热流已经网格收敛。区域名称来自解析几何真值。\n');
for i=1:numel(report.cases)
    e=report.cases(i);fprintf(fid,'\n## %s\n\n| 解析区域 | 面积 m² | 热流计算 | 热流缺失 |\n|---|---:|---:|---:|\n',e.name);
    for j=1:numel(e.regions)
        r=e.regions(j);fprintf(fid,'| %s | %.7g | %.2f%% | %.2f%% |\n',r.name,r.area,100*r.fields.heatFlux.computedFraction,100*r.fields.heatFlux.missingFraction);
    end
    for j=1:numel(e.heatingGaps)
        g=e.heatingGaps(j);fprintf(fid,'\n- %s：总面积的 %.2f%%。',g.reason,100*g.areaFraction);
    end
    fprintf(fid,'\n');
end
end
