function report=exportNasaBoattailStudy(folder)
%EXPORTNASABOATTAILSTUDY 将真实来源对照、名义几何和候选压力比接回统一科研工作台。
% 模块一：同窗五步骤切换三组已冻结研究；保留给定/未知及相关式身份。
% 模块二：仅在直尾锥映射条件候选压力比，原模型其余区域和全部热场保持未知。
% 模块三：导出逐点实验比较、来源、实际面元和公共色标；不将覆盖审计当物理解。
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'tests'));
assert(nargin==1&&~exist(folder,'dir'),'NasaBoattailView:Overwrite','须使用新的研究导出目录。');mkdir(folder);
comparison=compareNasaBoattailPressure();entries=struct([]);range=[.30 .75];
for k=1:3
    row=comparison.cases(k);angle=row.angleDeg;[m,t]=buildNasaBoattailBody(angle,48);
    relative=sprintf('cases/beta%d',angle);out=fullfile(folder,relative);mkdir(out);
    exportRocketTestGeometryView({m},fullfile(out,'geometry'),{sprintf('NASA名义15°鼻锥—柱体—%d°尾锥；无支杆/砂粒带',angle)});
    data=buildAerothermalViewerData(m);n=size(m.F,1);tail=strcmp(t.faceRegion,'boattail');
    data.fields=struct([]);data.defaultField='pressureRatio';
    candidateScope=['NASA TN D-2761声速流线压力近似；仅条件压力比，不解分离/横流/热流；' ...
        '名义几何不含支杆及砂粒带，未取得生产资格'];
    values=nan(n,1);values(tail)=row.predictionRatio;status=double(tail);
    f=field('pressureRatio','候选压力比 p_{bt}/p_c','1',values,status,row.candidate.method,candidateScope);
    f.stateNote='面积只表示该条件压力比的计算位置；不是生产合格压力/热流覆盖。';data.fields=f;
    for spec={{'pressure','绝对表面压力','kPa','未取得同工况当地柱面绝对静压；不以远场静压替代'}, ...
              {'heatFlux','对流热流密度','kW/m²','本研究仅比较压力关联，未执行适用热流求解'}, ...
              {'edgeTemperature','外缘温度','K','声速流线代理不是实际外缘热力状态'}, ...
              {'recoveryTemperature','恢复温度','K','没有受支持外缘/边界层解，不补造恢复温度'}, ...
              {'wallInput','实验壁温输入','K','原报告所选运行的壁温尚未重建'}, ...
              {'wallSolved','求解壁温','K','未执行独立壁温求解，不用给定值替代'}}
        s=spec{1};data.fields(end+1)=field(s{1},s{2},s{3},nan(n,1),zeros(n,1),['未求解：' s{4}],s{4}); %#ok<AGROW>
    end
    caption=sprintf('NASA TN D-2761鼻型1，α=0°，尾锥%d°；压力比关联对原图实验符号；绝对柱压/Tw未知',angle);
    data.caseNote=caption;data.resultCondition=struct('caption',caption,'angleOfAttackDeg',0, ...
        'freestreamMachNominalFromReport',5.98,'localCylinderMachFromPressureRatio',5.80, ...
        'wallTemperatureInputK',NaN,'absoluteCylinderPressurePa',NaN,'fieldPressureIsDimensionlessCandidate',true);
    data.flowDirection=[-1 0 0];data.coordinateNote='SI / m · X鼻尖正向 · 尾底X=0 · 灰色表示该场量未知';
    data.comparisonColorRanges=struct('pressureRatio',range);
    data.regions=struct([]);data.heatingGaps=struct([]);data.region=t.faceRegion;
    data.researchIdentity=struct('sourceManifestSHA256',comparison.sourceManifestSHA256, ...
        'method',row.candidate.method,'experimentalComparisonPerformed',true, ...
        'completeNSRANSExperimentalValidation',false,'productionQualified',false, ...
        'experimentalSymbolsAssignedToFacets',false,'newProductionSolvedAreaM2',0);
    data.candidateRatioEvaluatedFacetAreaM2=sum(m.area(tail));data.sourceGeometry=rmfield(t,'faceRegion');
    results=fullfile(out,'results');mkdir(results);
    for name={'index.html','viewer.css','viewer.js','section.js','export.js'},copyfile(fullfile(root,'web',name{1}),fullfile(results,name{1}));end
    json=jsonencode(data,'ConvertInfAndNaN',true);write(fullfile(results,'results.json'),json);write(fullfile(results,'case-data.js'),['window.ROCKET_CASE=' json ';' newline]);
    save(fullfile(results,'snapshot.mat'),'data');
    entry=struct('label',row.label,'angleDeg',angle,'geometryURL',[relative '/geometry/index.html'], ...
        'resultsURL',[relative '/results/index.html'],'resultsJSONURL',[relative '/results/results.json'], ...
        'geometry',rmfield(t,'faceRegion'),'comparison',row,'faceCount',n, ...
        'candidateRatioEvaluatedFacetAreaM2',data.candidateRatioEvaluatedFacetAreaM2);
    if isempty(entries),entries=entry;else,entries(end+1)=entry;end %#ok<AGROW>
end
payload=struct('kind','NASA-TND2761-normalized-pressure-study-v1','cases',entries, ...
    'comparison',comparison,'commonRatioRange',range,'readOnlyFrozenStudy',true, ...
    'sourceFigureURL','source_figure9.png','experimentalUncertaintyKnown',false, ...
    'newProductionSolvedAreaM2',0,'browserVisualOrWebGLRuntimeVerified',false);
json=jsonencode(payload,'ConvertInfAndNaN',true);write(fullfile(folder,'boattail-study-data.json'),json);write(fullfile(folder,'boattail-study-data.js'),['window.NASA_BOATTAIL_STUDY=' json ';' newline]);
for name={'nasa-boattail.html','nasa-boattail.js','rocket-sensitivity.css'}
    destination=name{1};if endsWith(destination,'.html'),destination='index.html';end
    copyfile(fullfile(root,'web',name{1}),fullfile(folder,destination));
end
copyfile(fullfile(root,'tests','data','nasa_tnd2761','raw','figure9.png'),fullfile(folder,'source_figure9.png'));
fid=fopen(fullfile(folder,'boattail_pressure_comparison.csv'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'point_id,beta_deg,xprime_over_D,pixel_x,pixel_y,measured_pbt_over_pc,read_half_width,candidate_pbt_over_pc,relative_difference,orifice_id,absolute_pressure_Pa,wall_temperature_K\n');
for k=1:3
    c=entries(k).comparison;
    for j=1:5,p=c.points(j);fprintf(fid,'%s,%g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,,,\n', ...
        p.id,c.angleDeg,p.sourceFigureXprimeOverD,p.pixelXY(1),p.pixelXY(2),p.pressureRatio,p.ratioDigitizationHalfWidth,c.predictionRatio,c.relativeDifference(j));end
end
report=struct('status','exported-three-normalized-pressure-correlation-experiment-cases', ...
    'htmlFile',fullfile(folder,'index.html'),'experimentalSymbols',15,'sourceFigureVisuallyChecked',true, ...
    'browserVisualOrWebGLRuntimeVerified',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'study_evidence.mat'),'report','comparison');write(fullfile(folder,'study_summary.json'),jsonencode(report));
end
function f=field(key,label,unit,values,status,method,scope)
% 模块用途：每个场量独立保存状态、方法、适用域及缺失原因。
n=numel(values);reasons=repmat({scope},n,1);reasons(status>0)={''};
f=struct('key',key,'label',label,'unit',unit,'values',values,'status',status, ...
    'methodNames',{{method}},'methodIndex',zeros(n,1),'scope',scope, ...
    'scopeNames',{{scope}},'scopeIndex',zeros(n,1),'missingReason',{reasons},'stateNote',scope);
end
function write(path,text)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
