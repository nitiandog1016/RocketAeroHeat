function report=exportRocketContactStudy(folder)
%EXPORTROCKETCONTACTSTUDY 同一五步骤工作台比较实际近壁分辨率与逐格二分。
% 模块一：四档工况逐链重放，仍用既有真实三角面投影、探针与周向曲线。
% 模块二：压力/热流/气体温度统一色标；未知Te/Tr/求解Tw不补零。
% 模块三：附同墙段分布、总热功率、实际壁距及y+代理，明确未达网格独立。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(fullfile(root,'tests'));
val=fullfile(root,'outputs','v78_validation');if ~exist(fullfile(val,'contact_independent.json'),'file'),val=fullfile(root,'validation','v78');end
if nargin<1,folder=fullfile(val,'final_workbench');end
assert(~exist(folder,'dir'),'RocketSensitivityView:Overwrite','禁止覆盖已导出对照。');
a=jsondecode(fileread(fullfile(val,'contact_independent.json')));assert(startsWith(a.status,'passed-independent'));
mkdir(folder);entries=struct([]);views={};ranges=struct;
kinds={'hlle_collar72','hlle_nested144','hllc_collar72','hllc_nested144'};
labels={'HLLE · 72层','HLLE · 144层','HLLC · 72层','HLLC · 144层'};
for k=1:4
    sourceKind=a.cases(k).sourceKind;
    if k<=2
        chain=fullfile(root,'outputs','v77_validation',sourceKind,'confirmed');
        if ~exist(fullfile(chain,'rocket_axisymmetric_evidence.mat'),'file'),chain=fullfile(root,'validation','v77',sourceKind,'confirmed');end
    else,chain=fullfile(val,sourceKind,'confirmed');end
    testRocketAxisymmetricControlledCheckpoint(chain,fullfile(folder,'audits',kinds{k}));
    out=fullfile(folder,'cases',kinds{k});
    view=exportRocketAxisymmetricWorkbench(chain,out,struct('name',labels{k},'title',['RocketAeroHeat · ' labels{k}]));
    snapshot=load(view.resultView.snapshotFile,'data');data=snapshot.data;
    for fi=1:numel(data.fields)
        f=data.fields(fi);
        v=f.values(f.status>0&isfinite(f.values));if isempty(v),continue;end
        pair=[min(v),max(v)];if isfield(ranges,f.key),pair=[min(pair(1),ranges.(f.key)(1)),max(pair(2),ranges.(f.key)(2))];end
        ranges.(f.key)=pair;
    end
    study=jsondecode(fileread(view.dataFile));prefix=['cases/' kinds{k} '/'];
    chainData=load(fullfile(chain,'rocket_axisymmetric_evidence.mat'),'caseDefinition','implicit');
    normalProfiles=buildRocketWallNormalProfiles(chainData.caseDefinition,chainData.implicit);
    metrics=rmfield(a.cases(k),'wallNormalProfiles');
    entry=struct('kind',kinds{k},'label',labels{k},'study',study,'metrics',metrics, ...
        'wallNormalProfiles',normalProfiles,'profileSourceMethod',chainData.implicit.finalInfo.method, ...
        'geometryURL',[prefix study.geometryURL],'resultsURL',[prefix study.resultsURL]);
    if isempty(entries),entries=entry;else,entries(end+1)=entry;end %#ok<AGROW>
    views{k}=view; %#ok<AGROW>
end
for key=fieldnames(ranges)'
    r=ranges.(key{1});if r(1)==r(2),padding=max(abs(r(1))*.01,.001);ranges.(key{1})=r+[-padding padding];end
end
for k=1:4
    path=views{k}.resultView.snapshotFile;snapshot=load(path,'data');data=snapshot.data;
    data.comparisonColorRanges=ranges;
    data.comparisonScope='四档同工况完整零攻角数值定常结果；各变量同单位共色标；网格独立与物理精度未验收';
    json=jsonencode(data,'ConvertInfAndNaN',true);destination=fileparts(path);
    write(fullfile(destination,'results.json'),json);write(fullfile(destination,'case-data.js'),['window.ROCKET_CASE=' json ';' newline]);
    save(path,'data','-append');
end
payload=struct('kind','complete-rocket-contact-discretization-workbench-v1','cases',entries, ...
    'meshDescription','同一真实72/144层网格与物理边界；只比较已有内部通量。HLLC细网格使用已收敛粗HLLC根初猜，直接HLLE细网格初猜的失败单列。', ...
    'failedTrials',{num2cell(publishFailedTrials(a.failedTrials))},'gridChanges',struct('hlle',a.gridChanges.hlle_einfeldt,'hllc',a.gridChanges.hllc), ...
    'comparisons',a.comparisons,'comparisonColorRanges',ranges,'yPlusDefinition',a.yPlusDefinition, ...
    'physicalSolution',false,'experimentalValidation',false,'gridIndependentClaimed',false, ...
    'streamwiseRefinementPerformed',false,'newProductionSolvedAreaM2',0, ...
    'productionHeatCoverageFraction',a.productionHeatCoverageFraction,'browserVisualOrWebGLRuntimeVerified',false);
json=jsonencode(payload,'ConvertInfAndNaN',true);write(fullfile(folder,'sensitivity-data.json'),json);
write(fullfile(folder,'sensitivity-data.js'),['window.ROCKET_SENSITIVITY=' json ';' newline]);
for name={'rocket-sensitivity.html','rocket-sensitivity.css','rocket-sensitivity.js'}
    if endsWith(name{1},'.html'),destination='index.html';else,destination=name{1};end
    copyfile(fullfile(root,'web',name{1}),fullfile(folder,destination));
end
copyfile(fullfile(val,'wall_sensitivity.csv'),fullfile(folder,'wall_sensitivity.csv'));
report=struct('status','exported-four-confirmed-cases-same-host-common-color-ranges', ...
    'htmlFile',fullfile(folder,'index.html'),'caseCount',4,'BEConfirmedEndpoints',12, ...
    'browserVisualOrWebGLRuntimeVerified',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
write(fullfile(folder,'sensitivity-view-summary.json'),jsonencode(report));
end
function write(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
function trials=publishFailedTrials(trials)
% 模块用途：错误标识中的冒号以字符串值发布，避免jsondecode字段规范化改变原标识。
known={'AxisFV:Temperature','AxisSA:State'};encoded=cellfun(@matlab.lang.makeValidName,known,'UniformOutput',false);
for k=1:numel(trials)
    counts=trials(k).guardRejectCounts;names=fieldnames(counts);items=struct('identifier',{},'count',{});
    for j=1:numel(names)
        id=find(strcmp(names{j},encoded));assert(isscalar(id),'RocketContactView:GuardName','未知保护标识不能猜测反向映射。');
        items(end+1)=struct('identifier',known{id},'count',counts.(names{j})); %#ok<AGROW>
    end
    trials(k).guardRejections=num2cell(items);
end
trials=rmfield(trials,'guardRejectCounts');
end
