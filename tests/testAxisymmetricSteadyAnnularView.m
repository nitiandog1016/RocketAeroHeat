function report=testAxisymmetricSteadyAnnularView(folder)
%TESTAXISYMMETRICSTEADYANNULARVIEW 新隐式终态接回原查看器，严查伪时间身份。
% 模块一：读取实际环隙稳态，保持原Q/gas/墙面空间重放与面积热账。
% 模块二：伪时间不写timeS或“流体推进时间”，未知温度保持NaN/null。
% 模块三：明确混入物理time或伪造时间标志必须拒绝；旧RK2页语义保留。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','annular_view');end
source=fullfile(root,'outputs','v63_validation','annular','annular_steady_evidence.mat');
if ~exist(source,'file'),source=fullfile(root,'validation','v63','annular','annular_steady_evidence.mat');end
s=load(source,'solutions','saved');piece=s.solutions{end};result=piece.result;mesh=piece.mesh;
options=struct('name','v63隐式稳态环隙数值校核','evidenceSource',source);
view=exportAxisymmetricAnnularStudy(mesh,result,fullfile(folder,'web'),options);
data=jsondecode(fileread(view.jsonFile));
assert(contains(data.caseNote,'求稳态伪时间（非物理时间）')&&~contains(data.caseNote,'流体推进时间'));
assert(strcmp(data.resultCondition.timeKind,'pseudo-time-steady-search')&& ...
    data.resultCondition.pseudoTimeS==result.pseudoTime&&data.resultCondition.physicalTimeAdvancedS==0&& ...
    ~isfield(data.resultCondition,'timeS')&&~data.researchIdentity.pseudoTimeIsPhysical&& ...
    data.researchIdentity.physicalTimeAdvancedS==0);
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    field=data.fields(strcmp({data.fields.key},key{1}));
    assert(all(isnan(field.values))&&all(field.status==0));
end
saved=load(view.snapshotFile,'result');assert(isequaln(saved.result,result)&&~isfield(saved.result,'time'));
refusals=0;bad=result;bad.time=bad.pseudoTime;reject(@()exportAxisymmetricAnnularStudy(mesh,bad,fullfile(folder,'bad1'),options));
bad=result;bad.pseudoTimeIsPhysical=true;reject(@()exportAxisymmetricAnnularStudy(mesh,bad,fullfile(folder,'bad2'),options));
bad=result;bad.physicalTimeAdvanced=1;reject(@()exportAxisymmetricAnnularStudy(mesh,bad,fullfile(folder,'bad3'),options));
legacy=s.saved.solutions{end};
legacyView=exportAxisymmetricAnnularStudy(legacy.mesh,legacy.result,fullfile(folder,'legacy_web'));
legacyData=jsondecode(fileread(legacyView.jsonFile));
assert(contains(legacyData.caseNote,'流体推进时间')&&legacyData.resultCondition.timeS==legacy.result.time&& ...
    ~isfield(legacyData.resultCondition,'pseudoTimeS'));
report=struct('status','passed-steady-pseudo-time-view-and-real-time-legacy-view', ...
    'triangles',numel(data.faceArea),'pseudoTimeS',result.pseudoTime,'physicalTimeAdvancedS',0, ...
    'view',view,'refusals',refusals,'physicalSolution',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'annular_steady_view_evidence.mat'),'report','view','legacyView');
fid=fopen(fullfile(folder,'annular_steady_view_summary.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Implicit steady view passed: %d triangles, pseudo time separated, real RK2 legacy retained, %d refusals.\n',numel(data.faceArea),refusals);
    function reject(action)
        try,action();catch e
            assert(strcmp(e.identifier,'AxisFV:AnnularViewTime'),'应由时间身份合同明确拒绝。');
            refusals=refusals+1;return;
        end
        error('AxisSteadyTest:MissingRefusal','伪造时间身份没有拒绝。');
    end
end
