function report=exportShockLayerStudy(parabolicFile,conicFile,translatedFile,folder)
%EXPORTSHOCKLAYERSTUDY 在统一科研框架中查看真实箭体及已保存的共同激波层列。
% 模块一：严格读取三个已保存MAT，不重新拟合、不重算激波或热流。
% 模块二：整箭与鼻部底图均来自原V/F，鼻部是明确的开放几何片。
% 模块三：层内状态、候选生成体表和真实面元距离分列；未知热量不补零。
% 模块四：导出五步骤、三维站点拾取、法向剖面、逐场量方法和CSV。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);
assert(~exist(folder,'dir'),'ShockLayer:StudyOverwrite','研究视图必须导出到新目录。');
p=load(parabolicFile);c=load(conicFile);t=load(translatedFile);
assert(p.report.actualFaces==19752&&c.report.actualFaces==19752&&t.report.originalFaces==19752);
assert(isequaln(p.geometry,t.geometry)&&isequaln(c.geometry,t.geometry));
assert(~p.bodyFit.productionQualified&&~c.bodyFit.productionQualified&&~t.bodyFit.productionQualified ...
    &&~t.bodyFit.boundaryLayerEdgeAvailable&&~t.report.heldOutBodyMatchingPassed, ...
    'ShockLayer:StudyQualification','此导出只允许尚未通过匹配的已保存研究，不可伪造正式资格。');
assert(p.report.newProductionSolvedAreaM2==0&&c.report.newProductionSolvedAreaM2==0&&t.report.newProductionSolvedAreaM2==0);

%% 一、逐列保留实际法向坐标与方法；绝不把这些线称为表面流线
groups={p.cases{2}.column,p.cases{1}.column,t.bodyFit.result.columns,t.heldOut.columns};
ids={'NASA7_profile','perfect_gas_profile','fitted_9','held_out_36'};
captions={'NASA7 · 层内法向剖面','定比热 · 层内法向剖面','七参数形状 · 9个拟合列','七参数形状 · 36个独立检查列'};
roles={'给定共同抛物激波的数学验证列；没有真实体形匹配', ...
    '相同几何与自由流的定比热对照；没有真实体形匹配', ...
    '原三角面参与拟合的列；最大残差仍超过10微米', ...
    '未参与任何形状优化的检查列；不能以拟合点改善代替整体匹配'};
viewCases=cell(4,1);stationCount=0;
for k=1:4
    columns=groups{k};if isstruct(columns),columns=num2cell(columns);end
    tubes=cell(numel(columns),1);
    for j=1:numel(columns)
        col=columns{j};assert(~col.productionQualified&&~col.boundaryLayerEdgeAvailable&&~col.bodyMatched);
        ns=repmat(-col.shock.normal',numel(col.eta),1);
        mapped=[];residual=[];face=[];bary=[];
        if k==3,mapped=t.bodyFit.result.geometricMapping;residual=t.bodyFit.result.residualM(j);end
        if k==4,mapped=t.heldOut.geometricMapping;residual=t.heldOut.residualM(j);end
        if ~isempty(mapped),face=mapped.originalFaceID(j);bary=mapped.barycentric(j,:);end
        methods=struct('pressure',col.state.pressureMethod,'temperature',col.state.temperatureMethod, ...
            'speed',col.state.speedMethod,'entropy','回溯同一共同激波的来源入口；沿真实流线保持熵', ...
            'distance','正向质量积分 + 完整非正交Jacobian法向距离反算');
        tubes{j}=struct('eta',col.eta,'position',col.position,'pickNormals',ns, ...
            'normalDistanceM',col.normalDistanceM,'pressurePa',col.state.pressure, ...
            'temperatureK',col.state.temperature,'densityKgM3',col.state.density, ...
            'speedMS',col.state.speed,'entropyJPerKgK',col.state.entropy,'jacobianRatio',col.jacobianRatio, ...
            'sourceShockCoordinateM',col.state.sourceShockCoordinate,'gasMode',col.shock.gas.thermalPropertyMode, ...
            'methods',methods,'bodyResidualM',residual,'originalFaceID',face,'barycentric',bary, ...
            'stateRole','层内及候选生成体表；不是正式壁面或边界层外缘', ...
            'heatFluxWPerM2',[],'wallTemperatureInputK',[],'wallTemperatureSolvedK',[], ...
            'edgeTemperatureK',[],'recoveryTemperatureK',[],'productionQualified',false);
        stationCount=stationCount+numel(col.eta);
    end
    viewCases{k}=struct('id',ids{k},'caption',captions{k},'scope',roles{k},'columns',{tubes}, ...
        'usedDuringFit',k==3,'bodyMatchToleranceM',t.bodyFit.bodyMatchToleranceM,'productionQualified',false);
end

%% 二、原整箭拓扑与鼻部视图分开；只选择原三角形，不新造或平滑实体
V=double(t.geometry.V);F=double(t.geometry.F);edges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
[~,~,index]=unique(sort(edges,2),'rows');counts=accumarray(index,1);directions=accumarray(index,sign(edges(:,2)-edges(:,1)));
assert(all(counts==2)&&all(directions==0),'ShockLayer:StudyGeometry','完整原箭体必须保持闭合定向。');
centers=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
originalIDs=find(centers(:,1)>=max(V(:,1))-.065);sub=F(originalIDs,:);[vertices,~,indices]=unique(sub(:));
focus=struct('vertices',V(vertices,:),'faces',reshape(indices,size(sub))-1,'originalFaceIDs',originalIDs, ...
    'closed',false,'scope','原三角面鼻部开放片，仅用于放大查看，不是新实体或新增解面积');
geometry=struct('vertices',V,'faces',F-1,'closed',true,'boundaryEdges',nnz(counts==1), ...
    'nonManifoldEdges',nnz(counts>2),'focus',focus);

%% 三、形状匹配账与场量身份随数据导出；来源由原MAT内容SHA绑定
fitSummary=struct('names',{{'初始抛物面','四参数抛物面','五参数二次曲线','七参数拟合点','七参数独立检查点'}}, ...
    'rmsM',[p.report.initialBodyRMSM p.report.finalBodyRMSM c.report.finalBodyRMSM t.report.fitRMSM t.report.heldOutRMSM], ...
    'maximumM',[max(abs(p.bodyFit.initial.residualM)) p.report.maximumBodyResidualM c.report.maximumBodyResidualM ...
        t.report.maximumFitResidualM t.report.maximumHeldOutResidualM],'toleranceM',1e-5,'experimentalValidation',false);
sources=struct('parabolic',fileRecord(parabolicFile),'conic',fileRecord(conicFile),'translated',fileRecord(translatedFile));
data=struct('kind','shock-layer-study-v1','geometry',geometry,'cases',{viewCases}, ...
    'freestream',struct('Mach',double(t.fs.Mach),'altitudeM',double(t.fs.altitude),'alphaDeg',double(t.fs.alpha), ...
    'pressurePa',double(t.fs.pressure),'temperatureK',double(t.fs.temperature)), ...
    'fitSummary',fitSummary,'sources',sources,'productionQualified',false,'newProductionSolvedAreaM2',0, ...
    'coordinateConvention','SI · 主轴X · 鼻尖+X · 零攻角来流−X', ...
    'methodScope','NASA TP-3018薄层近似；未通过真实体形匹配，驻点整体/跨声速/热流和壁温仍未解');
mkdir(folder);web=fullfile(folder,'workbench');mkdir(web);geoweb=fullfile(web,'geometry');mkdir(geoweb);
for name={'geometry.html','geometry.css','geometry.js'},copyfile(fullfile(root,'web',name{1}),fullfile(geoweb,name{1}));end
html=fileread(fullfile(geoweb,'geometry.html'));
html=strrep(html,'</body>','<script src="../shock-layer-data.js"></script><script src="shock-layer-geometry.js"></script></body>');
writeUtf8(fullfile(geoweb,'index.html'),html);
copyfile(fullfile(root,'web/shock-layer-geometry.js'),fullfile(geoweb,'shock-layer-geometry.js'));
for name={'rocket-sensitivity.css','shock-layer-study.css','shock-layer-study.js'},copyfile(fullfile(root,'web',name{1}),fullfile(web,name{1}));end
copyfile(fullfile(root,'web/shock-layer-study.html'),fullfile(web,'index.html'));
text=jsonencode(data);writeUtf8(fullfile(web,'shock-layer-data.json'),text);
writeUtf8(fullfile(web,'shock-layer-data.js'),['window.SHOCK_LAYER_STUDY=' text ';' newline]);
report=struct('status','exported-saved-unmatched-shock-layer-study-with-explicit-state-roles', ...
    'cases',4,'columns',sum(cellfun(@(x)numel(x.columns),viewCases)),'stations',stationCount, ...
    'originalFaces',size(F,1),'focusFaces',numel(originalIDs),'newPhysicalIntegrations',0, ...
    'newShapeFits',0,'productionQualified',false,'newProductionSolvedAreaM2',0, ...
    'dataSHA256',fileSHA(fullfile(web,'shock-layer-data.json')),'workbench',fullfile(web,'index.html'));
writeUtf8(fullfile(folder,'shock_layer_export.json'),jsonencode(report));disp(report);
end

function value=fileRecord(file)
value=struct('path',strrep(file,'\','/'),'sha256',fileSHA(file));
end
function value=fileSHA(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),b=fread(fid,1024^2,'*uint8');if ~isempty(b),md.update(typecast(b,'int8'));end,end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
function writeUtf8(file,value)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
