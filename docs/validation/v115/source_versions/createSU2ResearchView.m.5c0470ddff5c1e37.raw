function view=createSU2ResearchView(host,result)
%CREATESU2RESEARCHVIEW 同一主工作台内查看内嵌NS研究场，沿用离线WebGL。
% 模块一：保留原面元、原始q及独立保守热负荷，不将有限值算入合格面积。
% 模块二：固定色标可编辑；拾取同时给出p/q/温度角色和方法。
% 模块三：导出当前快照CSV和身份JSON，不重算且不给未知壁温补零。
grid=uigridlayout(host,[4 1]);grid.RowHeight={67,38,'1x',130};grid.Padding=[0 0 0 0];
% 模块四：固定显示当前场的保存时刻和真实终态，不能用末次历史冒充云图时刻。
auditText=uilabel(grid,'Text',strjoin(describeSU2Convergence(result),newline),'FontSize',11,'WordWrap','on');auditText.Layout.Row=1;
bar=uigridlayout(grid,[1 6]);bar.Layout.Row=2;bar.ColumnWidth={'1x',100,100,110,110,110};
bar.Padding=[0 0 0 0];bar.RowHeight={'1x'};
keys={'pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw','conservativeHeatLoadDensityRaw','wallInput','wallSolved'};
labels={'压力（研究）','物面平均热流（研究）','流体节点温度（研究）','保守热负荷密度（研究）','给定壁温（输入）','求解壁温（未知）'};
% 绝热的零热流是边界输入；不能把整面零云图称作自主热流预测。
isAdiabatic=isfield(result,'heatFluxIsBoundaryInput') && result.heatFluxIsBoundaryInput;
n=size(result.F,1);wallInput=result.wallTemperatureInput;
if isscalar(wallInput),wallInput=repmat(wallInput,n,1);end
assert(isequal(size(wallInput),[n 1]),'给定壁温与原面数量不一致。');
thermalMethod='用户给定等温壁面边界';thermalCaption='给定等温壁';
if isfield(result,'thermalBoundary'),thermalMethod=result.thermalBoundary.method;thermalCaption=result.thermalBoundary.mode;end
if isAdiabatic
    labels{2}='热流（绝热给定0）';labels{4}='热负荷密度（绝热给定0）';labels{5}='给定壁温（绝热未指定）';
elseif isfield(result,'thermalBoundary') && strcmp(result.thermalBoundary.mode,'axialTemperature')
    labels{5}='给定壁温分布（输入）';
end
units={'Pa','W/m²','K','W/m²','K','K'};
select=uidropdown(bar,'Items',labels,'ItemsData',keys,'Value',keys{1});
low=uieditfield(bar,'numeric');high=uieditfield(bar,'numeric');
low.Tooltip='固定色标下限（当前变量单位）';high.Tooltip='固定色标上限（当前变量单位）';
apply=uibutton(bar,'Text','应用固定色标','ButtonPushedFcn',@(~,~)render(false));
auto=uibutton(bar,'Text','按当前值范围','ButtonPushedFcn',@(~,~)render(true));
export=uibutton(bar,'Text','导出当前CSV','ButtonPushedFcn',@exportDialog);
panel=uipanel(grid,'BorderType','none');panel.Layout.Row=3;
geometry=createAerothermalGeometryView(panel,result,struct('su2Section',true));
probe=uitextarea(grid,'Editable','off');probe.Layout.Row=4;
oldCallback=geometry.html.DataChangedFcn;geometry.html.DataChangedFcn=@receive;
select.ValueChangedFcn=@(~,~)render(true);
sliceX=mean([min(result.V(:,1)),max(result.V(:,1))]);selectedFace=1;sectionState=[];
view=struct('geometry',geometry,'select',select,'low',low,'high',high,'apply',apply, ...
    'probe',probe,'export',export,'exportTo',@exportTo,'inspectFace',@inspectFace, ...
    'setSectionX',@setSectionX,'exportSectionTo',@exportSectionTo,'getSection',@getSection,'auditText',auditText);
render(true);inspectFace(1,[]);

    function render(resetRange)
        key=select.Value;i=find(strcmp(keys,key));n=size(result.F,1);
        if strcmp(key,'wallInput'),values=wallInput;
        elseif strcmp(key,'wallSolved'),values=nan(n,1);
        else,values=result.(key);end
        if resetRange
            v=values(isfinite(values));if isempty(v),range=[0 1];else,range=[min(v) max(v)];end
            if range(2)<=range(1),range=range+[-1 1]*max(1,abs(range(1))*.01);end
            low.Value=range(1);high.Value=range(2);
        end
        if low.Value>=high.Value,probe.Value={'色标上限必须大于下限。'};return;end
        caption=sprintf('SU2 8.6.0 · M %.4g / α %.4g° / β %.4g° · 研究快照，未获物理资格', ...
            result.flight.Mach,result.flight.alpha,result.flight.beta);
        updateAerothermalGeometryView(geometry,result,caption,struct('resetCamera',false));
        s=geometry.html.UserData;s.payload.faceField=struct('values',values, ...
            'range',[low.Value high.Value],'name',labels{i},'unit',units{i},'productionQualified',false);
        s.payload.coordinateNote='原三角面 · SI · laminar NS / ideal gas / Sutherland / Pr0.72 · 合格新增面积0';
        status=1;method=[result.solverMethod '；' result.projectionMethod];
        if strcmp(key,'wallInput'),status=2;method=thermalMethod;if isAdiabatic,status=0;end
        elseif strcmp(key,'wallSolved'),status=0;method='尚未求解壁体温度';end
        if isAdiabatic && ismember(key,{'heatFluxRaw','conservativeHeatLoadDensityRaw'}),status=2;method=thermalMethod;end
        s.payload.sectionView=struct('x',sliceX,'faceIndex',selectedFace-1,'fieldKey',key, ...
            'valueStatus',status,'method',method,'scope','完整原面研究快照；未获网格/物理资格', ...
            'condition',sprintf('M%.8g alpha%.8g beta%.8g; %s',result.flight.Mach,result.flight.alpha,result.flight.beta,thermalCaption));
        sectionState=[];
        geometry.html.UserData=s;geometry.html.Data=s.payload;
    end
    function receive(src,event)
        message=src.Data;
        if isstruct(message) && isfield(message,'kind') && strcmp(message.kind,'geometry-face-probe')
            if message.revision==src.UserData.revision,inspectFace(message.faceIndex+1,message.position);end
        elseif isstruct(message) && isfield(message,'kind') && strcmp(message.kind,'geometry-section-export')
            if acceptSection(message,src)
                [name,folder]=uiputfile('*.csv','导出当前截面研究值');if ~isequal(name,0),exportSectionTo(fullfile(folder,name));end
            end
        else
            oldCallback(src,event);
            if isstruct(message) && isfield(message,'kind') && strcmp(message.kind,'geometry-rendered'),acceptSection(message,src);end
        end
    end
    function yes=acceptSection(message,src)
        yes=isfield(message,'revision') && message.revision==src.UserData.revision && isfield(message,'sectionInfo');
        if ~yes,return;end
        info=message.sectionInfo;yes=isstruct(info) && info.revision==src.UserData.revision && strcmp(info.fieldKey,select.Value);
        if ~yes,return;end
        validateattributes(info.x,{'numeric'},{'scalar','finite','>=',min(result.V(:,1)),'<=',max(result.V(:,1))});
        sectionState=info;sliceX=info.x;
    end
    function inspectFace(i,position)
        validateattributes(i,{'numeric'},{'scalar','integer','>=',1,'<=',size(result.F,1)});
        selectedFace=i;
        if nargin<2,position=mean(result.V(result.F(i,:),:),1);end
        if ~isempty(position)
            validateattributes(position,{'numeric'},{'numel',3,'finite','real'});
            sliceX=position(1);render(false);
        end
        wallText=sprintf('给定Tw %.8g K',wallInput(i));
        if isnan(wallInput(i)),wallText='给定Tw：未指定（绝热壁）';end
        probe.Value={sprintf('原面元 %d · 面积 %.8g m² · %s',i,result.faceArea(i),result.resultStatus); ...
            sprintf('p %.8g Pa；物面平均 q %.8g W/m²；保守热负荷密度 %.8g W/m²（口径不同）', ...
                result.pressureRaw(i),result.heatFluxRaw(i),result.conservativeHeatLoadDensityRaw(i)); ...
            sprintf('流体节点平均T %.8g K；%s；求解壁温/Te/Tr 未知或未定义', ...
                result.fluidNodeTemperatureMeanRaw(i),wallText); ...
            ['算法：' result.solverMethod '；热边界：' thermalMethod '；映射：' result.projectionMethod]; ...
            sprintf('全物面原始q积分 %.8g W；原生SU2负荷 %.8g W。进程退出0不等于收敛/精度通过。', ...
                result.targetHeatPowerRaw,result.conservativeHeatLoadPowerRaw)};
    end
    function setSectionX(x)
        validateattributes(x,{'numeric'},{'scalar','real','finite','>=',min(result.V(:,1)),'<=',max(result.V(:,1))});
        sliceX=x;render(false);
    end
    function s=getSection()
        % 使用嵌套函数读取实时共享状态；匿名闭包会捕获创建时的空值。
        s=sectionState;
    end
    function exportSectionTo(file)
        assert(~isempty(sectionState) && sectionState.revision==geometry.html.UserData.revision && ...
            strcmp(sectionState.fieldKey,select.Value),'请等待当前截面完成绘制后再导出。');
        fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',sectionState.csv);clear guard;
        meta=struct('x',sectionState.x,'field',sectionState.fieldKey,'revision',sectionState.revision, ...
            'sourceRevision',result.sourceRevision,'solveId',result.solveId,'flight',result.flight, ...
            'productionQualified',false,'method',geometry.html.UserData.payload.sectionView.method);
        if isfield(result,'convergenceAudit'),meta.convergenceAudit=result.convergenceAudit;end
        meta.processStatus=struct('exitCode',result.receipt.exitCode,'timedOut',result.receipt.timedOut,'cancelled',result.receipt.cancelled);
        fid=fopen([file '.json'],'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(meta));
    end
    function exportDialog(~,~)
        [name,folder]=uiputfile('*.csv','导出当前研究快照');if isequal(name,0),return;end
        exportTo(fullfile(folder,name));
    end
    function exportTo(file)
        n=size(result.F,1);t=table((1:n)',result.faceArea,result.pressureRaw,result.heatFluxRaw, ...
            result.conservativeHeatLoadDensityRaw,result.fluidNodeTemperatureMeanRaw, ...
            wallInput,nan(n,1),false(n,1),repmat(isAdiabatic,n,1), ...
            'VariableNames',{'originalFaceID','area_m2','pressureRaw_Pa','heatFluxRaw_W_m2', ...
            'nativeLoadDensity_W_m2','fluidNodeTemperatureMean_K','givenWallTemperature_K', ...
            'solvedWallTemperature_K','productionQualified','heatFluxIsBoundaryInput'});
        writetable(t,file);meta=rmfield(result,{'V','F','pressureRaw','heatFluxRaw','conservativeHeatLoadDensityRaw', ...
            'fluidNodeTemperatureMeanRaw','faceArea'});
        fid=fopen([file '.json'],'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(meta));
    end
end
