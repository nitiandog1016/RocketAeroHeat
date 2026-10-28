function report=findRocketCandidateTemperatureFaces(candidate,cut,gas,options)
%FINDROCKETCANDIDATETEMPERATUREFACES 在原温域内批量定位候选内部Riemann采样越温面。
% 模块一：先列出原温域内合法单元对；非法单元相邻面单列，不能冒称已检查。
% 模块二：每批至多4096面；原内核温域拒绝时二分，不修改气体温域或裁剪状态。
% 模块三：保留任何原批次真实拒绝的面，即使极近边界在子批中有舍入差异也不丢失。
% 模块四：只做方向预览，0完整NS/流场更新；显式预算拒绝不能变成漏面。
if nargin<4,options=struct;end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),{'maximumBatchFaces'})), ...
    'InterfacePreview:Options','未知界面前瞻选项。');
if ~isfield(options,'maximumBatchFaces'),options.maximumBatchFaces=4096;end
validateattributes(options.maximumBatchFaces,{'numeric'},{'scalar','integer','positive','<=',4096});
n=size(candidate,1);assert(isequal(size(candidate),[n 5])&&all(isfinite(candidate),'all')&&isreal(candidate), ...
    'InterfacePreview:State','需要当前候选的有限五量状态。');
f=cut.openFaces;ids=find(f.area>0&f.cellMinus>0&f.cellPlus>0);ni=numel(ids);
arrayBound=8*(12*ni+10*n)+2048*options.maximumBatchFaces+65536;
assert(arrayBound<=256*1024^2,'InterfacePreview:Budget','界面预览数组上界超过256MiB，不裁面。');
rho=candidate(:,1);assert(all(rho>0),'InterfacePreview:State','候选密度非正，尚不能检查界面。');
u=candidate(:,2:4)./rho;cv=gas.R/(gas.gamma-1);internal=candidate(:,5)-.5*rho.*sum(u.^2,2);
assert(all(internal>0),'InterfacePreview:State','候选内能非正，尚不能检查界面。');
T=internal./(rho*cv);valid=T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2);
eligible=valid(f.cellMinus(ids))&valid(f.cellPlus(ids));eligibleSlots=find(eligible);blocked=ids(~eligible);
failed=false(ni,1);failedT=nan(ni,1);calls=0;faceInputs=0;failedBatchCount=0;
for first=1:options.maximumBatchFaces:numel(eligibleSlots)
    local=eligibleSlots(first:min(first+options.maximumBatchFaces-1,numel(eligibleSlots)));
    queue={[1 numel(local)]};
    while ~isempty(queue)
        span=queue{end};queue(end)=[];slots=local(span(1):span(2));faces=ids(slots);
        calls=calls+1;faceInputs=faceInputs+numel(faces);
        assert(calls<=2*max(1,ni)+1&&faceInputs<=16*max(1,ni), ...
            'InterfacePreview:WorkBudget','原内核界面检查超过显式工作预算，不漏报约束。');
        try
            cutCellExactEulerFlux(candidate(f.cellMinus(faces),:),candidate(f.cellPlus(faces),:),f.normal(faces,:),gas);
        catch e
            if ~strcmp(e.identifier,'CutExactEuler:Temperature'),rethrow(e);end
            failedBatchCount=failedBatchCount+1;
            % 面索引和温度来自原内核的真实拒绝文本，不放宽温域来读取非法采样。
            token=regexp(e.message,'(\d+)，T=([^ ]+) K','tokens','once');
            assert(~isempty(token),'InterfacePreview:Diagnostic','原内核温域拒绝缺少明确面号/温度。');
            position=str2double(token{1});temperature=str2double(token{2});
            assert(isfinite(position)&&position>=1&&position<=numel(slots)&&position==fix(position)&& ...
                isfinite(temperature)&&(temperature<gas.temperatureRange(1)||temperature>gas.temperatureRange(2)), ...
                'InterfacePreview:Diagnostic','原内核拒绝的面号或温度无效。');
            at=slots(position);if ~failed(at),failedT(at)=temperature;end;failed(at)=true;
            if numel(slots)>1
                middle=floor((span(1)+span(2))/2);queue{end+1}=[middle+1 span(2)];queue{end+1}=[span(1) middle]; %#ok<AGROW>
            end
        end
    end
end
report=struct('method','original-temperature-domain-batched-Riemann-preview-with-explicit-blocked-cell-pairs', ...
    'completeInteriorFaceCount',ni,'assessedFaceIDs',ids(eligible),'blockedFaceIDs',blocked, ...
    'blockedByOutOfRangeCellIDs',find(~valid),'allInteriorFacesAssessed',isempty(blocked), ...
    'allInputFacesAccountedFor',numel(eligibleSlots)+numel(blocked)==ni, ...
    'predictedViolatingInterfaceFaceIDs',ids(failed),'predictedInterfaceTemperatureK',failedT(failed), ...
    'maximumBatchFaces',options.maximumBatchFaces,'additionalArrayBoundBytes',arrayBound, ...
    'originalKernelCalls',calls,'kernelFaceInputsIncludingRechecks',faceInputs,'temperatureRejectedBatches',failedBatchCount, ...
    'originalTemperatureRangeK',gas.temperatureRange,'temperatureDomainRelaxed',false, ...
    'nonTemperatureFailuresSuppressed',false,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
end
