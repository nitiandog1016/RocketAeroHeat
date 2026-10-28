function report=inspectRocketTemperatureCapacity(folder,stateFile)
%INSPECTROCKETTEMPERATURECAPACITY 读取当前实际端点，完整构造超旧行限的温域导数用于容量评估。
% 模块一：绑定原状态/几何SHA，仅使用该端点已保存的完整NS采样，不更新流场。
% 模块二：按旧构造器的每批16面分批计算，全部面都保留，并核对当前采样温度。
% 模块三：逐项相同的H行共享表示，报告真实支撑/条件数；不据此宣称新NS或区域合格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
stateFile=char(java.io.File(stateFile).getCanonicalPath());
geometryFile=fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat');
sourceStateSHA256=rocketStreamingFileSHA256(stateFile);geometrySHA256=rocketStreamingFileSHA256(geometryFile);
s=load(stateFile,'result','gas');d=load(geometryFile,'domain');Q=s.result.Q;gas=s.gas;scale=s.result.residualScale;
audit=s.result.info.exactEulerFluxAudit;f=d.domain.openFaces;inside=find(f.area>0&f.cellMinus>0&f.cellPlus>0);
assert(isequal(inside,audit.faceID(:)));T=audit.samplePrimitive(:,5);upper=gas.temperatureRange(2)-T<=10;lower=T-gas.temperatureRange(1)<=10;
slots=find(upper|lower);assert(all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)));ids=inside(slots);assert(numel(ids)>16);
fields={'faceIDs','ownerCV','neighbourCV','normal','QL','QR','H','leftRowsK','rightRowsK', ...
    'baseSamplePrimitive','baseTemperatureK','remainingUpperTemperatureK','sampleRegion','rowTwoNormK'};
rows=[];
for first=1:16:numel(ids)
    part=buildRocketInterfaceTemperatureRows(Q,d.domain,gas,scale,ids(first:min(first+15,end)));
    if isempty(rows),rows=part;else
        for key=fields,rows.(key{1})=[rows.(key{1});part.(key{1})];end
    end
end
rows.temperatureBoundSigns=ones(numel(ids),1);rows.temperatureBoundSigns(lower(slots))=-1;
assert(isequal(rows.faceIDs,ids)&&max(abs(rows.baseTemperatureK-T(slots)))<=1e-10);
Ht=rows.H';representatives=[];mapping=zeros(numel(ids),1);
for j=1:numel(ids)
    found=0;for k=1:numel(representatives),if isequal(Ht(:,j),Ht(:,representatives(k))),found=k;break;end,end
    if found==0,representatives(end+1,1)=j;found=numel(representatives);end %#ok<AGROW>
    mapping(j)=found;
end
assert(isequal(rows.H,rows.H(representatives(mapping),:)));
support=find(any(Ht,2));compact=full(Ht(support,representatives));singularValues=svd(compact./vecnorm(compact,2,1));
report=struct('status','observed-current-over-limit-all-interface-rows-with-exact-representation', ...
    'sourceStateFile',stateFile,'sourceStateSHA256',sourceStateSHA256,'geometrySHA256',geometrySHA256, ...
    'originalFaceCount',numel(ids),'upperFaceCount',sum(upper(slots)),'lowerFaceCount',sum(lower(slots)), ...
    'exactDistinctRowCount',numel(representatives),'originalToExactRow',mapping,'representativeRows',representatives, ...
    'compactPhysicalUnknowns',numel(support),'normalizedCompactSingularValues',singularValues, ...
    'currentCompleteInteriorFaceCount',numel(inside),'allOriginalRowsExactlyRepresented',true, ...
    'QRCoordinateRankNotYetClaimed',true,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
assert(strcmp(sourceStateSHA256,rocketStreamingFileSHA256(stateFile))&&strcmp(geometrySHA256,rocketStreamingFileSHA256(geometryFile)));
save(fullfile(folder,'temperature_capacity.mat'),'rows','report','representatives','mapping','compact','support','-v7');
fid=fopen(fullfile(folder,'temperature_capacity.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
