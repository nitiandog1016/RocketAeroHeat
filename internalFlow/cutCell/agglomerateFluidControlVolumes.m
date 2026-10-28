function [group,report]=agglomerateFluidControlVolumes(source,options)
%AGGLOMERATEFLUIDCONTROLVOLUMES 相邻真实流体控制体的有界守恒聚合原型。
% 中文模块说明：只沿正开放共享面连接，保存全部流体体积及一阶体积矩、
% 全部原墙片；组内共享面严格消失，组外仍是一面一条共享通量记录。
% 不删除tiny体积、不移动质心、不改变墙面积/法向/多边形，不推进NS。
% 对无法改善的组显式报告拒绝；聚合几何通过不等于时间稳定证明。
if nargin<2,options=struct;end
defaults=struct('minimumVolumeFraction',.1,'maximumMembers',12,'maximumIndexExtent',[4 4 4], ...
    'maximumPasses',20,'failOnUnresolved',false);
names=fieldnames(defaults);for k=1:numel(names),if ~isfield(options,names{k}),options.(names{k})=defaults.(names{k});end,end
assert(options.minimumVolumeFraction>0 && options.minimumVolumeFraction<1);
assert(options.maximumMembers>=1 && options.maximumMembers==fix(options.maximumMembers));
assert(numel(options.maximumIndexExtent)==3 && all(options.maximumIndexExtent>=1));
v=source.cells.volume(:); c=source.cells.centroid; N=numel(v); h=source.grid.spacing;
if isscalar(h),h=repmat(h,1,3);end
baseVolume=prod(h); targetVolume=options.minimumVolumeFraction*baseVolume;
fluid=v>0; assert(any(fluid),'CutMerge:NoFluid','输入须含流体控制体。');
faces=source.openFaces; walls=source.walls;
assert(isfield(walls,'normal'),'CutMerge:MissingNormal','须使用冻结原面可靠单位法向。');
assert(max(abs(sqrt(sum(walls.normal.^2,2))-1))<=16*eps);
scale=max(1,max(abs(source.inputGeometry.V),[],'all'));
distanceRoundoff=64*eps(scale);
active=fluid; map=zeros(N,1); map(fluid)=find(fluid);
members=cell(N,1); wallIDs=cell(N,1); volumes=zeros(N,1); centers=nan(N,3);
for id=find(fluid)',members{id}=id;volumes(id)=v(id);centers(id,:)=c(id,:);end
for w=1:numel(walls.area),id=walls.ownerCell(w);if id>0,wallIDs{id}(end+1)=w;end,end
goodFace=faces.area>0 & faces.cellMinus>0 & faces.cellPlus>0;
am=faces.cellMinus(goodFace);ap=faces.cellPlus(goodFace);aa=faces.area(goodFace);
adj=sparse([am;ap],[ap;am],[aa;aa],N,N);
metric=cell(N,1);for id=find(fluid)',metric{id}=wallMetric(centers(id,:),wallIDs{id},walls,distanceRoundoff);end
before=quality(source,v,c,wallIDs,fluid,baseVolume,targetVolume,distanceRoundoff);
history=struct('leftRepresentative',{},'rightRepresentative',{},'members',{},'volume',{},'centroid',{}, ...
    'oldNonpositiveWalls',{},'newNonpositiveWalls',{},'oldMinimumDistance',{},'newMinimumDistance',{},'sharedAperture',{});
mergeCount=0; totalCandidates=0; rejectedByBounds=0; rejectedByQuality=0;twoNeighbourMerges=0;

%% 贪心只接受局部壁距改善或不破坏壁距的小体积改善；每次合并都保持图连接
for pass=1:options.maximumPasses
    ids=find(active); bad=zeros(numel(ids),1);
    for k=1:numel(ids),bad(k)=metric{ids(k)}.nonpositiveCount;end
    wanted=bad>0 | volumes(ids)<targetVolume; ids=ids(wanted); bad=bad(wanted);
    if isempty(ids),break;end
    [~,order]=sortrows([volumes(ids)>=targetVolume,-bad,volumes(ids)],[1 2 3]);ids=ids(order);
    passMerges=0;
    for g=ids'
        if ~active(g) || (metric{g}.nonpositiveCount==0 && volumes(g)>=targetVolume),continue;end
        neighbourOriginal=find(any(adj(members{g},:),1)); neighbours=unique(map(neighbourOriginal));
        neighbours=neighbours(neighbours>0 & neighbours~=g); best=[];bestScore=[];
        for other=neighbours(:)'
            if ~active(other),continue;end
            totalCandidates=totalCandidates+1; joined=[members{g} members{other}];
            extent=max(source.cells.index(joined,:),[],1)-min(source.cells.index(joined,:),[],1)+1;
            if numel(joined)>options.maximumMembers || any(extent>options.maximumIndexExtent)
                rejectedByBounds=rejectedByBounds+1;continue;
            end
            [vv,cc]=conservedCenter(v,c,joined); ww=[wallIDs{g} wallIDs{other}];
            mm=wallMetric(cc,ww,walls,distanceRoundoff);
            oldBad=metric{g}.nonpositiveCount+metric{other}.nonpositiveCount;
            oldMin=min(metric{g}.minimumDistance,metric{other}.minimumDistance);
            if oldBad>0
                improves=mm.nonpositiveCount<oldBad || ...
                    (mm.nonpositiveCount<=oldBad && mm.minimumDistance>oldMin+distanceRoundoff);
            else
                improves=mm.nonpositiveCount==0 && min(vv,targetVolume)>min(volumes(g),targetVolume)+eps(targetVolume);
            end
            if ~improves,rejectedByQuality=rejectedByQuality+1;continue;end
            shared=full(sum(sum(adj(members{g},members{other}))));
            minimum=mm.minimumDistance;if isinf(minimum),minimum=10*max(h);end
            score=[mm.nonpositiveCount==0 && vv>=targetVolume,oldBad-mm.nonpositiveCount, ...
                minimum/max(h),shared/max(h)^2,-numel(joined)];
            if isempty(best) || lexGreater(score,bestScore)
                best=struct('others',other,'members',joined,'volume',vv,'centroid',cc,'wallIDs',ww, ...
                    'metric',mm,'oldBad',oldBad,'oldMin',oldMin,'shared',shared);bestScore=score;
            end
        end
        % 单邻居无法改善时，只前瞻两条真实相邻通路的联合组；不接受坏中间组。
        % 这处理尾缘等相互制约墙面需同时向轴向/径向移动守恒质心的局部情形。
        if isempty(best) && numel(neighbours)>=2
            for left=1:numel(neighbours)-1
                for right=left+1:numel(neighbours)
                    others=[neighbours(left) neighbours(right)];
                    if ~all(active(others)),continue;end
                    totalCandidates=totalCandidates+1;
                    joined=[members{g} members{others(1)} members{others(2)}];
                    extent=max(source.cells.index(joined,:),[],1)-min(source.cells.index(joined,:),[],1)+1;
                    if numel(joined)>options.maximumMembers || any(extent>options.maximumIndexExtent)
                        rejectedByBounds=rejectedByBounds+1;continue;
                    end
                    [vv,cc]=conservedCenter(v,c,joined);ww=[wallIDs{g} wallIDs{others(1)} wallIDs{others(2)}];
                    mm=wallMetric(cc,ww,walls,distanceRoundoff);
                    oldBad=metric{g}.nonpositiveCount+metric{others(1)}.nonpositiveCount+metric{others(2)}.nonpositiveCount;
                    oldMin=min([metric{g}.minimumDistance metric{others(1)}.minimumDistance metric{others(2)}.minimumDistance]);
                    if oldBad>0
                        improves=mm.nonpositiveCount<oldBad || ...
                            (mm.nonpositiveCount<=oldBad && mm.minimumDistance>oldMin+distanceRoundoff);
                    else
                        improves=mm.nonpositiveCount==0 && min(vv,targetVolume)>min(volumes(g),targetVolume)+eps(targetVolume);
                    end
                    if ~improves,rejectedByQuality=rejectedByQuality+1;continue;end
                    shared=full(sum(sum(adj(members{g},[members{others(1)} members{others(2)}]))));
                    minimum=mm.minimumDistance;if isinf(minimum),minimum=10*max(h);end
                    score=[mm.nonpositiveCount==0 && vv>=targetVolume,oldBad-mm.nonpositiveCount, ...
                        minimum/max(h),shared/max(h)^2,-numel(joined)];
                    if isempty(best) || lexGreater(score,bestScore)
                        best=struct('others',others,'members',joined,'volume',vv,'centroid',cc,'wallIDs',ww, ...
                            'metric',mm,'oldBad',oldBad,'oldMin',oldMin,'shared',shared);bestScore=score;
                    end
                end
            end
        end
        if isempty(best),continue;end
        selected=[g best.others];rep=min(selected);gone=selected(selected~=rep);active(gone)=false;
        members{rep}=best.members;wallIDs{rep}=best.wallIDs;volumes(rep)=best.volume;centers(rep,:)=best.centroid;
        metric{rep}=best.metric;map(best.members)=rep;
        for goneID=gone,members{goneID}=[];wallIDs{goneID}=[];end
        mergeCount=mergeCount+1;passMerges=passMerges+1;
        if numel(best.others)==2,twoNeighbourMerges=twoNeighbourMerges+1;end
        history(mergeCount)=struct('leftRepresentative',g,'rightRepresentative',best.others,'members',best.members, ...
            'volume',best.volume,'centroid',best.centroid,'oldNonpositiveWalls',best.oldBad, ...
            'newNonpositiveWalls',best.metric.nonpositiveCount,'oldMinimumDistance',best.oldMin, ...
            'newMinimumDistance',best.metric.minimumDistance,'sharedAperture',best.shared);
    end
    fprintf('  aggregation pass %d: %d merges, %d groups\n',pass,passMerges,sum(active));
    if passMerges==0,break;end
end

%% 压紧组编号，原共享面逐条映射；内部面没有任何外露通量副本
roots=find(active); G=numel(roots);compact=zeros(N,1);compact(roots)=1:G;
originalCellToGroup=zeros(N,1);originalCellToGroup(fluid)=compact(map(fluid));
groupMembers=members(roots);gv=volumes(roots);gc=centers(roots,:);
fm=zeros(numel(faces.area),1);fp=fm;
positiveMinus=faces.cellMinus>0;positivePlus=faces.cellPlus>0;
fm(positiveMinus)=originalCellToGroup(faces.cellMinus(positiveMinus));
fp(positivePlus)=originalCellToGroup(faces.cellPlus(positivePlus));
removedInternal=faces.area>0 & fm>0 & fm==fp;
keep=faces.area>0 & ~removedInternal;
assert(all((fm(keep)>0 | fp(keep)>0) & fm(keep)~=fp(keep)));
originalFaceToGroupFace=zeros(numel(faces.area),1);originalFaceToGroupFace(keep)=1:sum(keep);
gf=selectRows(faces,keep);gf.originalCartesianFaceID=find(keep);
gf.originalCellMinus=faces.cellMinus(keep);gf.originalCellPlus=faces.cellPlus(keep);
gf.cellMinus=fm(keep);gf.cellPlus=fp(keep);gf.lowerCellID=gf.cellMinus;gf.upperCellID=gf.cellPlus;
gw=walls;gw.originalOwnerCell=walls.ownerCell;activeWall=walls.ownerCell>0;
gw.ownerCell=zeros(size(walls.ownerCell));gw.ownerCell(activeWall)=originalCellToGroup(walls.ownerCell(activeWall));
gw.fluidOwnerCellID=gw.ownerCell;gw.cellID=gw.ownerCell;
gw.ownershipPolicy='conservative connected agglomerate containing original fluid owner; originalOwnerCell retained';
gw.polygonVerticesGroupLocal=cell(size(walls.polygonVertices));
for w=1:numel(walls.area)
    id=gw.ownerCell(w);if id>0,gw.polygonVerticesGroupLocal{w}=walls.polygonVerticesGridLocal{w}-(gc(id,:)-source.grid.origin);end
end
state=ones(G,1);for id=1:G,if numel(groupMembers{id})==1,state(id)=source.cells.state(groupMembers{id});end,end
group=struct('grid',source.grid,'cells',struct('volume',gv,'centroid',gc,'center',gc,'state',state, ...
    'index',source.cells.index(roots,:),'members',{groupMembers},'memberCount',cellfun(@numel,groupMembers)), ...
    'openFaces',gf,'walls',gw,'inputGeometry',source.inputGeometry,'originalCellToGroup',originalCellToGroup, ...
    'originalFaceToGroupFace',originalFaceToGroupFace,'removedInternalOriginalFaces',removedInternal, ...
    'closedOriginalFaces',faces.area==0,'members',{groupMembers});
[closure,boundaryArea]=groupClosure(group);group.closure=closure;
group.closureAreaVector=closure.areaVectorResidual;group.closureTensor=closure.firstMomentResidual;
group.cells.boundaryArea=boundaryArea;group.cells.volumeOverBoundaryArea=gv./boundaryArea;

%% 守恒、原墙身份、连通和有限扩展：没有以几何重构篡改源控制体
firstBefore=sum(v(fluid).*c(fluid,:),1);firstAfter=sum(gv.*gc,1);
volumeError=abs(sum(gv)-sum(v));firstError=max(abs(firstAfter-firstBefore));
assert(volumeError<=128*eps(sum(v)),'CutMerge:VolumeConservation','体积守恒失败。');
assert(firstError<=512*eps(max(1,max(abs(firstBefore)))),'CutMerge:FirstMomentConservation','体积一阶矩守恒失败。');
assert(isequaln(gw.area,walls.area) && isequaln(gw.centroid,walls.centroid) && isequaln(gw.normal,walls.normal));
assert(isequaln(gw.areaVector,walls.areaVector) && isequaln(gw.originalFaceID,walls.originalFaceID));
assert(isequaln(gw.polygonVertices,walls.polygonVertices),'CutMerge:WallChanged','原曲面多边形不可改变。');
connected=true(G,1);extents=zeros(G,3);diameters=zeros(G,1);
for id=1:G
    mem=groupMembers{id};connected(id)=isConnected(adj,mem);
    extents(id,:)=max(source.cells.index(mem,:),[],1)-min(source.cells.index(mem,:),[],1)+1;
    diameters(id)=norm(extents(id,:).*h);
end
assert(all(connected) && max(cellfun(@numel,groupMembers))<=options.maximumMembers);
assert(all(extents<=options.maximumIndexExtent,'all'));
assert(closure.maxRelativeFluidTensorError<1e-8 && closure.normalizedAreaError<1e-10, ...
    'CutMerge:Closure','聚合逐组九矩/面积向量闭合失败。');
afterWallIDs=wallIDs(roots);after=quality(group,gv,gc,afterWallIDs,true(G,1),baseVolume,targetVolume,distanceRoundoff);
unresolved=[];unresolvedReasons=cell(G,1);
for id=1:G
    if after.nonpositiveWallCountByCell(id)>0 || gv(id)<targetVolume
        unresolved(end+1,1)=id; %#ok<AGROW>
        if after.nonpositiveWallCountByCell(id)>0
            unresolvedReasons{id}='rejected: no accepted adjacent merge improves all-wall criterion within fixed member/extent limits';
        else,unresolvedReasons{id}='rejected: no accepted adjacent merge reaches volume target within fixed member/extent limits';end
    end
end
group.unresolvedGroupIDs=unresolved;group.unresolvedReasons=unresolvedReasons;
group.mergeHistory=history;group.options=options;group.physicalSolution=false;group.eligibleForSurfaceCoverage=false;
group.algorithm='bounded connected conservative fluid-control-volume agglomeration across original positive shared apertures';
group.transportStabilityProven=false;group.allRequestedCellsResolved=isempty(unresolved);
group.cells.indexExtent=extents;group.cells.boundingBoxDiameter=diameters;
report=struct('before',before,'after',after,'mergeCount',mergeCount,'groups',G,'maximumMembers',max(cellfun(@numel,groupMembers)), ...
    'maximumIndexExtent',max(extents,[],1),'maximumBoundingBoxDiameter',max(diameters),'allGroupsConnected',all(connected), ...
    'volumeConservationError',volumeError,'firstVolumeMomentConservationError',firstError, ...
    'removedInternalFaces',sum(removedInternal),'remainingUniqueSharedFaces',sum(keep), ...
    'unresolvedGroups',numel(unresolved),'allRequestedCellsResolved',isempty(unresolved), ...
    'distanceRoundoff',distanceRoundoff,'totalCandidateEvaluations',totalCandidates, ...
    'candidatesRejectedByBounds',rejectedByBounds,'candidatesRejectedByQuality',rejectedByQuality,'twoNeighbourMerges',twoNeighbourMerges, ...
    'closure',rmfield(closure,{'areaVectorResidual','firstMomentResidual','relativeTensorErrorPerFluidVolume'}), ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'transportStabilityProven',false);
group.report=report;
if options.failOnUnresolved && ~isempty(unresolved)
    error('CutMerge:NoFeasibleNeighbour','%d个问题组无可接受相邻聚合；保留原体积、质心及墙，不宣称已修复。',numel(unresolved));
end
end

function m=wallMetric(center,ids,walls,tol)
if isempty(ids),m=struct('nonpositiveCount',0,'minimumDistance',Inf);return;end
d=sum((center-walls.centroid(ids,:)).*(-walls.normal(ids,:)),2);
m=struct('nonpositiveCount',sum(d<=tol),'minimumDistance',min(d));
end

function yes=lexGreater(a,b)
k=find(a~=b,1);yes=~isempty(k) && a(k)>b(k);
end

function [volume,centroid]=conservedCenter(v,c,ids)
volume=sum(v(ids));reference=c(ids(1),:);
centroid=reference+sum(v(ids).*(c(ids,:)-reference),1)/volume;
end

function out=selectRows(in,keep)
out=in;names=fieldnames(in);N=numel(keep);
for k=1:numel(names)
    x=in.(names{k});if (isnumeric(x)||islogical(x)) && size(x,1)==N,out.(names{k})=x(keep,:);end
end
end

function [cl,boundaryArea]=groupClosure(g)
N=numel(g.cells.volume);av=zeros(N,3);mt=zeros(N,3,3);boundaryArea=zeros(N,1);
f=g.openFaces;
for k=1:numel(f.area)
    an=f.area(k)*f.normal(k,:);
    for side=1:2
        if side==1,id=f.cellMinus(k);sgn=1;else,id=f.cellPlus(k);sgn=-1;end
        if id==0,continue;end
        r=f.centroid(k,:)-g.cells.centroid(id,:);av(id,:)=av(id,:)+sgn*an;
        mt(id,:,:)=mt(id,:,:)+reshape(r'*(sgn*an),1,3,3);boundaryArea(id)=boundaryArea(id)+f.area(k);
    end
end
w=g.walls;
for k=1:numel(w.area)
    id=w.ownerCell(k);if id==0,continue;end
    an=w.areaVector(k,:);r=w.centroid(k,:)-g.cells.centroid(id,:);
    av(id,:)=av(id,:)+an;mt(id,:,:)=mt(id,:,:)+reshape(r'*an,1,3,3);boundaryArea(id)=boundaryArea(id)+w.area(k);
end
for a=1:3,mt(:,a,a)=mt(:,a,a)-g.cells.volume;end
relative=max(reshape(abs(mt),N,9),[],2)./g.cells.volume;
cl=struct('areaVectorResidual',av,'firstMomentResidual',mt,'relativeTensorErrorPerFluidVolume',relative, ...
    'maxRelativeFluidTensorError',max(relative),'normalizedAreaError',max(abs(av),[],'all')/max(g.grid.spacing)^2, ...
    'normalizedTensorError',max(abs(mt),[],'all')/prod(g.grid.spacing));
end

function yes=isConnected(adj,mem)
seen=false(1,numel(mem));seen(1)=true;
while true
    next=any(adj(mem(seen),mem),1);new=logical(next)&~seen;
    if ~any(new),break;end
    seen=seen|new;
end
yes=all(seen);
end

function out=quality(g,v,c,wallIDs,fluid,baseVolume,target,tol)
N=numel(v);bad=zeros(N,1);minimum=inf(N,1);
for id=find(fluid)'
    m=wallMetric(c(id,:),wallIDs{id},g.walls,tol);bad(id)=m.nonpositiveCount;minimum(id)=m.minimumDistance;
end
[inside,ambiguous]=centroidSolidDiagnostics(c,fluid,g.inputGeometry);
boundary=zeros(N,1);f=g.openFaces;
for k=1:numel(f.area)
    if f.area(k)==0,continue;end
    if f.cellMinus(k)>0,boundary(f.cellMinus(k))=boundary(f.cellMinus(k))+f.area(k);end
    if f.cellPlus(k)>0,boundary(f.cellPlus(k))=boundary(f.cellPlus(k))+f.area(k);end
end
for k=1:numel(g.walls.area),id=g.walls.ownerCell(k);if id>0,boundary(id)=boundary(id)+g.walls.area(k);end,end
acoustic=v(fluid)./boundary(fluid);
out=struct('fluidControlVolumes',sum(fluid),'minimumPositiveVolume',min(v(fluid)), ...
    'minimumVolumeFraction',min(v(fluid))/baseVolume,'belowVolumeTarget',sum(v(fluid)<target), ...
    'nonpositiveWallPatches',sum(bad),'nonpositiveWallCountByCell',bad,'minimumWallSignedDistance',min(minimum), ...
    'centroidsInsideSolid',sum(inside),'centroidInsideSolidIDs',find(inside),'ambiguousCentroidWindingIDs',find(ambiguous), ...
    'minimumAcousticVolumeOverArea',min(acoustic),'maximumAcousticVolumeOverArea',max(acoustic));
end

function [inside,ambiguous]=centroidSolidDiagnostics(c,fluid,geometry)
% 中文模块说明：立体角仅作质心位置诊断，不改变任何组几何。
V=geometry.V;F=geometry.F;lo=min(V,[],1);hi=max(V,[],1);inside=false(numel(fluid),1);ambiguous=inside;
candidate=find(fluid & all(c>=lo & c<=hi,2));
for id=candidate'
    a=V(F(:,1),:)-c(id,:);b=V(F(:,2),:)-c(id,:);d=V(F(:,3),:)-c(id,:);
    la=sqrt(sum(a.^2,2));lb=sqrt(sum(b.^2,2));ld=sqrt(sum(d.^2,2));
    numerator=sum(a.*cross(b,d,2),2);
    denominator=la.*lb.*ld+sum(a.*b,2).*ld+sum(b.*d,2).*la+sum(d.*a,2).*lb;
    winding=sum(2*atan2(numerator,denominator))/(4*pi);
    if abs(winding-1)<1e-8,inside(id)=true;elseif abs(winding)>1e-8,ambiguous(id)=true;end
end
end
