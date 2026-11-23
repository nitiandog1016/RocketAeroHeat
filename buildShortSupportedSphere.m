function g=buildShortSupportedSphere(extensionD)
% 模块一：由 Short 原始离散球面生成带两根赤道支管的试件，尺寸使用 SI。
% extensionD 是截断站到理想球极切平面的距离，以球径为单位。
if nargin<1,extensionD=1;end
validateattributes(extensionD,{'numeric'},{'scalar','real','finite','positive'});
root=fileparts(mfilename('fullpath'));folder=fullfile(root,'cases','short_sphere');
source=jsondecode(fileread(fullfile(folder,'test65.json')));
key=jsondecode(fileread(fullfile(folder,'surface65.json')));
archive=fullfile(folder,'surface65.zip');
assert(strcmp(cm011FileSHA256(archive),key.archiveSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(folder,'test65.json')),key.caseSHA256));
temp=tempname;mkdir(temp);cleanup=onCleanup(@()clearSource(temp)); %#ok<NASGU>
unzip(archive,temp);path=fullfile(temp,'body.json');bodySHA=cm011FileSHA256(path);
assert(strcmp(bodySHA,'05738f3c56c74b8419468b777089da3555c202ca5f2624063a4c45c09008af68'));
body=jsondecode(fileread(path));D=source.diameterIn*.0254;V=body.V/D;F=body.F;
n=size(F,1);nSide=64;radius=source.supportDiameterIn/(2*source.diameterIn);
phase=pi/nSide;phi=phase+(0:nSide-1)'*2*pi/nSide;
corner=radius*[cos(phi) sin(phi)];foot=polyshape(corner);
normal=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
original=vecnorm(normal,2,2)/2;
covered=zeros(n,1);allV=V;allF=F;parent=(1:n)';keep=true(n,1);
% 模块二：在每个受影响父三角形的平面内精确裁去正多边形柱的投影。
for i=1:n
    p=V(F(i,:),:);xz=p(:,[1 3]);
    if any(min(xz)>radius) || any(max(xz)<-radius),continue;end
    tri=polyshape(xz);shadow=intersect(tri,foot);
    if area(shadow)==0,continue;end
    assert(abs(normal(i,2))>0.85*norm(normal(i,:)));
    covered(i)=area(shadow)*norm(normal(i,:))/abs(normal(i,2));keep(i)=false;
    rest=subtract(tri,foot);
    if area(rest)==0,continue;end
    parts=triangulation(rest);uv=parts.Points;
    y=p(1,2)-(uv-xz(1,:))*normal(i,[1 3])'/normal(i,2);
    v=[uv(:,1) y uv(:,2)];f=parts.ConnectivityList;
    orientation=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2)*normal(i,:)';
    flip=orientation<0;f(flip,[2 3])=f(flip,[3 2]);
    allF=[allF;size(allV,1)+f];allV=[allV;v]; %#ok<AGROW>
    parent=[parent;repmat(i,size(f,1),1)]; %#ok<AGROW>
end
keep=[keep;true(size(allF,1)-n,1)];sphereF=allF(keep,:);parent=parent(keep);
% 仅合并相同几何交点；保留原坐标，不把点投影到平滑球上。
[~,first,index]=unique(round(allV(sphereF(:),:),12),'rows');
v=allV(sphereF(first),:);f=reshape(index,size(sphereF));part=ones(size(f,1),1);
edges=freeBoundary(triangulation(f,v));
assert(all(sign(v(edges(:,1),2))==sign(v(edges(:,2),2))));
% 模块三：两个真实接头共享球面的自由边。第一排管面消化接缝中额外的交点。
joins=cell(2,1);endY=.5+extensionD;
for side=1:2
    direction=3-2*side;
    near=edges(sign(v(edges(:,1),2))==direction,:);
    node=unique(near(:));angles=mod(atan2(v(node,3),v(node,1))-phase,2*pi);
    angles(abs(angles-2*pi)<1e-10)=0;
    [~,order]=sort(angles);node=node(order);
    loop=[node circshift(node,-1)];
    assert(isequal(sortrows(sort(loop,2)),sortrows(sort(near,2))));
    joins{side}=loop;
    [v,newF,newPart]=tube(v,node,corner,endY,direction,side);
    f=[f;newF];part=[part;newPart];parent=[parent;zeros(size(newF,1),1)]; %#ok<AGROW>
end
% 模块四：验收闭合、相邻面方向、体积与逐父面的面积，七站按原测线落在外露面上。
v=v*D;vector=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
areas=vecnorm(vector,2,2)/2;assert(all(areas>0));
edge=[f(:,[1 2]);f(:,[2 3]);f(:,[3 1])];
[pairs,~,eid]=unique(sort(edge,2),'rows');
assert(all(accumarray(eid,1)==2) && all(accumarray(eid,sign(edge(:,2)-edge(:,1)))==0));
assert(size(v,1)-size(pairs,1)+size(f,1)==2 && numel(unique(conncomp(graph(pairs(:,1),pairs(:,2)))))==1);
assert(sum(sum(vector.*v(f(:,1),:),2))/6>0);
sphere=part==1;
exposed=accumarray(parent(sphere),areas(sphere),[n 1]);original=original*D^2;covered=covered*D^2;
assert(max(abs(original-covered-exposed))<1e-16);
g=struct('V',v,'F',f,'regionID',part,'parentID',parent,'areaM2',areas,'joinPositive',joins{1},'joinNegative',joins{2});
g.parents=table((1:n)',original,exposed,covered,'VariableNames',{'parentID','originalM2','exposedM2','coveredM2'});
g.stations=locateStations(v,f,sphere,parent,source.angleDeg,D);
% 先成对累加各分区的小面积，避免几万个三角形因加入顺序造成总账尾数漂移。
regionAreas=arrayfun(@(k)pairSum(areas(part==k)),(1:5)');
g.summary=struct('test',source.test,'bodySHA256',bodySHA,'caseSHA256',key.caseSHA256,'diameterM',D, ...
    'tubeDiameterM',2*radius*D,'extensionD',extensionD,'endYM',endY*D,'polygonSides',nSide, ...
    'regionNames',{{'exposed_silver','positive_tube','negative_tube','positive_cut','negative_cut'}}, ...
    'regionAreaM2',regionAreas,'regionFaceCount',accumarray(part,1,[5 1]), ...
    'smoothGridM2',pairSum(original),'coveredM2',pairSum(covered),'analyticCoveredM2',pi*D^2*(1-sqrt(1-(2*radius)^2)), ...
    'publishedAreaFactor',source.supportAreaFactor,'meshAreaFactor',pairSum(original)/regionAreas(1), ...
    'solverWallAreaM2',sum(regionAreas(1:3)),'artificialCutAreaM2',sum(regionAreas(4:5)), ...
    'pressureHeatSolvedM2',0,'lengthStatus','nominal','supportTemperatureStatus','unknown');
end

function total=pairSum(values)
while numel(values)>1
    if mod(numel(values),2),values(end+1)=0;end
    values=values(1:2:end)+values(2:2:end);
end
total=values;
end

function [v,f,type]=tube(v,seam,corner,endY,sgn,side)
% 模块五：支管侧壁沿 ±Y 延伸，每个断面用同一规则多边形；远端圆盘只作几何截断。
n=size(corner,1);spacing=norm(corner(1,:)-corner(2,:));
% 以实际接缝外侧作为第一排的起点，保证首排的完整可见长度。
joinY=max(sgn*v(seam,2));layers=max(1,ceil((endY-joinY)/spacing));
levels=joinY+(1:layers)'*(endY-joinY)/layers;
ring=reshape(size(v,1)+(1:n*layers),n,layers);
v=[v;repmat([corner(:,1) zeros(n,1) corner(:,2)],layers,1)];
v(ring(:),2)=sgn*repelem(levels,n);
anchor=zeros(n,1);
for k=1:n
    [gap,anchor(k)]=min(vecnorm(v(seam,[1 3])-corner(k,:),2,2));assert(gap<1e-10);
end
assert(anchor(1)==1 && all(diff(anchor)>0));
first=zeros(numel(seam)+n,3);used=0;
for k=1:n
    next=mod(k,n)+1;last=anchor(next);
    if next==1,last=numel(seam)+1;end
    chain=seam(mod((anchor(k):last)-1,numel(seam))+1);count=numel(chain)-1;
    first(used+(1:count),:)=[chain(1:end-1) chain(2:end) repmat(ring(k,1),count,1)];
    first(used+count+1,:)=[chain(end) ring(next,1) ring(k,1)];used=used+count+1;
end
a=ring(:,1:end-1);b=circshift(a,-1,1);c=circshift(ring(:,2:end),-1,1);d=ring(:,2:end);
f=[first;a(:) b(:) c(:);a(:) c(:) d(:)];sideCount=size(f,1);
v(end+1,:)=[0 sgn*endY 0];cap=ring(:,end);
f=[f;cap circshift(cap,-1) repmat(size(v,1),n,1)];
normals=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
radial=v(f(1:sideCount,1),[1 3]);
flip=[sum(normals(1:sideCount,[1 3]).*radial,2)<0;sgn*normals(sideCount+1:end,2)<0];
f(flip,[2 3])=f(flip,[3 2]);type=[repmat(1+side,sideCount,1);repmat(3+side,n,1)];
end

function t=locateStations(v,f,sphere,parent,angles,D)
% 模块六：用测站射线与原球面子三角形求交。共边或共点时选编号最小的命中面。
ids=find(sphere);a=v(f(ids,1),:);e=v(f(ids,2),:)-a;h=v(f(ids,3),:)-a;
normal=cross(e,h,2);nn=sum(normal.^2,2);numerator=sum(normal.*a,2);
rows=zeros(numel(angles),10);
for k=1:numel(angles)
    direction=[cosd(angles(k)) 0 sind(angles(k))];
    den=normal*direction';ok=den>1e-10*sqrt(nn);distance=zeros(size(den));distance(ok)=numerator(ok)./den(ok);
    delta=distance.*direction-a;
    u=sum(cross(delta,h,2).*normal,2)./nn;w=sum(cross(e,delta,2).*normal,2)./nn;
    weights=[1-u-w u w];hit=find(ok & distance>0 & all(weights>=-1e-10,2));
    assert(~isempty(hit));i=hit(1);p=distance(i)*direction;
    assert(norm(p-weights(i,:)*v(f(ids(i),:),:))<1e-12 && abs(distance(i)-D/2)<1e-5);
    rows(k,:)=[ids(i) parent(ids(i)) p weights(i,:) numel(hit) D/2-distance(i)];
end
t=array2table([angles(:) rows],'VariableNames',{'angleDeg','faceID','parentID','xM','yM','zM','w1','w2','w3','incidentFaces','radialChordOffsetM'});
end

function clearSource(folder)
% 清理本次临时展开的四个源文件。
for name={'body.json','stations.csv','wall.csv','summary.json'}
    p=fullfile(folder,name{1});if isfile(p),delete(p);end
end
rmdir(folder);
end
