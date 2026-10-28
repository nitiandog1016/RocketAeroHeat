function [delta,report]=assembleRocketWlsNormalCorrection(cache,maximumBlockBytes,assemblyBatchSize)
%ASSEMBLEROCKETWLSNORMALCORRECTION 提取完整WLS与两点近似之差的同原面法向线块。
% 模块一：保守量按当前EOS提升到速度/温度，保留特征样本的当前状态导数。
% 模块二：直接传播原WLS邻居及面位移修正；含应力、导热和两项黏性功。
% 模块三：只提取同原面线内的五条块带，不删除原CV或修改完整空间方程。
% 其他原面之间的WLS耦合未装入本近似；最终方向必须由原NS验收。
Q=cache.Q;cut=cache.geometry;gas=cache.gas;scale=cache.scale;f=cut.openFaces;n=size(Q,1);
if nargin<3,assemblyBatchSize=2048;end
validateattributes(assemblyBatchSize,{'numeric'},{'real','finite','scalar','integer','positive'});
nf=size(cut.inputGeometry.F,1);nl=n/nf;
assert(strcmp(cut.kind,'original-triangle-body-fitted-three-dimensional-shell')&&nl==floor(nl)&&all(f.cellMinus>0)&&all(cut.walls.ownerCell>0), ...
    'WlsNormal:Geometry','本模块要求全部正CV的完整原面层域。');
assert(strcmp(cutCellViscousWorkPolicy(cache.boundary),'arithmetic')&& ...
    strcmp(cutCellWallVelocityPolicy(cache.boundary.solidWall),'dirichlet-tangent')&& ...
    strcmp(cutCellWallThermalPolicy(cache.boundary.solidWall),'dirichlet-tangent')&& ...
    all(cache.boundary.solidWall.velocity==0)&&all(~cache.state.wall),'WlsNormal:Policy','需要已验证的算术功和静止Dirichlet墙策略。');
validateattributes(maximumBlockBytes,{'numeric'},{'real','finite','scalar','positive'});
blockBytes=5*5*5*n*8;assert(blockBytes<=maximumBlockBytes,'WlsNormal:Budget','五带块数据超预算，拒绝构造。');
v=cut.cells.volume(:);c=cut.cells.centroid;rho=Q(:,1);u=Q(:,2:4)./rho;cv=gas.R/(gas.gamma-1);
T=(Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);H=zeros(4,5,n);
for a=1:3,H(a,1,:)=reshape(-u(:,a)./rho*scale(1),1,1,n);H(a,a+1,:)=reshape(scale(a+1)./rho,1,1,n);end
dT=[(.5*sum(u.^2,2)-cv*T)./(rho*cv),-u./(rho*cv),1./(rho*cv)].*scale;
H(4,:,:)=permute(dT,[3 2 1]);primitiveLift=zeros(5,5,n);primitiveLift(1,1,:)=scale(1);primitiveLift(2:5,:,:)=H;
incoming=cache.gradient.incomingFaceIDs;incomingOwner=double(cache.gradient.incomingOwner);
assert(numel(unique(incomingOwner))==numel(incomingOwner),'WlsNormal:Boundary','当前完整层域每个外边界CV最多一个入流样本。');
[found,where]=ismember(incoming,cache.outsideFaceIDs);assert(all(found));
assert(numel(unique(cache.outsideOwner))==numel(cache.outsideOwner),'WlsNormal:Boundary','当前完整层域每个外边界CV最多一个外样本。');
allGivenLift=pagemtimes(cache.boundaryJacobian,primitiveLift(:,:,cache.outsideOwner));allGivenLift=allGivenLift(2:5,:,:);
givenLift=allGivenLift(:,:,where);outsideIndex=zeros(n,1);outsideIndex(cache.outsideOwner)=1:numel(cache.outsideOwner);
incomingIndex=zeros(n,1);incomingIndex(incomingOwner)=1:numel(incoming);clear primitiveLift;
delta=zeros(5,5,5,n);kappa=gas.mu*gas.gamma*cv/gas.Pr;inside=f.cellPlus>0;
sourceFace=mod((1:n)'-1,nf)+1;layers=floor(((1:n)'-1)/nf)+1;
sameLineMaximumOffset=0;
for d=1:3
    [row,col]=find(cache.gradient.gradientOperator{d});same=sourceFace(row)==sourceFace(col);
    if any(same),sameLineMaximumOffset=max(sameLineMaximumOffset,max(abs(layers(row(same))-layers(col(same)))));end
end
assert(sameLineMaximumOffset<=1,'WlsNormal:Stencil','WLS同线一环超出已证明的五带作用范围。');
active=find(f.area>0);timer=tic;

%% 模块四：按真实共享面及两端owner逐个积累同线五量块，不构造细层完整稠密Jacobian。
for endpoint=1:2
    faces=active;if endpoint==2,faces=faces(inside(faces));end
    for first=1:assemblyBatchSize:numel(faces)
        batch=faces(first:min(first+assemblyBatchSize-1,numel(faces)));m=double(f.cellMinus(batch));p=double(f.cellPlus(batch));
        if endpoint==1,rows=m;sign=1;else,rows=p;sign=-1;end
        for offset=-2:2
            columns=rows+offset*nf;valid=columns>0&columns<=n;ids=batch(valid);mi=m(valid);pi=p(valid);r=rows(valid);j=columns(valid);count=numel(j);
            if count==0,continue;end
            in=pi>0;normal=f.normal(ids,:);distanceVector=f.centroid(ids,:)-c(mi,:);
            distanceVector(in,:)=c(pi(in),:)-c(mi(in),:);distance=sum(distanceVector.*normal,2);assert(all(distance>0));
            L=localLift(mi,j);R=zeros(4,5,count);R(:,:,in)=localLift(pi(in),j(in));
            G=gradientBlock(mi,j);G(:,:, :,in)=.5*(G(:,:,:,in)+gradientBlock(pi(in),j(in)));
            if any(in),G(:,:,:,in)=correct(G(:,:,:,in),distanceVector(in,:),R(:,:,in)-L(:,:,in));end
            outside=find(~in);given=zeros(4,5,count);velocity=u(mi,:);velocity(in,:)=.5*(u(mi(in),:)+u(pi(in),:));
            dvelocity=.5*(L(1:3,:,:)+R(1:3,:,:));
            if ~isempty(outside)
                owner=mi(outside);which=outsideIndex(owner);same=j(outside)==owner;pick=outside(which>0&same);
                given(:,:,pick)=allGivenLift(:,:,which(which>0&same));
                incomingMask=cache.state.incoming(ids(outside));pick=outside(incomingMask);
                if ~isempty(pick),G(:,:,:,pick)=correct(G(:,:,:,pick),distanceVector(pick,:),given(:,:,pick)-L(:,:,pick));end
                pick=outside(~incomingMask);
                if ~isempty(pick)
                    nn=reshape(normal(pick,:)',1,3,1,[]);G(:,:,:,pick)=G(:,:,:,pick)-sum(G(:,:,:,pick).*nn,2).*nn;
                end
                velocity(outside,:)=cache.state.primitive(ids(outside),2:4);dvelocity(:,:,outside)=given(1:3,:,outside);
            end
            fullFlux=fluxBlock(G,normal,velocity,dvelocity,cache.baseInfo.viscousFaceFlux(ids,2:4));
            difference=R-L;nn=reshape(normal',1,3,1,count);
            twoPoint=reshape(difference,4,1,5,count).*nn./reshape(distance,1,1,1,count);
            frozenVelocity=u(mi,:);frozenVelocity(in,:)=.5*(u(mi(in),:)+u(pi(in),:));
            frozenFlux=fluxBlock(twoPoint,normal,frozenVelocity,zeros(3,5,count),zeros(count,3));
            contribution=(fullFlux-frozenFlux).*reshape(sign*f.area(ids)./v(r),1,1,count)./reshape(scale,5,1,1);
            [owners,~,group]=unique(r);
            for a=1:5,for b=1:5
                selected=accumarray(group,reshape(contribution(a,b,:),[],1),[numel(owners) 1]);
                delta(a,b,offset+3,owners)=delta(a,b,offset+3,owners)+reshape(selected,1,1,1,[]);
            end,end
        end
    end
end
assert(all(isfinite(delta),'all'),'WlsNormal:Nonfinite','同线黏性修正非有限，不提交裁剪结果。');
report=struct('kind','exact-selected-WLS-minus-frozen-two-point-viscous-same-original-facet-five-band-correction', ...
    'sourceFaces',nf,'layers',nl,'cellCount',n,'maximumSameLineGradientLayerOffset',sameLineMaximumOffset, ...
    'blockStorageUpperBoundBytes',blockBytes,'maximumBlockBytes',maximumBlockBytes,'assemblyBatchSize',assemblyBatchSize,'assemblySeconds',toc(timer), ...
    'fullSelectedModelApplications',0,'newActualNSEvaluations',0,'crossOriginalFacetWlsBlocksIncluded',false, ...
    'stressHeatAndBothWorkTermsIncluded',true,'wallCorrectionExactlyZeroForSelectedPolicy',true, ...
    'physicalEquationChanged',false,'newProductionSolvedAreaM2',0);
    function lift=localLift(cells,columns)
        lift=zeros(4,5,numel(cells));same=cells==columns;lift(:,:,same)=H(:,:,columns(same));
    end
    function G=gradientBlock(cells,columns)
        gCount=numel(cells);G=zeros(4,3,5,gCount);gIndex=sub2ind([n n],cells,columns);
        for gDimension=1:3
            gCoefficient=full(cache.gradient.gradientOperator{gDimension}(gIndex));
            G(:,gDimension,:,:)=reshape(H(:,:,columns).*reshape(gCoefficient,1,1,[]),4,1,5,gCount);
            gLocal=incomingIndex(cells);gPick=find(cells==columns&gLocal>0);
            if ~isempty(gPick)
                gExtra=givenLift(:,:,gLocal(gPick)).*reshape(cache.gradient.boundaryCoefficientXYZ(gLocal(gPick),gDimension),1,1,[]);
                G(:,gDimension,:,gPick)=G(:,gDimension,:,gPick)+reshape(gExtra,4,1,5,[]);
            end
        end
    end
    function G=correct(G,displacement,difference)
        cCount=size(G,4);cDx=reshape(displacement',1,3,1,cCount);cD2=sum(displacement.^2,2);
        cMismatch=reshape(difference,4,1,5,cCount)-sum(G.*cDx,2);G=G+cMismatch.*cDx./reshape(cD2,1,1,1,cCount);
    end
    function out=fluxBlock(G,normal,velocity,dvelocity,baseTraction)
        fCount=size(G,4);out=zeros(5,5,fCount);
        for fColumn=1:5
            fVelocityGradient=reshape(G(1:3,:,fColumn,:),3,3,fCount);fStress=gas.mu*(fVelocityGradient+pagetranspose(fVelocityGradient));
            fTrace=fVelocityGradient(1,1,:)+fVelocityGradient(2,2,:)+fVelocityGradient(3,3,:);
            for fAxis=1:3,fStress(fAxis,fAxis,:)=fStress(fAxis,fAxis,:)-(2/3)*gas.mu*fTrace;end
            fTraction=reshape(pagemtimes(fStress,permute(normal,[2 3 1])),3,fCount)';out(2:4,fColumn,:)=permute(fTraction,[2 3 1]);
            fHeat=kappa*sum(reshape(G(4,:,fColumn,:),3,fCount)'.*normal,2);
            fPower=sum(velocity.*fTraction,2)+sum(reshape(dvelocity(:,fColumn,:),3,fCount)'.*baseTraction,2);
            out(5,fColumn,:)=reshape(fHeat+fPower,1,1,fCount);
        end
    end
end
