function report=diagnoseObliqueStreamlineConnectivity(out)
%DIAGNOSEOBLIQUESTREAMLINECONNECTIVITY 比较原中点跨边和边内采样的连通性。
% 本诊断只测有向图可达性，不重新播种，不改变物理判据或宣称热流有效。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(tempdir,'RocketAeroHeat_connectivity');end
if ~exist(out,'dir'),mkdir(out);end
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=20;cfg.flight.beta=0;cfg.flight.Twall=300;
report=struct([]);
for kind={'hemisphereCylinder','bluntedCone'}
    for resolution=[24 48]
        [mesh,meta]=buildCanonicalBody(kind{1},struct('resolution',resolution));
        evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
        mesh=analyzeLocalGeometry(mesh,cfg);cfg.stagnation.physicalNoseRadius=meta.noseRadius;
        mesh=calcFlightStationState(mesh,cfg.flight,cfg);
        s=mesh.surfaceStreamlines;n=size(mesh.F,1);C=mesh.localGeometry.center;D=s.tangentDirection;
        F=mesh.F;V=mesh.V;
        [edges,~,group]=unique(sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2),'rows');
        count=accumarray(group,1);face=repmat((1:n)',3,1);
        lo=accumarray(group,face,[],@min);hi=accumarray(group,face,[],@max);ok=count==2;
        from=[lo(ok);hi(ok)];to=[hi(ok);lo(ok)];edge=[edges(ok,:);edges(ok,:)];
        A=V(edge(:,1),:);B=V(edge(:,2),:);cost=inf(numel(from),1);bestParameter=nan(size(cost));
        for u=linspace(0,1,81)
            P=A+u*(B-A);a=P-C(from,:);b=C(to,:)-P;
            la=vecnorm(a,2,2);lb=vecnorm(b,2,2);
            ca=sum(D(from,:).*a,2)./la;cb=sum(D(to,:).*b,2)./lb;
            keep=ca>=cosd(s.maxStepDeviationDeg)&cb>=cosd(s.maxStepDeviationDeg);
            value=la./ca+lb./cb;improve=keep&value<cost;
            cost(improve)=value(improve);bestParameter(improve)=u;
        end
        keep=isfinite(cost);source=n+1;
        graphSample=digraph([from(keep);repmat(source,numel(s.seedFaceIDs),1)], ...
            [to(keep);s.seedFaceIDs],[cost(keep);zeros(numel(s.seedFaceIDs),1)],source);
        d=distances(graphSample,source);sampleConnected=isfinite(d(1:n))';
        frontier=s.connectedFace(from)&~s.connectedFace(to);
        row=struct('name',sprintf('%s_n%d_a20',kind{1},resolution),'faces',n, ...
            'seedIDs',s.seedFaceIDs,'seedSource',s.seedSource,'candidateFaces',nnz(mesh.stagnation.candidateFace), ...
            'candidateValid',nnz(mesh.stagnation.candidateFace & mesh.boundaryLayerEdge.algebraicallyValid), ...
            'connectedMidpoint',nnz(s.connectedFace),'connectedSampled',nnz(sampleConnected), ...
            'frontierEdges',nnz(frontier),'sampledFrontierEdges',nnz(frontier&keep), ...
            'heatAreaFraction',mesh.heatingCoverage.summary.computedAreaFraction);
        if isempty(report),report=row;else,report(end+1)=row;end
        fprintf('%s: candidates %d, seeds %d, midpoint %d, sampled %d; frontier pass %d/%d\n', ...
            row.name,row.candidateFaces,numel(row.seedIDs),row.connectedMidpoint,row.connectedSampled, ...
            row.sampledFrontierEdges,row.frontierEdges);
        save(fullfile(out,[row.name '.mat']),'mesh','cfg','meta','from','to','edge','cost','bestParameter','sampleConnected');
    end
end
fid=fopen(fullfile(out,'connectivity.json'),'w','n','UTF-8');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
end
