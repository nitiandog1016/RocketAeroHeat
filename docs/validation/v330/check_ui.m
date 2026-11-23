function check_ui(out)
% 模块一：真实侧栏打开近壁对照，从五区差值探针切换到两档原场，再返回压力截面。
fig=rocketAeroHeatWorkbench();done=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supportResolution';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v) && isfield(v,'comparison'),'%s',fig.UserData.status.Text);
html=v.web.html;t=tic;while isempty(html.UserData.rendered) && toc(t)<65,drawnow;pause(.15);end
gl=html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==32480);
e=v.web.evidenceView;d=e.getData();p=v.mesh.shortSupportGridPair;
assert(isequal(d.comparisonColorRanges.heatFlux,d.comparisonColorRanges.supportCoarseHeat));
assert(isequal(d.comparisonColorRanges.pressureGauge,d.comparisonColorRanges.supportCoarseGauge));
for key={'heatFlux','supportCoarseHeat','supportCoarseGauge','supportGaugeChange','supportHeatChange'}
    e.selectField(key{1});
    for r=v.report.regions(:)'
        e.inspectFace(r.peakOutFace,[]);state=e.getState();
        assert(strcmp(state.fieldKey,key{1}) && state.selectedFaceID==r.peakOutFace);
    end
end
e.selectField('pressureGauge');x=.00137;e.setSectionX(x);t=tic;
while toc(t)<50
    drawnow;cut=e.getSection();
    if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressureGauge')&&abs(cut.x-x)<1e-13,break;end
    pause(.1);
end
assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressureGauge')&&abs(cut.x-x)<1e-13);
e.exportSectionTo(fullfile(out,'section.csv'));w=e.exportWallPairTo(fullfile(out,'wall.csv'));
assert(w.rows==32480 && w.joint.pressureHeatFluxQualifiedArea==0);
v.tabs.SelectedTab=v.comparison.tab;drawnow;assert(height(v.comparison.table.Data)==5 && numel(v.comparison.line.YData)==7);
page=exportAerothermalWeb(v.mesh,fullfile(out,'workbench'));
proof=struct('coarseFrame',p.coarseFrame,'fineFrame',p.fineFrame,'faces',gl.faces,'pixels',gl.visiblePixelSamples, ...
    'commonHeatScale',p.fluxRange,'commonGaugeScale',p.gaugeRange,'page',page);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
