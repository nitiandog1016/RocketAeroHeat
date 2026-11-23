function check_ui(out)
% 模块一：用普通 STL 入口生成体网，把重剖子面与文件真值逐行核对。
root=fileparts(which('rocketAeroHeatWorkbench'));fixture=fullfile(root,'outputs','v320_validation','fixture');
truth=readmatrix(fullfile(fixture,'face_truth.csv'));ids=find(truth(:,2)>0);
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(fixture,'probe-partition-mm.stl');c.scale.Value=.001;press(c.load);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'thermistorB');
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh),'%s',fig.UserData.status.Text);b=s.mesh.thermistorImport;
assert(isequal(b.referenceFaceID,truth(:,1)) && isequal(ids,b.childIDs) && isequal(s.mesh.F,s.rawMesh.F));
proof=struct('stlSHA256',b.source.sha256,'prepared',s.mesh.thermistorPrepared.directory,'faces',size(truth,1),'children',ids);
press(fig.UserData.controls.nextBoundary);
% 模块二：实测三个热边界及各自的求解来源，保存完整面的科学 CSV。
modes={'longWire','resistance','resistance','uniform','uniform'};
actions={'archive','solve','archive','solve','archive'};
for k=1:numel(modes)
    name=[modes{k} '_' actions{k}];dest=fullfile(out,name);mkdir(dest);
    choose(fig.UserData.controls.wireBoundary,modes{k});press(fig.UserData.controls.confirm);
    if strcmp(actions{k},'solve'),press(fig.UserData.controls.run);else,press(fig.UserData.controls.readBead);end
    s=fig.UserData.workflowState;assert(s.currentStep==5,'%s',fig.UserData.status.Text);
    res=s.result;v=fig.UserData.webUI.evidenceView;html=fig.UserData.webUI.html;
    timer=tic;while isempty(html.UserData.rendered) && toc(timer)<40,drawnow;pause(.1);end
    gl=html.UserData.rendered;assert(gl.faces==size(truth,1) && gl.webglError==0 && gl.visiblePixelSamples>0);
    v.selectField('heatFlux');
    for i=ids',v.inspectFace(i,[]);u=v.getState();assert(u.selectedFaceID==i);end
    wall=v.exportWallPairTo(fullfile(dest,'wall.csv'));assert(wall.rows==size(truth,1) && wall.joint.pressureHeatFluxQualifiedArea==0);
    if k==1
        % 分别穿过前接头和前引线，导出实际浏览器传回的两条截线。
        for item=[183 184]
            x=.000574;if item==184,x=.000582;end
            v.setSectionX(x);timer=tic;
            while toc(timer)<40
                drawnow;cut=v.getSection();
                if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-13,break;end
                pause(.1);
            end
            assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'heatFlux')&&abs(cut.x-x)<1e-13);
            v.exportSectionTo(fullfile(dest,sprintf('section%d.csv',item)));
        end
        exportAerothermalWeb(res,fullfile(dest,'workbench'));
    end
    proof.(name)=struct('frame',res.evidenceIdentity.sourceFrame,'shapeError',res.thermistorSubfaces.shapeError, ...
        'regions',fig.UserData.caseReport.regions,'powerW',fig.UserData.caseReport.wholeSurfaceHeatPowerW, ...
        'glPixels',gl.visiblePixelSamples,'glError',gl.webglError);
    if k<numel(modes),press(fig.UserData.workflowButtons(3));end
end
% 模块三：缺口和等面积重叠必须分别拒绝，不允许跨越接头和引线。
[ref,~]=thermistorWorkbenchModel('leads');g=s.orientedMesh;
y=g;y.F(ids(1),:)=[];proof.gap=reject(y,ref,'RAH:Bead:Coverage');
pair=ids(truth(ids,1)==183);y=ref;corners=ref.F(183,:);new=size(ref.V,1)+1;
y.V(new,:)=mean(ref.V(corners,:),1);y.F(183,:)=[];
fan=[corners' corners([2 3 1])' repmat(new,3,1)];y.F=[y.F;fan([1 1 3],:)];
proof.overlap=reject(y,ref,'RAH:Bead:Coverage');
y=g;far=ids(truth(ids,1)==184);y.F(pair(1),:)=[y.F(pair(1),[1 3]) y.F(far(1),3)];proof.cross=reject(y,ref,'RAH:Bead:Triangle');
y=g;y.F(pair(1),[2 3])=y.F(pair(1),[3 2]);proof.winding=reject(y,ref,'RAH:Bead:NormalArea');
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function message=reject(mesh,reference,code)
try,mapThermistorSTL(mesh,reference);
catch ex,assert(strcmp(ex.identifier,code),'%s: %s',ex.identifier,ex.message);message=ex.message;return;end
error('负例未被拒绝。');
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
