function testHybridWallCalculation(folder)
%TESTHYBRIDWALLCALCULATION 单工况公共入口、局部回退、真实稳态和独立解析对照。
[r,mesh,cfg]=demoHybridWallCalculation(folder);c=cfg.hybridSurface.condition;
w=readtable(fullfile(folder,'wall.csv'),'TextType','string');
s=load(fullfile(folder,'states.mat'),'result');ns=s.result;
assert(isequal(w.faceID,(0:7)') && r.faceCount==8 && r.totalArea==4 && r.engineeringFirstAreaM2==2);
assert(r.regionalAddedAreaM2==2 && r.jointComputedAreaM2==4 && r.fullSurfaceComputed && numel(r.attempts)==1);
assert(isequal(r.attempts(1).targetFaceIDs,(4:7)') && ns.acceptedSteps>0 && ~isequal(ns.Q,ns.initialState));
assert(all(ns.history.residualNorm(end-7:end)<=2e-6) && ns.globalBalanceRelative<1e-10);
initialFile=load(fullfile(folder,'initial.mat'),'initial');assert(isequal(initialFile.initial,ns.initialState));
for step=[1 50:50:ns.acceptedSteps]
    frame=load(fullfile(folder,sprintf('step%04d.mat',step)));
    assert(frame.step==step && frame.time==ns.history.time(step+1) && isequal(size(frame.Q),size(ns.Q)));
end
% 模块一：独立用质量积分反推两域压力，使用剪切机械功及线性导热核对逐面符号。
k=c.mu*c.gamma*c.R/((c.gamma-1)*c.Pr);
p0=c.meanDensity*c.R/integral(@(z)1./(300+30*z),0,1);
p1=c.meanDensity*c.R/integral(@(z)1./(300+c.mu*c.upperSpeedMps^2/(2*k)*z.*(1-z)),0,1);
pExpected=[repmat(p0,4,1);repmat(p1,4,1)];
shearPower=c.mu*c.upperSpeedMps^2/c.heightM;
qExpected=[k*30/c.heightM*[1;1;-1;-1];repmat(shearPower/2,4,1)];
pErr=max(abs(w.pressurePa./pExpected-1));qErr=max(abs(w.heatFluxWm2./qExpected-1));
assert(pErr<1e-4 && qErr<1e-5 && abs(r.wholeSurfaceHeatPowerW/shearPower-1)<1e-5);
% 模块二：核实送给科研工作台的值、给定壁温、独立验收及来源身份。
js=strtrim(fileread(fullfile(folder,'workbench','case-data.js')));
prefix='window.ROCKET_CASE=';assert(startsWith(js,prefix)&&endsWith(js,';'));
d=jsondecode(js(numel(prefix)+1:end-1));
p=d.fields(strcmp({d.fields.key},'pressure'));q=d.fields(strcmp({d.fields.key},'heatFlux'));
assert(max(abs(1000*p.values-w.pressurePa))<1e-9 && max(abs(1000*q.values-w.heatFluxWm2))<1e-9);
assert(all(p.status==1) && all(q.status==1) && numel(q.methodNames)==2 && strcmp(d.evidenceIdentity.origin,'research'));
tw=d.fields(strcmp({d.fields.key},'wallInput'));assert(all(tw.status==2)&&isequal(tw.values,cfg.hybridSurface.wallTemperatureK));
% 模块三：未收敛结果保持缺项；合法零、已有工程半解和乱序编号分别检验。
cfg.hybridSurface.localSolver=@(varargin)struct('accepted',false,'reason','bounded rejection');
[~,~,partial]=runAeroHeatSingleCase(mesh,cfg);
assert(~partial.fullSurfaceComputed && isnan(partial.wholeSurfaceHeatPowerW) && partial.jointComputedAreaM2==2);
onlyP=table(4,315,NaN,"test computed pressure",'VariableNames',w.Properties.VariableNames);
eng=cfg.hybridSurface.engineeringSolver;
cfg.hybridSurface.engineeringSolver=@(g,c,t)[eng(g,c,t);onlyP];
cfg.hybridSurface.localSolver=@(varargin)struct('accepted',true,'reason','partial field', ...
    'wall',table(4,NaN,0,"test computed flux",'VariableNames',w.Properties.VariableNames));
[m,~,partialFields]=runAeroHeatSingleCase(mesh,cfg);
assert(m.hybridWall.pressurePa(5)==315 && isnan(m.hybridWall.heatFluxWm2(5)) && partialFields.jointComputedAreaM2==2);
localZero=table(4,301,0,"test zero flux",'VariableNames',w.Properties.VariableNames);
cfg.hybridSurface.localSolver=@(varargin)struct('accepted',true,'reason','test','wall',localZero);
[m,~,zero]=runAeroHeatSingleCase(mesh,cfg);
assert(m.hybridWall.heatFluxWm2(5)==0 && m.hybridWall.pressurePa(5)==301 && zero.jointComputedAreaM2==2.5);
badRows=localZero; badRows.faceID=0;
cfg.hybridSurface.localSolver=@(varargin)struct('accepted',true,'reason','invalid region','wall',badRows);
try
    runAeroHeatSingleCase(mesh,cfg);
    error('RAH:Test:ExpectedRejection','混入其他区域的编号应被拒绝。');
catch ex
    assert(strcmp(ex.identifier,'RAH:Hybrid:FaceID'),ex.message);
end
summary=struct('fullSurfaceComputed',r.fullSurfaceComputed,'totalWettedAreaM2',r.totalArea, ...
    'engineeringFirstAreaM2',r.engineeringFirstAreaM2,'regionalAddedAreaM2',r.regionalAddedAreaM2, ...
    'jointComputedAreaM2',r.jointComputedAreaM2,'wholeSurfaceHeatPowerW',r.wholeSurfaceHeatPowerW, ...
    'acceptedSteps',ns.acceptedSteps,'initialResidual',ns.initialResidualNorm,'lastEightResidual',ns.history.residualNorm(end-7:end), ...
    'globalBalanceRelative',ns.globalBalanceRelative,'maxPressureRelativeError',pErr,'maxHeatRelativeError',qErr, ...
    'modelScope','two isolated periodic laboratory gaps, constant-property test gas', ...
    'experimentalValidation',false,'rocketCoverageGainM2',0,'negativeControlsPassed',true);
fid=fopen(fullfile(folder,'result.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));
end
