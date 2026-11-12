function path=marchPhaseAwareWavePath(p,T,U,direction,tangent,planeNormal)
%MARCHPHASEAWAREWAVEPATH 同模型二维折线波系逐段传递，禁止每段重置自由流。
% 模块用途：适用于具名上游超声速状态及同一平面内的折线/表面剖面。
% planeNormal规定正压缩方向；不含轴对称几何源、三维横流及黏性入口。
% 未知/脱体/超域一旦发生，全部下游保持未知；不能跳到后续有效段重启。
validateattributes(tangent,{'numeric'},{'real','finite','ncols',3,'nonempty'});
validateattributes(direction,{'numeric'},{'real','finite','size',[1 3]});
validateattributes(planeNormal,{'numeric'},{'real','finite','size',[1 3]});
assert(abs(norm(direction)-1)<1e-10&&abs(norm(planeNormal)-1)<1e-10 && ...
    all(abs(vecnorm(tangent,2,2)-1)<1e-10),'RAH:PhasePath:Unit','方向必须为单位向量。');
assert(abs(dot(direction,planeNormal))<1e-10&&all(abs(tangent*planeNormal')<1e-10), ...
    'RAH:PhasePath:Planar','输入不在同一波系平面，不能把三维偏转当二维简单波。');
initial=phaseAwareAirState('PT',p,T,'linearTable');initialSpeed=U;
n=size(tangent,1);rows=cell(n,1);reason=repmat({'上游受阻，不允许重置自由流'},n,1);
values=nan(n,3);turn=nan(n,1);valid=false(n,1);last=direction;
for i=1:n
    angle=atan2d(dot(planeNormal,cross(last,tangent(i,:))),dot(last,tangent(i,:)));
    turn(i)=angle;
    if abs(angle)>=90,reason{i}='转角达到或超过90度，入口需要脱体/亚声速鼻部解';break;end
    kind='oblique';if angle<0,kind='expansion';end
    w=solvePhaseAwareWave(kind,p,T,U,abs(angle));rows{i}=w;reason{i}=w.status;
    if ~w.available,break;end
    p=w.pressurePa;T=w.temperatureK;U=w.speedMps;values(i,:)=[p T U];
    valid(i)=true;last=tangent(i,:);
end
path=struct('schema','rah-phase-wave-path-v1','pressurePa',values(:,1),'temperatureK',values(:,2), ...
    'speedMps',values(:,3),'velocityMps',values(:,3).*tangent,'turnDeg',turn, ...
    'valid',valid,'reason',{reason},'waves',{rows},'model',initial.model, ...
    'transportModel',initial.transportModel,'energyReference',initial.energyReference, ...
    'initialState',initial,'initialSpeedMps',initialSpeed, ...
    'direction',direction,'planeNormal',planeNormal,'tangent',tangent, ...
    'physicalQualification',false,'boundaryLayerOriginProvided',false, ...
    'scope','二维局部顺序波系候选；无三维横流/曲率源及边界层入口，不可直接授热资格');
end
