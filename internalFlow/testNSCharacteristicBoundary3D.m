function testNSCharacteristicBoundary3D(outputFolder)
%TESTNSCHARACTERISTICBOUNDARY3D 以独立不变量和通量检验开放边界草稿。
% 校核范围：恒定理想气体法向边界、旋转/反向、回流与外推合同；
% 不代表非反射性质、真实箭体尾迹或试验热流精度已验证。
addpath(fileparts(mfilename('fullpath')),'-begin');
gas=struct('gamma',1.4,'R',287,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]);
g=gas.gamma;n=[1 0 0];c=[1.2 100 12 -8 300];
b=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',c);
maximumUniformError=0;

%% 1. 入射/出射特征数量和均匀状态保持；两侧方向均检查
for un=[-600 -100 0 100 600]
    local=c;local(2)=un;b.reservoirPrimitive=local;
    for normal=[n;-n]'
        [face,info]=resolveNSCharacteristicBoundary3D(local,b,normal',gas);
        error=max(abs(face-local)./max(1,abs(local)));
        maximumUniformError=max(maximumUniformError,error);
        assert(error<2e-13);
        unExpected=dot(local(2:4),normal');a=sqrt(g*gas.R*local(5));
        count=(unExpected-a<0)+3*(unExpected<0)+(unExpected+a<0);
        assert(info.incomingCharacteristics==count && ~info.pressureInputApplied);
    end
end

%% 2. 压力出口保持内部熵、切向速度和出射声学不变量
b.mode='pressureOutlet';b.pressure=1.1*c(1)*gas.R*c(5);b.reservoirPrimitive=c;
[face,out]=resolveNSCharacteristicBoundary3D(c,b,n,gas);
p=c(1)*gas.R*c(5);pb=face(1)*gas.R*face(5);
assert(abs(pb/b.pressure-1)<2e-13 && out.pressureInputApplied);
assert(abs((pb/face(1)^g)/(p/c(1)^g)-1)<2e-13);
Jc=c(2)+2*sqrt(g*gas.R*c(5))/(g-1);
assert(abs((face(2)+2*sqrt(g*gas.R*face(5))/(g-1))/Jc-1)<2e-13);
assert(max(abs(face(3:4)-c(3:4)))==0);

%% 3. 明确储库回流；pressure输入不参与该分支，不双重强加边界
in=c;in(2)=-30;reservoir=[.9 -40 -2 7 350];b.reservoirPrimitive=reservoir;
b.pressure=1e3;[back,backInfo]=resolveNSCharacteristicBoundary3D(in,b,n,gas);
b.pressure=1e8;[back2,backInfo2]=resolveNSCharacteristicBoundary3D(in,b,n,gas);
assert(isequal(back,back2) && ~backInfo.pressureInputApplied && ~backInfo2.pressureInputApplied);
assert(backInfo.isInflow && backInfo.incomingCharacteristics==4);
Jpi=in(2)+2*sqrt(g*gas.R*in(5))/(g-1);
Jmr=reservoir(2)-2*sqrt(g*gas.R*reservoir(5))/(g-1);
ab=sqrt(g*gas.R*back(5));pb=back(1)*gas.R*back(5);pr=reservoir(1)*gas.R*reservoir(5);
assert(abs(back(2)+2*ab/(g-1)-Jpi)<1e-10);
assert(abs(back(2)-2*ab/(g-1)-Jmr)<1e-10);
assert(abs((pb/back(1)^g)/(pr/reservoir(1)^g)-1)<2e-13);
assert(isequal(back(3:4),reservoir(3:4)));

%% 4. 按面流向重分类，以及不自洽分支/声速/非正声速明确拒绝
c20=c;c20(2)=20;b.pressure=1.2*p;b.reservoirPrimitive=[c(1),-50,c(3:5)];
[flipped,flipInfo]=resolveNSCharacteristicBoundary3D(c20,b,n,gas);
assert(flipped(2)<0 && flipInfo.reclassifiedByFaceVelocity && ~flipInfo.pressureInputApplied);
b.reservoirPrimitive=[c(1),80,c(3:5)];
reject(@()resolveNSCharacteristicBoundary3D(c20,b,n,gas),'CharacteristicAmbiguous');
sonic=c;sonic(2)=sqrt(g*gas.R*c(5));
reject(@()resolveNSCharacteristicBoundary3D(sonic,b,n,gas),'CharacteristicSonic');
bad=c;bad(5)=0;reject(@()resolveNSCharacteristicBoundary3D(bad,b,n,gas),'CharacteristicState');
missing=rmfield(b,'reservoirPrimitive');
reject(@()resolveNSCharacteristicBoundary3D(c,missing,n,gas),'CharacteristicInput');

%% 5. 任意坐标旋转下法向、速度与标量状态协变
theta=.73;rotation=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
b.mode='riemannFarfield';b.reservoirPrimitive=reservoir;
base=in;[face,~]=resolveNSCharacteristicBoundary3D(base,b,n,gas);
rotated=base;rotated(2:4)=base(2:4)*rotation';br=b;
br.reservoirPrimitive(2:4)=reservoir(2:4)*rotation';
[rotFace,~]=resolveNSCharacteristicBoundary3D(rotated,br,n*rotation',gas);
expected=face;expected(2:4)=face(2:4)*rotation';
assert(max(abs(rotFace-expected)./max(1,abs(expected)))<5e-13);

%% 6. 实际三维空间内核的开放边界：均匀流及五量全域平衡
shape=[3 2 2];grid=struct('spacing',[.02 .01 .03]);solid=false(shape);
cv=gas.R/(g-1);uniformResidual=0;uniformBalance=0;
for speed=[-100 100]
    input=c;input(2)=speed;input(3:4)=0;
    q=[input(1),input(1)*input(2:4),input(1)*(cv*input(5)+.5*sum(input(2:4).^2))];
    Q=repmat(reshape(q,1,1,1,5),[shape,1]);
    far=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',input);
    outlet=far;outlet.mode='pressureOutlet';outlet.pressure=input(1)*gas.R*input(5);
    periodic=struct('type','periodic');
    boundary=struct('xmin',far,'xmax',outlet,'ymin',periodic,'ymax',periodic,'zmin',periodic,'zmax',periodic);
    [rate,meta]=compressibleNavierStokesResidual3D(Q,grid,solid,gas,boundary);
    reference=max(abs(q))*(abs(speed)+sqrt(g*gas.R*input(5)))*sum(1./grid.spacing);
    relative=max(abs(rate(:)))/reference;assert(relative<5e-13);
    assert(numel(meta.openFaces)==8 && all(meta.faceKind{1}([1 end],:,:,:)==6,'all'));
    assert(nnz(arrayfun(@(f)f.normal(1)<0,meta.openFaces))==4);
    assert(nnz(arrayfun(@(f)f.normal(1)>0,meta.openFaces))==4);
    assert(all(arrayfun(@(f)isequal(f.normal,f.boundaryState.normal),meta.openFaces)), ...
        '外边界面法向记录必须与实际特征构造法向一致。');
    % 两开放面的净通量接近零，不能用净值作舍入误差尺度。
    % 独立均匀流解析穿流量为F(q)*端面面积，以双端绝对穿流归一各守恒量。
    pressure=input(1)*gas.R*input(5);planeArea=prod(shape(2:3).*grid.spacing(2:3));
    exactFlux=[q(2),q(2)*speed+pressure,0,0,(q(5)+pressure)*speed];
    balance=max(abs(meta.balanceResidual)./max(1,2*planeArea*abs(exactFlux)));
    assert(balance<5e-13);uniformResidual=max(uniformResidual,relative);uniformBalance=max(uniformBalance,balance);
end
fprintf('Characteristic boundary passed: uniform %.3g, spatial %.3g, global balance %.3g; explicit backflow, pressure identity and rotation checked.\n', ...
    maximumUniformError,uniformResidual,uniformBalance);
if nargin>0
    record=struct('uniformStateRelativeError',maximumUniformError, ...
        'spatialUniformRelativeError',uniformResidual,'globalBalanceRelativeError',uniformBalance, ...
        'physicalRegionSolved',false,'experimentalValidation',false);
    fid=fopen(fullfile(outputFolder,'characteristic_test_summary.json'),'w','n','UTF-8');
    closer=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(record),'char'); %#ok<NASGU>
end
end

function reject(f,identifier)
failed=false;try,f();catch err,failed=endsWith(err.identifier,identifier);end
assert(failed,['Expected ',identifier]);
end
