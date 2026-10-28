function report=testCutCellExactEulerFlux(folder)
%TESTCUTCELLEXACTEULERFLUX 解析波系/经典Sod及真实收缩面独立Python参考，保留温域/真空拒绝。
% 模块一：SI缩放Sod、对称强稀疏波/激波、接触/均匀通量和三维旋转/反向一致。
% 模块二：直接读取v81真实两侧Q/n及公开曲线参考，不以自画场替代真实点。
% 模块三：标明解析/数值校核，不认证完整NS、热流、实验或物理解面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CutExactTest:Overwrite','不得覆盖真实测试。');mkdir(folder);
gas=struct('gamma',1.4,'R',287,'temperatureRange',[80 2000]);g=gas.gamma;normal=[1 0 0];p0=gas.R*300;
QL=cons([1 0 50 -20 300]);QR=cons([.125 0 -30 70 240]);[sod,sa]=cutCellExactEulerFlux(QL,QR,normal,gas);
assert(abs(sa.pressureStarPa/p0-.30313017805064685)<2e-12&&abs(sa.normalVelocityStarMs/sqrt(p0)-.9274526200489499)<2e-12);
equal=cons([1.2 -100 30 75 400]);[fe,~]=cutCellExactEulerFlux(equal,equal,normal,gas);expected=physical(equal,normal);
assert(max(abs(fe-expected))/max(1,max(abs(expected)))<2e-13);
QL=cons([1 35 10 20 300]);QR=cons([.4 35 -30 90 750]);[contact,ca]=cutCellExactEulerFlux(QL,QR,normal,gas);
assert(max(abs(contact-physical(QL,normal)))/max(abs(contact))<2e-12&&abs(ca.normalVelocityStarMs-35)<1e-10);
QL=cons([1 -600 0 0 300]);QR=cons([1 600 0 0 300]);[rare,ra]=cutCellExactEulerFlux(QL,QR,normal,gas);
a=sqrt(g*p0);pRare=p0*(1-(g-1)*1200/(4*a))^(2*g/(g-1));assert(abs(ra.pressureStarPa/pRare-1)<2e-12&&abs(rare(2)/pRare-1)<2e-12);
QL=cons([1 600 0 0 300]);QR=cons([1 -600 0 0 300]);[shock,sh]=cutCellExactEulerFlux(QL,QR,normal,gas);
A=2/(g+1);B=(g-1)/(g+1)*p0;v2=600^2;coefficient=2*p0+v2/A;
pShock=.5*(coefficient+sqrt(coefficient^2-4*(p0^2-v2*B/A)));
assert(abs(sh.pressureStarPa/pShock-1)<2e-12&&abs(shock(2)/pShock-1)<2e-12);

%% 真实受限面及独立Python/公开压力曲线参考
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
s=load(fullfile(root,'outputs','v81_validation','arithmetic02','steady_krylov_evidence.mat'),'result','gas');
ref=jsondecode(fileread(fullfile(root,'outputs','v82_validation','actual_exact_riemann_comparison.json')));
face=70995;f=d.domain.openFaces;m=f.cellMinus(face);p=f.cellPlus(face);n=f.normal(face,:);ql=s.result.Q(m,:);qr=s.result.Q(p,:);
[actual,aa]=cutCellExactEulerFlux(ql,qr,n,s.gas);difference=max(abs(actual-ref.actualExactRiemann.flux(:)'))/max(abs(actual));
assert(difference<2e-12&&abs(aa.pressureStarPa/ref.actualExactRiemann.pressureStar-1)<2e-12&&aa.sampleRegion==3);
[reverse,~]=cutCellExactEulerFlux(qr,ql,-n,s.gas);assert(max(abs(reverse+actual))/max(abs(actual))<2e-12);
angle=.47;rotation=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];qll=ql;qrr=qr;qll(2:4)=(rotation*ql(2:4)')';qrr(2:4)=(rotation*qr(2:4)')';
[rotated,~]=cutCellExactEulerFlux(qll,qrr,(rotation*n')',s.gas);expected=actual;expected(2:4)=(rotation*actual(2:4)')';
assert(max(abs(rotated-expected))/max(abs(expected))<2e-12);
[batch,ba]=cutCellExactEulerFlux([ql;qr],[qr;ql],[n;-n],s.gas);assert(max(abs(batch-[actual;reverse]),[],'all')/max(abs(actual))<2e-12);
bad={{cons([1 -2000 0 0 300]),cons([1 2000 0 0 300]),normal,'CutExactEuler:Vacuum'}, ...
     {ql,qr,2*n,'CutExactEuler:Geometry'}, ...
     {[1 0 0 0 -1],qr,n,'CutExactEuler:InternalEnergy'}, ...
     {cons([1 10 20 30 70]),qr,n,'CutExactEuler:Temperature'}};
for j=1:numel(bad)
    z=bad{j};rejected=false;try,cutCellExactEulerFlux(z{1},z{2},z{3},s.gas);catch ex,rejected=strcmp(ex.identifier,z{4});end
    assert(rejected,'CutExactTest:Rejection','非法/未支持输入未被明确拒绝。');
end
report=struct('status','passed-SI-Sod-strong-rarefaction-shock-contact-and-real-contraction-exact-flux-reference', ...
    'actualFace',face,'actualPressureStarPa',aa.pressureStarPa,'actualSampleTemperatureK',aa.samplePrimitive(5), ...
    'actualFiveQuantityFluxRelativeDifference',difference,'actualRootIterations',aa.pressureRootIterations, ...
    'actualSampleRegion','left-rarefaction-fan','explicitRejections',numel(bad),'reverseRotationBatchChecked',true, ...
    'fullThreeDimensionalNSValidated',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'exact_flux_evidence.mat'),'report','actual','aa','ql','qr','n','gas','sod','sa','rare','ra','shock','sh','contact','ca','batch','ba','-v7');
fid=fopen(fullfile(folder,'exact_flux_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function q=cons(P)
        q=[P(1),P(1)*P(2:4),P(1)*(gas.R/(gas.gamma-1)*P(5)+.5*sum(P(2:4).^2))];
    end
    function f=physical(q,n)
        rho=q(1);u=q(2:4)/rho;p=(gas.gamma-1)*(q(5)-.5*rho*sum(u.^2));un=dot(u,n);
        f=[rho*un,q(2:4)*un+p*n,(q(5)+p)*un];
    end
end
