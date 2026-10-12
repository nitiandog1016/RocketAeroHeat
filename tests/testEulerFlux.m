function testEulerFlux()
%TESTEULERFLUX 用常状态、旋转/反向守恒和独立激波管解检验通量。
g=1.4;normal=[.6,.8];area=2.3;
for method={'HLL','HLLenthalpy'}
    q=state(1.7,.3,-.8,1.2,g); un=.3*normal(1)-.8*normal(2);
    expected=area*[1.7*un,1.7*un*.3+1.2*normal(1), ...
        1.7*un*(-.8)+1.2*normal(2),(q(4)+1.2)*un];
    f=axisymmetricEulerFlux(q,q,area*normal(1),area*normal(2),g,method{1});
    assert(max(abs(reshape(f,1,4)-expected))<1e-12);
    l=state(1.7,.8,-.3,.7,g);r=state(.4,-.1,.7,.3,g);
    f=axisymmetricEulerFlux(l,r,normal(1),normal(2),g,method{1});
    back=axisymmetricEulerFlux(r,l,-normal(1),-normal(2),g,method{1});
    assert(max(abs(f+back),[],'all')<1e-12);
    angle=.37;rot=[cos(angle),-sin(angle);sin(angle),cos(angle)];n=rot*normal.';
    lr=l;rr=r;lr(2:3)=rot*reshape(l(2:3),2,1);rr(2:3)=rot*reshape(r(2:3),2,1);
    fr=axisymmetricEulerFlux(lr,rr,n(1),n(2),g,method{1});
    target=f;target(2:3)=rot*reshape(f(2:3),2,1);
    assert(max(abs(fr-target),[],'all')<1e-12);
    % 反射滑移壁的质量/能量通量严格为零。
    ghost=l;mn=l(2)*normal(1)+l(3)*normal(2);
    ghost(2:3)=reshape(l(2:3),2,1)-2*mn*normal.';
    wall=axisymmetricEulerFlux(l,ghost,normal(1),normal(2),g,method{1});
    assert(abs(wall(1))<1e-12 && abs(wall(4))<1e-12);
end
% Sod 激波管：用独立精确 Riemann 状态核对激波、接触间断和稀疏波。
for method={'HLL','HLLenthalpy'}
    errors=zeros(2,1);
    for j=1:2,errors(j)=sod(200*j,g,method{1});end
    assert(errors(2)<.8*errors(1) && errors(2)<.02);
    fprintf('Sod %s density L1 errors: %.6g -> %.6g\n',method{1},errors);
end
fprintf('Euler flux passed: consistent states, rotation, conservation, slip wall and exact Sod refinement.\n');
end
function q=state(rho,u,v,p,g)
q=reshape([rho,rho*u,rho*v,p/(g-1)+.5*rho*(u*u+v*v)],1,1,4);
end
function errorL1=sod(n,g,method)
% 此处使用物理时间的独立一维积分器；不调用二维局部伪时间推进。
x=((1:n)-.5)'/n;dx=1/n;Q=repmat(state(1,0,0,1,g),n,1,1);
Q(x>.5,:,:)=repmat(state(.125,0,0,.1,g),nnz(x>.5),1,1);
initial=reshape(sum(Q,1)*dx,1,4);time=0;finish=.15;
while time<finish
    rho=Q(:,:,1);u=Q(:,:,2)./rho;p=(g-1)*(Q(:,:,4)-.5*rho.*u.^2);
    dt=min(finish-time,.45*dx/max(abs(u)+sqrt(g*p./rho)));
    ghosts=cat(1,Q(1,:,:),Q,Q(end,:,:));
    f=axisymmetricEulerFlux(ghosts(1:end-1,:,:),ghosts(2:end,:,:),ones(n+1,1),zeros(n+1,1),g,method);
    Q=Q-dt/dx*diff(f,1,1);time=time+dt;
    rho=Q(:,:,1);p=(g-1)*(Q(:,:,4)-.5*Q(:,:,2).^2./rho);
    assert(all(rho>0) && all(p>0));
end
final=reshape(sum(Q,1)*dx,1,4);
assert(max(abs(final([1 4])-initial([1 4])))<1e-9);
assert(abs(final(2)-.9*finish)<1e-9);
% 标准 Sod 问题的独立精确解常数，gamma=1.4。
pstar=.303130178050647;ustar=.92745262004895;cL=sqrt(g);
rhoLstar=pstar^(1/g);rhoRstar=.265573711705307;
cStar=cL*pstar^((g-1)/(2*g));shock=1.75215573203018;
xi=(x-.5)/finish;exact=ones(n,1);
fan=xi>-cL & xi<ustar-cStar;
exact(fan)=(2/(g+1)*(cL-.5*(g-1)*xi(fan))/cL).^(2/(g-1));
exact(xi>=ustar-cStar & xi<ustar)=rhoLstar;
exact(xi>=ustar & xi<shock)=rhoRstar;exact(xi>=shock)=.125;
errorL1=mean(abs(Q(:,:,1)-exact));
end
