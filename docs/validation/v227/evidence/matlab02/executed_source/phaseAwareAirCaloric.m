function c=phaseAwareAirCaloric(T)
%PHASEAWAREAIRCALORIC 新波系的数学量热积分，不是经过压力相态域审核的气体状态。
% 模块用途：从已绑定的Cp表返回h、s0及声速，用于守恒方程的有界寻根。
% 仅60—2000K；最终物理p/T必须另由phaseAwareAirState审核，不延拓表格。
persistent x cp b h s R digest
if isempty(x)
    file=fullfile(fileparts(mfilename('fullpath')),'data','thermo','phase_aware_dry_air_v1.json');
    fid=fopen(file,'rb');assert(fid>=0);clean=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
    digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    assert(strcmp(digest,'7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'), ...
        'RAH:PhaseWave:Table','Cp源表身份错误。');
    d=jsondecode(native2unicode(bytes.','UTF-8'));x=d.temperatureK(:);cp=d.cpJkgK(:);R=d.gasConstantJkgK;
    b=diff(cp)./diff(x);h=[0;cumsum((cp(1:end-1)+cp(2:end)).*diff(x)/2)];
    s=[0;cumsum((cp(1:end-1)-b.*x(1:end-1)).*log(x(2:end)./x(1:end-1))+b.*diff(x))];
    h=h-h(x==298.15);s=s-s(x==298.15);
end
validateattributes(T,{'numeric'},{'real','finite','nonempty','>=',60,'<=',2000});
i=discretize(T,x);xi=reshape(x(i),size(T));ci=reshape(cp(i),size(T));bi=reshape(b(i),size(T));dt=T-xi;
C=ci+bi.*dt;
c=struct('cp',C,'h',reshape(h(i),size(T))+ci.*dt+bi.*dt.^2/2, ...
    's0',reshape(s(i),size(T))+(ci-bi.*xi).*log(T./xi)+bi.*dt, ...
    'a',sqrt(C./(C-R).*R.*T),'R',R,'sourceSHA256',digest);
end
