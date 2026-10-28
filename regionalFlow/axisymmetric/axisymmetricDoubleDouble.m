function out=axisymmetricDoubleDouble(operation,a,b)
%AXISYMMETRICDOUBLEDOUBLE 用两个double保存约106位有效二进制的显式算术。
% 模块一：TwoSum及Dekker分裂记录加/乘形成的舍入余项，不使用符号工具箱。
% 模块二：除法/平方根以余量校正；所有数组按MATLAB逐元素规则广播。
% 模块三：只供明确选择的通量内核使用；不修改Q、不裁值、不静默降回普通double。
% 输入为数值或struct(hi,lo)，输出同结构；value操作返回最终舍入后的double。
if nargin<3,b=[];end
a=pair(a);if ~isempty(b),b=pair(b);end
switch operation
    case 'from',out=a;
    case 'value',out=a.hi+a.lo;
    case 'add'
        [s,e]=twoSum(a.hi,b.hi);[t,f]=twoSum(a.lo,b.lo);
        [s,e]=twoSum(s,e+t);[h,l]=twoSum(s,e+f);out=struct('hi',h,'lo',l);
    case 'sub',b.hi=-b.hi;b.lo=-b.lo;out=axisymmetricDoubleDouble('add',a,b);
    case 'mul'
        p=a.hi.*b.hi;split=134217729;
        ca=split*a.hi;cb=split*b.hi;
        assert(all(isfinite(ca)&isfinite(cb),'all'),'AxisDD:Range','Dekker分裂溢出，未退回低精度。');
        ah=ca-(ca-a.hi);al=a.hi-ah;bh=cb-(cb-b.hi);bl=b.hi-bh;
        error=((ah.*bh-p)+ah.*bl+al.*bh)+al.*bl;
        tail=error+(a.hi.*b.lo+a.lo.*b.hi)+a.lo.*b.lo;
        [h,l]=twoSum(p,tail);out=struct('hi',h,'lo',l);
    case 'div'
        assert(all(b.hi~=0,'all'),'AxisDD:Division','除数为零。');
        q=a.hi./b.hi;out=pair(q);
        for k=1:2
            r=axisymmetricDoubleDouble('sub',a,axisymmetricDoubleDouble('mul',b,out));
            correction=(r.hi+r.lo)./b.hi;out=axisymmetricDoubleDouble('add',out,correction);
        end
    case 'sqrt'
        assert(all(a.hi>=0,'all'),'AxisDD:Sqrt','平方根输入为负。');
        q=sqrt(a.hi+a.lo);out=pair(q);positive=q>0;
        % 全部真实Euler声学输入严格正；零数组也明确支持而不以epsilon代替。
        for k=1:2
            r=axisymmetricDoubleDouble('sub',a,axisymmetricDoubleDouble('mul',out,out));
            correction=zeros(size(q));correction(positive)=(r.hi(positive)+r.lo(positive))./(2*out.hi(positive));
            out=axisymmetricDoubleDouble('add',out,correction);
        end
    otherwise,error('AxisDD:Operation','未知双double算术操作。');
end
if isstruct(out),assert(all(isfinite(out.hi)&isfinite(out.lo),'all'),'AxisDD:Arithmetic','双double结果非有限。');end
end
function out=pair(a)
if isstruct(a)
    assert(isscalar(a)&&isequal(sort(fieldnames(a)),{'hi';'lo'})&&isequal(size(a.hi),size(a.lo)), ...
        'AxisDD:Input','双double结构无效。');out=a;
else
    assert(isnumeric(a)&&isreal(a)&&all(isfinite(a),'all'),'AxisDD:Input','需要有限实数。');
    out=struct('hi',double(a),'lo',zeros(size(a)));
end
end
function [s,e]=twoSum(a,b)
s=a+b;back=s-a;e=(a-(s-back))+(b-back);
end
