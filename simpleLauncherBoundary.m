function c=simpleLauncherBoundary(mode)
% 模块用途：按用户所选的低速风洞状态冻结来流、等温壁以及瞬态初值身份。
root=fileparts(mfilename('fullpath'));
if nargin==0,mode='steady';end
mode=validatestring(mode,{'steady','wake'});
folder=fullfile(root,'cases','simple_launcher');
if strcmp(mode,'wake'),file=fullfile(folder,'wake','condition.json');else,file=fullfile(folder,'tunnel.json');end
c=jsondecode(fileread(file));
c.velocityMps=reshape(c.velocityMps,1,3);
end
