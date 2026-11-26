function c=simpleLauncherBoundary()
% 模块用途：从项目内固定低速风洞工况读取来流、物性和等温壁；不继承飞行入口默认值。
root=fileparts(mfilename('fullpath'));
c=jsondecode(fileread(fullfile(root,'cases','simple_launcher','tunnel.json')));
c.velocityMps=reshape(c.velocityMps,1,3);
end
