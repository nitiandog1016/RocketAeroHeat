function [flux,leftDerivative,rightDerivative]=layeredWallFlux(grid,temperature)
%LAYEREDWALLFLUX 相邻节点热流及其温度导数；正号表示热从内侧流向外侧。
% 材料内部按 Kirchhoff 势差求通量；接触层 Q=(T右-T左)/Rc。
% 所有节点深度按外到内排列；热阻连接没有几何厚度或附加储热。
nf=size(temperature,1); n=grid.nodeCount;
flux=zeros(nf,n-1); leftDerivative=flux; rightDerivative=flux;
for j=1:numel(grid.layers)
    ids=grid.layerNodes{j}; edges=grid.layerLinks{j}; layer=grid.layers{j};
    p=wallMaterialProperties(layer,temperature(:,ids));
    flux(:,edges)=integrateWallMaterialProperty(layer,'thermalConductivity', ...
        temperature(:,ids(1:end-1)),temperature(:,ids(2:end)))/grid.spacing(j);
    leftDerivative(:,edges)=-p.thermalConductivity(:,1:end-1)/grid.spacing(j);
    rightDerivative(:,edges)=p.thermalConductivity(:,2:end)/grid.spacing(j);
end
for j=1:numel(grid.contactResistance)
    if grid.contactResistance(j)>0
        edge=grid.contactLinks(j); g=1/grid.contactResistance(j);
        flux(:,edge)=g*(temperature(:,edge+1)-temperature(:,edge));
        leftDerivative(:,edge)=-g; rightDerivative(:,edge)=g;
    end
end
end
