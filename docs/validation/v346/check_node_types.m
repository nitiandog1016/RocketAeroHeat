function check_node_types()
% 模块用途：以小型单元连接覆盖五种温界分支，检查共享边与错误归属。
s=(1:25)';d.vertices=[s s.^2 s.^3];d.outer=[8 9 10];
d.tetrahedra=[3 4 2 1;5 6 7 8;11 12 13 14;15 16 17 18;18 19 20 21;22 23 24 25];
d.wall=[1 2 3;5 6 7;11 12 13;11 12 14;15 16 17;18 19 20;22 23 24;22 23 25];
d.wallOwner=[1;2;3;3;4;5;6;6];
faceID=(1:8)';cadFace=faceID;areaM2=ones(8,1);givenWallK=[300;300;300;300;300;315;300;315];
temperatureSource=repmat("外壁",8,1);temperatureSource([4 6 8])="内壁";
wall=table(faceID,cadFace,areaM2,givenWallK,temperatureSource);
r=buildWallNodeConstraints(d,wall);
assert(isequal(r.faces.oppositeNodeID,[4;8;14;13;18;21;25;24]));
assert(isequal(double(r.faces.nodalConstraint),[1;2;3;3;4;1;5;5]));
assert(isequal(r.wallNodes.vertexID(isnan(r.wallNodes.uniqueGivenK)),[22;23]));
assert(all(r.faces.mixedWallNodesInCell(7:8)==2));
assert(string(r.wallNodes.temperatureSource(r.wallNodes.vertexID==11))=="外壁 / 内壁");
d.wallOwner(1)=5;
try
    buildWallNodeConstraints(d,wall);
catch ex
    assert(strcmp(ex.identifier,'RAH:WallNodes:Owner'));return
end
error('错误的单元归属未被拦截。');
end
