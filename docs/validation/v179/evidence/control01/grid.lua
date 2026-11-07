-- Run14物理外壁与官方两块外域；坐标为本验证算例局部轴向+X。
config.dimensions=2
config.axisymmetric=true
config.grid_format="gziptext"
local L1,L2,R=0.1017,0.220,0.0325
local rise=math.tan(math.rad(30))*(L2-L1)
local a0=Vector3:new{x=0,y=R};local a1=Vector3:new{x=0,y=R+0.005}
local b0=Vector3:new{x=L1,y=R};local b1=Vector3:new{x=L1-0.005,y=R+0.020}
local c0=Vector3:new{x=L2,y=R+rise};local c1=Vector3:new{x=L2,y=R+rise+0.025}
local wallCluster=RobertsFunction:new{end0=true,end1=false,beta=1.01}
local noseCluster=RobertsFunction:new{end0=true,end1=false,beta=1.2}
registerFluidGrid{grid=StructuredGrid:new{
  psurface=CoonsPatch:new{p00=a0,p10=b0,p11=b1,p01=a1},niv=161,njv=97,
  cfList={north=noseCluster,south=noseCluster,east=wallCluster,west=wallCluster}},
  fsTag="inflow",bcTags={north="inflow",west="inflow",south="wall"}}
registerFluidGrid{grid=StructuredGrid:new{
  psurface=CoonsPatch:new{p00=b0,p10=c0,p11=c1,p01=b1},niv=161,njv=97,
  cfList={east=wallCluster,west=wallCluster}},
  fsTag="initial",bcTags={north="inflow",east="outflow",south="wall"}}
identifyGridConnections()
