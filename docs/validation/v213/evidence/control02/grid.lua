-- 模块二：原生只读取已有棱柱控制，保持楔体节点顺序。
config.dimensions=3
config.grid_format="rawbinary"
registerFluidGrid{grid=UnstructuredGrid:new{filename="mesh.su2",fmt="su2text",scale=1.0,expect_gmsh_order_for_wedges=false},fsTag="initial"}
