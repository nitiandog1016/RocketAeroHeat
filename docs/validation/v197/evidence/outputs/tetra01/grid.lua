-- 新挤出控制网格，非真实20度箭体。
config.dimensions = 3
config.grid_format = "rawbinary"
registerFluidGrid{grid=UnstructuredGrid:new{filename="mesh.su2",fmt="su2text",
 scale=1.0,expect_gmsh_order_for_wedges=false},fsTag="initial"}
