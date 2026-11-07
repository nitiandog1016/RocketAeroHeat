-- 原始网格只读引用；不修改壁面、层数或尾缘。
config.dimensions = 3
config.grid_format = "rawbinary"
registerFluidGrid{grid=UnstructuredGrid:new{filename="/mnt/e/RocketAeroHeat/outputs/v178_validation/fixture03/mixed_original.su2",
    fmt="su2text", scale=1.0, expect_gmsh_order_for_wedges=false}, fsTag="initial"}
