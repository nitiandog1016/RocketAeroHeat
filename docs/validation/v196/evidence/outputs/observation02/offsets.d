import lmr.fluidfvcell; import lmr.flowstate; import lmr.fvinterface; import gas;
pragma(msg, "RAH-OFFSET lengths=", FluidFVCell.lengths.offsetof);
pragma(msg, "RAH-OFFSET L_min=", FluidFVCell.L_min.offsetof);
pragma(msg, "RAH-OFFSET fs=", FluidFVCell.fs.offsetof);
pragma(msg, "RAH-OFFSET iface=", FluidFVCell.iface.offsetof);
pragma(msg, "RAH-OFFSET vel=", FlowState.vel.offsetof);
pragma(msg, "RAH-OFFSET gas=", FlowState.gas.offsetof);
pragma(msg, "RAH-OFFSET a=", GasState.a.offsetof);
pragma(msg, "RAH-OFFSET normal=", FVInterface.n.offsetof);