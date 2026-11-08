// RocketAeroHeat：仅适配区域边界数据来源；原通量、梯度计算和推进内核不变。
module lmr.bc.ghost_cell_effect.rah_frozen_mapped;

import std.math : isFinite;
import std.conv : to;
import geom;
import nm.number;
import gas;
import lmr.globalconfig;
import lmr.fluidfvcell;
import lmr.flowstate;
import lmr.flowgradients;
import lmr.lsqinterp;
import lmr.bc;
import lmr.bc.ghost_cell_effect.mapped_cell_copy;

// 模块一：接口只保存原生交换使用的数据，不持有donor块/单元引用或时间推进数组。
// 这是已检查状态的内存接口；磁盘文件的身份/轮次验证由调用侧另行绑定。
class FrozenMappedDonor {
    private FlowState state_;
    private LSQInterpGradients recon_;
    private FlowGradients viscous_;
    private Vector3 position_;
    private number volume_, area_, il_, jl_, kl_, min_, max_;
    immutable size_t sourceCellID;

    this(FluidFVCell source, size_t originalID) {
        if (source.myConfig.n_grid_time_levels != 1 || source.myConfig.n_species != 1 ||
            source.myConfig.n_modes != 0 || source.myConfig.turb_model.nturb != 1)
            throw new Exception("RAH_FROZEN_CONFIG");
        // 不裁复数扰动：冻结邻域必须来自外迭代间保存的实数状态。
        foreach (q; [source.fs.gas.rho, source.fs.gas.p, source.fs.gas.T, source.fs.gas.u,
                     source.fs.vel.x, source.fs.vel.y, source.fs.vel.z, source.fs.turb[0],
                     source.volume[0], source.iLength, source.jLength, source.kLength]) {
            if (!isFinite(q.re) || q.im != 0.0) throw new Exception("RAH_FROZEN_NONREAL_STATE");
        }
        if (source.volume[0].re <= 0.0 || source.iLength.re <= 0.0 ||
            source.jLength.re <= 0.0 || source.kLength.re <= 0.0)
            throw new Exception("RAH_FROZEN_GEOMETRY");
        sourceCellID=originalID;
        state_=FlowState(*(source.fs),cast(GasModel)source.myConfig.gmodel);
        recon_=LSQInterpGradients(*(source.gradients));
        viscous_=FlowGradients(*(source.grad));
        position_=source.pos[0];volume_=source.volume[0];area_=source.areaxy[0];
        il_=source.iLength;jl_=source.jLength;kl_=source.kLength;min_=source.L_min;max_=source.L_max;
    }
    void geometry(FluidFVCell target) const {
        target.pos[0]=position_;target.volume[0]=volume_;target.areaxy[0]=area_;
        target.iLength=il_;target.jLength=jl_;target.kLength=kl_;target.L_min=min_;target.L_max=max_;
    }
    void state(FluidFVCell target) const {target.fs.copy_values_from(state_);}
    void reconstruction(FluidFVCell target) const {target.gradients.copy_values_from(recon_);}
    void viscous(FluidFVCell target) const {target.grad.copy_values_from(viscous_);}
    void transport(FluidFVCell target) const {target.fs.mu_t=state_.mu_t;target.fs.k_t=state_.k_t;}
    void shock(FluidFVCell target) const {target.fs.S=state_.S;}
}

// 模块二：继承原生mapped类型，使交换调度、非结构梯度及Jacobian边界分支保持原义。
class RAHFrozenMappedCellCopy : GhostCellMappedCellCopy {
    private FrozenMappedDonor[] frozen;
    private size_t[] faceIDs;
    this(int id, int boundary, FrozenMappedDonor[] donors, size_t[] expectedFaceIDs) {
        Vector3 zero=Vector3(0,0,0),axis=Vector3(0,0,1);
        double[] rotation=[1.,0.,0.,0.,1.,0.,0.,0.,1.];
        super(id,boundary,false,"none",false,false,zero,axis,0.,zero,false,false,rotation);
        if (GlobalConfig.in_mpi_context || blk.grid_type != Grid_t.unstructured_grid ||
            blk.myConfig.n_grid_time_levels!=1 || donors.length==0 || donors.length!=expectedFaceIDs.length)
            throw new Exception("RAH_FROZEN_UNSUPPORTED_ENTRY");
        frozen=donors.dup;faceIDs=expectedFaceIDs.dup;
        foreach (d;frozen) if (d is null) throw new Exception("RAH_FROZEN_MISSING_DONOR");
        desc="RAHFrozenMappedCellCopy";
    }
    override void set_up_cell_mapping() {
        auto boundary=blk.bc[which_boundary];
        if (boundary.faces.length!=frozen.length || ghost_cells.length!=0)
            throw new Exception("RAH_FROZEN_BOUNDARY_COUNT");
        foreach (i,face;boundary.faces) {
            if (face.id!=faceIDs[i]) throw new Exception("RAH_FROZEN_FACE_ID");
            auto ghost=boundary.outsigns[i]==1 ? face.right_cell : face.left_cell;
            if (!ghost.is_ghost) throw new Exception("RAH_FROZEN_NOT_GHOST");
            ghost.is_interior_to_domain=true;
            ghost_cells~=ghost;
        }
        // mapped_cells故意保持空；不存在整个供体块的驻留引用。
    }
    // 模块三：原生三个交换阶段的顺序保持，单进程只需phase2拷贝冻结数据。
    override void exchange_geometry_phase0() {}
    override void exchange_geometry_phase1() {}
    override void exchange_geometry_phase2() {foreach(i,g;ghost_cells) frozen[i].geometry(g);}
    override void exchange_flowstate_phase0(double t,int gtl,int ftl) {}
    override void exchange_flowstate_phase1(double t,int gtl,int ftl) {}
    override void exchange_flowstate_phase2(double t,int gtl,int ftl) {foreach(i,g;ghost_cells) frozen[i].state(g);}
    override void exchange_convective_gradient_phase0(double t,int gtl,int ftl) {}
    override void exchange_convective_gradient_phase1(double t,int gtl,int ftl) {}
    override void exchange_convective_gradient_phase2(double t,int gtl,int ftl) {foreach(i,g;ghost_cells) frozen[i].reconstruction(g);}
    override void exchange_viscous_gradient_phase0(double t,int gtl,int ftl) {}
    override void exchange_viscous_gradient_phase1(double t,int gtl,int ftl) {}
    override void exchange_viscous_gradient_phase2(double t,int gtl,int ftl) {foreach(i,g;ghost_cells) frozen[i].viscous(g);}
    override void exchange_turbulent_transprops_phase0() {}
    override void exchange_turbulent_transprops_phase1() {}
    override void exchange_turbulent_transprops_phase2() {foreach(i,g;ghost_cells) frozen[i].transport(g);}
    override void exchange_shock_phase0(double t,int gtl,int ftl) {}
    override void exchange_shock_phase1(double t,int gtl,int ftl) {}
    override void exchange_shock_phase2(double t,int gtl,int ftl) {foreach(i,g;ghost_cells) frozen[i].shock(g);}
}
