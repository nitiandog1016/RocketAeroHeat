// 模块职责：使用完整原生残差流程对照冻结边界；没有推进/牛顿迭代入口。
module rah_frozen_mapped_probe;
import std.stdio;
import std.file;
import std.json;
import std.conv;
import std.math;
import std.algorithm : sort, max;
import std.array : array;
import std.string : join;
import geom;
import gas;
import nm.number;
import ntypes.complex;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.newtonkrylovsolver;
import lmr.commands.prepsim;
import lmr.coredata;
import lmr.fluidfvcell;
import lmr.fluidblock;
import lmr.flowstate;
import lmr.flowgradients;
import lmr.lsqinterp;
import lmr.bc;
import lmr.bc.ghost_cell_effect.rah_frozen_mapped;
import lmr.bc.ghost_cell_effect.mapped_cell_copy;
import lmr.bc.boundary_interface_effect;

FluidBlock block;
Vector3 origin;
double[3] thermalGradient=[3.,-5.,7.];
double[3][3] velocityGradient=[[20.,-10.,3.],[-5.,8.,12.],[9.,4.,-7.]];
double[3] saGradient=[1e-5,-2e-5,1e-5];
FluidCellData donorData;
FluidFVCell[] donors;
GhostCellMappedCellCopy resident;
RAHFrozenMappedCellCopy frozen;
int bindex;
size_t[] interfaceIDs;

Vector3 vec(JSONValue row) {return Vector3(row[0].floating,row[1].floating,row[2].floating);}
string key(size_t[] ids) {auto sorted=ids.dup;sort(sorted);string[] s;foreach(i;sorted)s~=to!string(i);return s.join(":");}
double coordinate(ref const(Vector3) p,size_t i) {return i==0?p.x.re:i==1?p.y.re:p.z.re;}

void fill(ref FlowState fs,ref const(Vector3) pos,int regime,number perturb=0) {
    // 模块一：可微、合法的非均匀状态；仅作离散算子控制，不声称为物理解。
    double[3] x=[pos.x.re-origin.x.re,pos.y.re-origin.y.re,pos.z.re-origin.z.re];
    fs.gas.p=2500.+perturb;
    fs.gas.T=300.+perturb*.2;
    foreach(i;0..3)fs.gas.T+=thermalGradient[i]*x[i];
    number[3] u=[to!number(regime==0?-1000.:regime==1?1000.:0.),to!number(0.),to!number(regime==2?0.:200.)];
    foreach(i;0..3)foreach(j;0..3)u[i]+=velocityGradient[i][j]*x[j];
    fs.vel=Vector3(u[0]+perturb*.3,u[1]-perturb*.4,u[2]+perturb*.2);
    fs.turb[0]=1e-5+perturb*1e-7;fs.turb[1]=0;
    foreach(i;0..3)fs.turb[0]+=saGradient[i]*x[i];
    fs.mu_t=0;fs.k_t=0;fs.S=0;
    auto gm=block.myConfig.gmodel;
    gm.update_thermo_from_pT(fs.gas);gm.update_sound_speed(fs.gas);gm.update_trans_coeffs(fs.gas);
}

void donorGradients(FluidFVCell c) {
    auto g=c.gradients;
    g.velx=Vector3(20,-10,3);g.vely=Vector3(-5,8,12);g.velz=Vector3(9,4,-7);
    g.T=Vector3(3,-5,7);g.p=Vector3(0,0,0);
    g.rho=g.T*(-c.fs.gas.rho/c.fs.gas.T);g.u=g.T*block.myConfig.gmodel.Cv(c.fs.gas);
    g.turb[0]=Vector3(1e-5,-2e-5,1e-5);g.turb[1]=Vector3(0,0,0);g.rho_s[0]=g.rho;
    // 全部使用量都有有界极值；未使用第二湍流槽的phi明确为-1。
    g.velxPhi=1;g.velyPhi=1;g.velzPhi=1;g.TPhi=1;g.pPhi=1;g.rhoPhi=1;g.uPhi=1;
    g.velxMin=-2e3;g.velxMax=2e3;g.velyMin=-2e3;g.velyMax=2e3;g.velzMin=-2e3;g.velzMax=2e3;
    g.TMin=200;g.TMax=500;g.pMin=1000;g.pMax=10000;g.rhoMin=.001;g.rhoMax=1;g.uMin=-1e6;g.uMax=1e6;
    g.turbPhi[0]=1;g.turbPhi[1]=-1;g.turbMin[0]=0;g.turbMax[0]=.01;g.turbMin[1]=0;g.turbMax[1]=0;
    g.rho_sPhi[0]=-1;g.rho_sMin[0]=.001;g.rho_sMax[0]=1;
    foreach(i;0..3){c.grad.T[i]=thermalGradient[i];c.grad.turb[0][i]=saGradient[i];c.grad.turb[1][i]=0;c.grad.massf[0][i]=0;
        foreach(j;0..3)c.grad.vel[i][j]=velocityGradient[i][j];}
}

void setupDonors() {
    auto data=parseJSON(readText("donors.json"));origin=vec(data["origin"]);
    foreach(i,b;block.bc)if(b.label=="exchange")bindex=to!int(i);
    auto boundary=block.bc[bindex];auto n=boundary.faces.length;
    assert(n==data["donors"].array.length && n>0);
    donorData.allocate_minimal(n,0,block.myConfig);
    foreach(i;0..n)donorData.lsqgradients~=LSQInterpGradients(1,0,1);
    JSONValue[string] records;
    foreach(row;data["donors"].array){size_t[] ids;foreach(id;row["faceVertices"].array)ids~=to!size_t(id.integer);records[key(ids)]=row;}
    Vector3 zero=Vector3(0,0,0),axis=Vector3(0,0,1);double[] rotation=[1.,0.,0.,0.,1.,0.,0.,0.,1.];
    resident=new GhostCellMappedCellCopy(block.id,bindex,false,"none",false,false,zero,axis,0.,zero,false,false,rotation);
    foreach(i,face;boundary.faces){
        size_t[] ids;foreach(v;face.vtx)ids~=v.id;auto row=records[key(ids)];Vector3[] p;
        foreach(x;row["vertices"].array)p~=vec(x);
        auto c=new FluidFVCell(block.myConfig,&donorData,to!int(i));
        wedge_properties(p[0],p[1],p[2],p[3],p[4],p[5],false,c.pos[0],c.volume[0],c.L_min);
        c.iLength=c.L_min;c.jLength=c.L_min;c.kLength=c.L_min;c.L_max=c.L_min;c.areaxy[0]=0;
        auto ghost=boundary.outsigns[i]==1?face.right_cell:face.left_cell;
        ghost.is_interior_to_domain=true;resident.ghost_cells~=ghost;resident.mapped_cells~=c;
        donors~=c;interfaceIDs~=face.id;
    }
    boundary.type="exchange_using_mapped_cells";
    boundary.convective_flux_computed_in_bc=false;
    boundary.preSpatialDerivActionAtBndryFaces=[new BIE_UpdateThermoTransCoeffs(block.id,bindex,block.myConfig.thermo_interpolator)];
    boundary.preReconAction=[resident];resident.exchange_geometry_phase2();
    block.compute_least_squares_setup(0);
}

void reset(int regime,number perturb=0) {
    foreach(c;block.cells){fill(*(c.fs),c.pos[0],regime,perturb);c.encode_conserved(0,0,0);}
    foreach(f;block.faces)fill(*(f.fs),f.pos,regime);
}

number[] sample(bool useFrozen,int regime,number perturb=0) {
    block.bc[bindex].preReconAction=useFrozen?[cast(GhostCellEffect)frozen]:[cast(GhostCellEffect)resident];
    reset(regime,perturb);
    evalResidualWorker(0);
    number[] values;
    foreach(c;block.cells){foreach(q;c.dUdt[0])values~=q;foreach(q;c.grad.T)values~=q;
        foreach(row;c.grad.vel)foreach(q;row)values~=q;foreach(q;c.grad.turb[0])values~=q;}
    foreach(f;block.faces)foreach(q;f.F)values~=q*f.area[0];
    foreach(q;values)if(!isFinite(q.re)||!isFinite(q.im))throw new Exception("RAH_NONFINITE_OPERATOR");
    return values;
}

double scaledDifference(number[] a,number[] b,bool imaginary=false,double epsilon=1e-30) {
    assert(a.length==b.length);double error=0;
    foreach(i;0..a.length){double av=imaginary?a[i].im/epsilon:a[i].re,bv=imaginary?b[i].im/epsilon:b[i].re;
        error=max(error,fabs(av-bv)/max(1.,max(fabs(av),fabs(bv))));}return error;
}

int main(string[] args) {
    GlobalConfig.in_mpi_context=false;GlobalConfig.is_master_task=true;
    if(args.length>1 && args[1]=="prepare") return lmr.commands.prepsim.main_([args[0],"--job=job.lua"]);
    initNewtonKrylovSimulation(0,1,1,"60");
    assert(localFluidBlocks.length==1);block=localFluidBlocks[0];assert(block.myConfig.cqi.n==6);
    readNewtonKrylovConfig();activePhase=nkPhases[0];SimState.time=0;
    setupDonors();
    // 模块二：真实原棱柱face点云的原生weighted-QR，一般三维线性场及法向梯度。
    reset(0);block.flow_property_spatial_derivatives(0);double gradientError=0,normalError=0;
    auto gradientFile=File("gradient.csv","w");
    gradientFile.writeln("cell,x,y,z,nx,ny,nz,Tx,Ty,Tz,ux,uy,uz,vx,vy,vz,wx,wy,wz");
    foreach(c;block.cells){foreach(j;0..3){gradientError=max(gradientError,fabs(c.grad.T[j].re-thermalGradient[j])/max(1.,fabs(thermalGradient[j])));
        foreach(i;0..3)gradientError=max(gradientError,fabs(c.grad.vel[i][j].re-velocityGradient[i][j])/max(1.,fabs(velocityGradient[i][j])));}
        auto n=c.iface[0].n;double expected=thermalGradient[0]*n.x.re+thermalGradient[1]*n.y.re+thermalGradient[2]*n.z.re;
        double actual=c.grad.T[0].re*n.x.re+c.grad.T[1].re*n.y.re+c.grad.T[2].re*n.z.re;
        normalError=max(normalError,fabs(actual-expected)/max(1.,fabs(expected)));
        gradientFile.writef("%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",c.id,c.pos[0].x.re,c.pos[0].y.re,c.pos[0].z.re,n.x.re,n.y.re,n.z.re);
        foreach(q;c.grad.T)gradientFile.writef(",%.17g",q.re);
        foreach(row;c.grad.vel)foreach(q;row)gradientFile.writef(",%.17g",q.re);
        gradientFile.writeln();}
    gradientFile.close();
    auto output=File("operator_results.jsonl","w");
    output.writefln("{\"stage\":\"gradient\",\"cells\":%s,\"faces\":%s,\"donors\":%s,\"error\":%.17g,\"normalError\":%.17g}",
        block.cells.length,block.faces.length,donors.length,gradientError,normalError);output.flush();
    if(gradientError>1e-7 || normalError>1e-7)throw new Exception("RAH_GRADIENT_GATE");
    // 模块三：完整原生面通量/六残差/梯度与仅受体方向导数，不做流动发展。
    // 分量输出用于独立核对；全部预定方向完成后统一判定，门槛保持不变。
    auto raw=File("samples.csv","w");
    raw.writeln("regime,index,kind,id,component,resident,frozen,residentDerivative,frozenDerivative,plus,minus");
    auto geometry=File("geometry.csv","w");
    geometry.writeln("kind,id,volumeOrArea,x,y,z,nx,ny,nz,wallDistance");
    foreach(c;block.cells)geometry.writefln("cell,%s,%.17g,%.17g,%.17g,%.17g,0,0,0,%.17g",c.id,c.volume[0].re,c.pos[0].x.re,c.pos[0].y.re,c.pos[0].z.re,c.dwall.re);
    foreach(f;block.faces)geometry.writefln("face,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,0",f.id,f.area[0].re,f.pos.x.re,f.pos.y.re,f.pos.z.re,f.n.x.re,f.n.y.re,f.n.z.re);
    geometry.close();bool accepted=true;
    foreach(regime;0..3){
        FrozenMappedDonor[] snapshots;
        foreach(i,c;donors){fill(*(c.fs),c.pos[0],regime);donorGradients(c);snapshots~=new FrozenMappedDonor(c,i);}
        frozen=new RAHFrozenMappedCellCopy(block.id,bindex,snapshots,interfaceIDs);frozen.set_up_cell_mapping();
        frozen.exchange_geometry_phase2();assert(frozen.mapped_cells.length==0);
        auto ref0=sample(false,regime),snap0=sample(true,regime);
        auto refJ=sample(false,regime,complex(0.,1e-30)),snapJ=sample(true,regime,complex(0.,1e-30));
        auto plus=sample(true,regime,to!number(1e-4)),minus=sample(true,regime,to!number(-1e-4));
        double error=scaledDifference(ref0,snap0),deriv=scaledDifference(refJ,snapJ,true),fd=0;
        foreach(i;0..snapJ.length){double av=snapJ[i].im/1e-30,bv=(plus[i].re-minus[i].re)/2e-4;
            fd=max(fd,fabs(av-bv)/max(1.,max(fabs(av),fabs(bv))));
            string kind;size_t id,component;
            if(i<21*block.cells.length){auto j=i%21;id=block.cells[i/21].id;
                kind=j<6?"residual":j<9?"temperatureGradient":j<18?"velocityGradient":"saGradient";
                component=j<6?j:j<9?j-6:j<18?j-9:j-18;
            } else {auto j=i-21*block.cells.length;id=block.faces[j/6].id;component=j%6;kind="integratedFlux";}
            raw.writefln("%s,%s,%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",regime,i,kind,id,component,
                ref0[i].re,snap0[i].re,refJ[i].im/1e-30,snapJ[i].im/1e-30,plus[i].re,minus[i].re);}
        raw.flush();
        output.writefln("{\"stage\":\"operator\",\"regime\":%s,\"components\":%s,\"error\":%.17g,\"derivativeError\":%.17g,\"finiteDifferenceError\":%.17g,\"residentDonorReferences\":%s}",
            regime,snap0.length,error,deriv,fd,frozen.mapped_cells.length);output.flush();
        if(error>1e-11 || deriv>1e-11 || fd>1e-5)accepted=false;
    }
    raw.close();output.close();
    if(!accepted)throw new Exception("RAH_OPERATOR_GATE");
    writeln("RAH_NATIVE_FROZEN_OPERATOR_PASS");return 0;
}
