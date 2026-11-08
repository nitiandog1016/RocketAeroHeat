// 模块职责：使用完整原生残差流程对照冻结边界；没有推进/牛顿迭代入口。
module rah_regional_stencil_probe;
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
size_t[] donorMotherIDs;
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
        donors~=c;interfaceIDs~=face.id;donorMotherIDs~=to!size_t(row["motherCellID"].integer);
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



// 模块四：只读指定位置的实际局部输入与切向，必须在原始量虚部清理之前观察。
double[] direction;
void need(bool ok,string msg){if(!ok)throw new Exception(msg);}
void row(ref File f,int regime,int sample,int cell,string kind,int index,int component,number q){
    need(isFinite(q.re)&&isFinite(q.im),"RAH_OBSERVE_NONFINITE");
    f.writefln("%s,%s,%s,%s,%s,%s,%.17g,%.17g",regime,sample,cell,kind,index,component,q.re,q.im/1e-30);
}
void state(ref File f,int regime,int sample,int cell,string kind,int index,FlowState* fs){
    number[] q=[fs.gas.T,fs.vel.x,fs.vel.y,fs.vel.z,fs.turb[0],fs.gas.rho,fs.gas.u,fs.gas.p,fs.gas.mu,fs.gas.k];
    foreach(j,v;q)row(f,regime,sample,cell,kind,index,cast(int)j,v);
}
void readDirection(int regime){
    import std.string : split;
    direction.length=6*block.cells.length;int count=0;
    foreach(line;File("directions.csv").byLineCopy()){
        auto columns=line.split(",");if(columns[0]=="regime"||to!int(columns[0])!=regime)continue;
        auto id=to!size_t(columns[1]);need(id<block.cells.length,"RAH_DIRECTION_ID");
        foreach(k;0..6)direction[6*id+k]=to!double(columns[11+k]);count++;
    }
    need(count==block.cells.length,"RAH_DIRECTION_COUNT");
}
void observe(ref File output,int regime,int sample,int primary){
    int[] ids=[primary];
    foreach(f;block.cells[primary].iface){
        auto other=(f.left_cell is block.cells[primary])?f.right_cell:f.left_cell;
        if(other.id>=0 && other.id<block.cells.length && other is block.cells[other.id]){
            bool found=false;foreach(id;ids)found=found||id==other.id;if(!found)ids~=other.id;
        }
    }
    foreach(id;ids){auto c=block.cells[id];auto ws=c.ws_grad;
        need(ws.n==6 && ws.loop_init==1 && !ws.compute_about_mid,"RAH_OBSERVE_CLOUD");
        foreach(i,fs;c.cloud_fs){
            state(output,regime,sample,id,"cloud",cast(int)i,fs);
            auto p=c.cloud_pos[i];number[3] xyz=[p.x,p.y,p.z];
            foreach(j,q;xyz)row(output,regime,sample,id,"position",cast(int)i,cast(int)j,q);
            if(i>0){number[3] w=[ws.wx[i],ws.wy[i],ws.wz[i]];
                foreach(j,q;w)row(output,regime,sample,id,"weight",cast(int)i,cast(int)j,q);}
        }
        foreach(j;0..3){
            row(output,regime,sample,id,"gradient",0,j,c.grad.T[j]);
            foreach(k;0..3)row(output,regime,sample,id,"gradient",k+1,j,c.grad.vel[k][j]);
            row(output,regime,sample,id,"gradient",4,j,c.grad.turb[0][j]);
        }
        foreach(j;0..6){
            row(output,regime,sample,id,"residual",0,j,c.dUdt[1][j]);
            row(output,regime,sample,id,"source",0,j,c.Q[j]);
            row(output,regime,sample,id,"U0",0,j,c.U[0][j]);
            row(output,regime,sample,id,"U1",0,j,c.U[1][j]);
        }
        row(output,regime,sample,id,"geometry",0,0,c.volume[0]);
        row(output,regime,sample,id,"geometry",0,1,c.dwall);
        foreach(i,f;c.iface){
            row(output,regime,sample,id,"faceID",cast(int)i,0,to!number(f.id));
            row(output,regime,sample,id,"outsign",cast(int)i,0,to!number(c.outsign[i]));
            foreach(j;0..6)row(output,regime,sample,id,"flux",cast(int)i,j,f.F[j]*f.area[0]);
            state(output,regime,sample,id,"left",cast(int)i,f.left_cell.fs);
            state(output,regime,sample,id,"right",cast(int)i,f.right_cell.fs);
        }
    }
    output.flush();
}
int main(string[] args){
    GlobalConfig.in_mpi_context=false;GlobalConfig.is_master_task=true;
    initNewtonKrylovSimulation(0,1,1,"60");
    need(localFluidBlocks.length==1,"RAH_OBSERVE_BLOCK");block=localFluidBlocks[0];
    need(block.myConfig.cqi.n==6 && block.cells.length==416,"RAH_OBSERVE_LAYOUT");
    readNewtonKrylovConfig();activePhase=nkPhases[0];SimState.time=0;setupDonors();
    need(activePhase.jacobianInterpolationOrder==1 && activePhase.residualInterpolationOrder==2 &&
        !activePhase.frozenLimiterForJacobian && !activePhase.frozenLimiterForResidual,"RAH_OBSERVE_SETTINGS");
    auto output=File("stencil.csv","w");output.writeln("regime,sample,cell,kind,index,component,re,imScaled");
    int[3] primary=[23,22,5];
    foreach(regime;0..3){
        FrozenMappedDonor[] snapshots;
        foreach(i,c;donors){fill(*(c.fs),c.pos[0],regime);donorGradients(c);snapshots~=new FrozenMappedDonor(c,donorMotherIDs[i]);}
        frozen=new RAHFrozenMappedCellCopy(block.id,bindex,snapshots,interfaceIDs);frozen.set_up_cell_mapping();frozen.exchange_geometry_phase2();
        block.bc[bindex].preReconAction=[cast(GhostCellEffect)frozen];readDirection(regime);
        foreach(sample;0..4){
            reset(regime);setJacobianEvalSettings(false);block.clear_fluxes_of_conserved_quantities();
            foreach(c;block.cells)c.clear_source_vector();
            foreach(i,c;block.cells){
                c.U[1][]=c.U[0][];
                foreach(k;0..6){
                    if(sample==3)c.U[1][k]+=complex(0.,1e-30*direction[6*i+k]);
                    else c.U[1][k]+=(sample==1?1e-4:sample==2?-1e-4:0.)*direction[6*i+k];
                }
                need(c.decode_conserved(0,1,0)==0,"RAH_OBSERVE_DECODE");
            }
            evalResidual(1,0,0);observe(output,regime,sample,primary[regime]);
            foreach(c;block.cells)c.fs.clear_imaginary_components();
            foreach(bc;block.bc)foreach(c;bc.ghostcells)c.fs.clear_imaginary_components();
            foreach(f;block.faces)f.fs.clear_imaginary_components();setResidualEvalSettings();
        }
    }
    output.close();writeln("RAH_STENCIL_OBSERVATION_COMPLETE");return 0;
}
