// 模块职责：使用完整原生残差流程对照冻结边界；没有推进/牛顿迭代入口。
module rah_native_nk_probe;
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


// 模块四：独立有向量控制；正式Jv与实差分共用原生守恒变量解码和一阶设置。
double[] direction;
void need(bool ok,string msg) {if(!ok)throw new Exception(msg);}

void choose(bool snapshot) {
    block.bc[bindex].preReconAction=snapshot?[cast(GhostCellEffect)frozen]:[cast(GhostCellEffect)resident];
}

number[] capture() {
    number[] values;
    foreach(c;block.cells){foreach(q;c.dUdt[1])values~=q;foreach(q;c.grad.T)values~=q;
        foreach(row;c.grad.vel)foreach(q;row)values~=q;foreach(q;c.grad.turb[0])values~=q;}
    foreach(f;block.faces)foreach(q;f.F)values~=q*f.area[0];
    foreach(q;values)need(isFinite(q.re)&&isFinite(q.im),"RAH_NK_NONFINITE");
    return values;
}

void makeDirection(int regime,ref File file) {
    direction.length=6*block.cells.length;size_t j=0;
    foreach(c;block.cells){
        fill(*(c.fs),c.pos[0],regime,complex(0.,1e-30));c.encode_conserved(0,0,0);
        file.writef("%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",regime,c.id,
            c.fs.gas.p.re,c.fs.gas.T.re,c.fs.gas.rho.re,c.fs.gas.u.re,c.fs.vel.x.re,c.fs.vel.y.re,c.fs.vel.z.re,
            c.fs.turb[0].re,block.myConfig.gmodel.Cv(c.fs.gas).re);
        foreach(q;c.U[0]){direction[j++]=q.im/1e-30;file.writef(",%.17g",q.im/1e-30);}
        file.writeln();
    }
    reset(regime);file.flush();
}

number[] realSample(bool snapshot,int regime,double step) {
    choose(snapshot);reset(regime);setJacobianEvalSettings(false);
    need(block.myConfig.interpolation_order==1,"RAH_NK_ORDER");size_t j=0;
    foreach(c;block.cells){
        c.U[1][]=c.U[0][];
        foreach(k;0..6)c.U[1][k]+=step*direction[j++];
        need(c.decode_conserved(0,1,0)==0,"RAH_NK_DECODE");
    }
    evalResidual(1,0,0);auto result=capture();setResidualEvalSettings();return result;
}

number[] nativeJv(bool snapshot,int regime,ref double zedError) {
    choose(snapshot);reset(regime);block.zed[]=direction[];
    // 使用正式分派入口，包括原生时间层、设置切换、清理与恢复；不替代Jv。
    evalJacobianVectorProduct(1e-30,0,0,false);
    auto result=capture();
    foreach(i,c;block.cells)foreach(k;0..6){
        double a=block.zed[6*i+k],b=result[21*i+k].im/1e-30;
        zedError=max(zedError,fabs(a-b)/max(1.,max(fabs(a),fabs(b))));
    }
    need(block.myConfig.interpolation_order==2,"RAH_NK_RESTORE_ORDER");
    return result;
}

int main(string[] args) {
    GlobalConfig.in_mpi_context=false;GlobalConfig.is_master_task=true;
    initNewtonKrylovSimulation(0,1,1,"60");
    need(localFluidBlocks.length==1,"RAH_NK_BLOCK_COUNT");block=localFluidBlocks[0];
    need(block.myConfig.cqi.n==6 && block.cells.length==416,"RAH_NK_LAYOUT");
    readNewtonKrylovConfig();activePhase=nkPhases[0];SimState.time=0;
    need(activePhase.jacobianInterpolationOrder==1 && activePhase.residualInterpolationOrder==2 &&
        !activePhase.frozenLimiterForJacobian && !activePhase.frozenLimiterForResidual &&
        !nkCfg.useRealValuedFrechetDerivative && nkCfg.numberOfPhases==1,"RAH_NK_CONFIG");
    block.zed.length=6*block.cells.length;setupDonors();
    auto directions=File("directions.csv","w");
    directions.writeln("regime,cell,p,T,rho,u,vx,vy,vz,nuhat,Cv,drho,dmx,dmy,dmz,dE,dsa");
    auto output=File("nk_results.jsonl","w");auto raw=File("nk_samples.csv","w");
    raw.writeln("regime,index,kind,id,component,resident,frozen,residentDerivative,frozenDerivative,residentPlus,residentMinus,frozenPlus,frozenMinus");
    bool accepted=true;
    foreach(regime;0..3){
        FrozenMappedDonor[] snapshots;
        foreach(i,c;donors){fill(*(c.fs),c.pos[0],regime);donorGradients(c);snapshots~=new FrozenMappedDonor(c,donorMotherIDs[i]);}
        frozen=new RAHFrozenMappedCellCopy(block.id,bindex,snapshots,interfaceIDs);frozen.set_up_cell_mapping();
        frozen.exchange_geometry_phase2();need(frozen.mapped_cells.length==0,"RAH_NK_DONOR_REFERENCE");
        makeDirection(regime,directions);
        auto r0=realSample(false,regime,0),f0=realSample(true,regime,0);
        double zedError=0;auto rj=nativeJv(false,regime,zedError),fj=nativeJv(true,regime,zedError);
        auto rp=realSample(false,regime,1e-4),rm=realSample(false,regime,-1e-4);
        auto fp=realSample(true,regime,1e-4),fm=realSample(true,regime,-1e-4);
        double equivalent=0,fdR=0,fdF=0;
        need(f0.length==15918,"RAH_NK_COMPONENTS");
        foreach(i;0..f0.length){
            double ar=rj[i].im/1e-30,af=fj[i].im/1e-30;
            foreach(pair;[[r0[i].re,f0[i].re],[ar,af],[rp[i].re,fp[i].re],[rm[i].re,fm[i].re]])
                equivalent=max(equivalent,fabs(pair[0]-pair[1])/max(1.,max(fabs(pair[0]),fabs(pair[1]))));
            double br=(rp[i].re-rm[i].re)/2e-4,bf=(fp[i].re-fm[i].re)/2e-4;
            fdR=max(fdR,fabs(ar-br)/max(1.,max(fabs(ar),fabs(br))));
            fdF=max(fdF,fabs(af-bf)/max(1.,max(fabs(af),fabs(bf))));
            string kind;size_t id,component;
            if(i<21*block.cells.length){auto j=i%21;id=block.cells[i/21].id;
                kind=j<6?"residual":j<9?"temperatureGradient":j<18?"velocityGradient":"saGradient";
                component=j<6?j:j<9?j-6:j<18?j-9:j-18;
            } else {auto j=i-21*block.cells.length;id=block.faces[j/6].id;component=j%6;kind="integratedFlux";}
            raw.writefln("%s,%s,%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",regime,i,kind,id,component,
                r0[i].re,f0[i].re,ar,af,rp[i].re,rm[i].re,fp[i].re,fm[i].re);
        }
        raw.flush();bool pass=equivalent<=1e-11 && zedError<=1e-11 && fdR<=1e-5 && fdF<=1e-5;
        output.writefln("{\"regime\":%s,\"components\":%s,\"equivalenceError\":%.17g,\"zedError\":%.17g,\"residentFDError\":%.17g,\"frozenFDError\":%.17g,\"passed\":%s}",
            regime,f0.length,equivalent,zedError,fdR,fdF,pass?"true":"false");output.flush();accepted=accepted&&pass;
    }
    directions.close();raw.close();output.close();
    need(accepted,"RAH_NK_DERIVATIVE_GATE");writeln("RAH_NATIVE_NK_OPERATOR_PASS");return 0;
}
