// 模块职责：使用完整原生残差流程对照冻结边界；没有推进/牛顿迭代入口。
module rah_native_file_probe;
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


// 模块四：真实文件边界的独立进程控制，仅调用六组交换，无残差、Jv或时间推进。
import lmr.bc.ghost_cell_effect.rah_frozen_file;
import std.exception : collectException;

FileFaceBinding[] bindings() {
    auto data=parseJSON(readText("donors.json"));JSONValue[string] rows;
    foreach(row;data["donors"].array){size_t[] ids;foreach(x;row["faceVertices"].array)ids~=to!size_t(x.integer);rows[key(ids)]=row;}
    foreach(i,b;block.bc)if(b.label=="exchange")bindex=to!int(i);
    FileFaceBinding[] result;
    foreach(face;block.bc[bindex].faces){size_t[] ids;foreach(v;face.vtx)ids~=v.id;
        auto row=rows[key(ids)];FileFaceBinding b;b.receiverFaceID=face.id;b.motherCellID=to!size_t(row["motherCellID"].integer);
        foreach(x;row["motherFaceNodes"].array)b.motherFaceNodes~=to!size_t(x.integer);
        Vector3[] p;foreach(x;row["vertices"].array)p~=vec(x);Vector3 pos;number vol,len;
        wedge_properties(p[0],p[1],p[2],p[3],p[4],p[5],false,pos,vol,len);
        b.geometry=[pos.x,pos.y,pos.z,vol,to!number(0),len,len,len,len,len];result~=b;
    }
    return result;
}
void exchanges(RAHFrozenMappedCellCopy bc) {
    bc.exchange_geometry_phase2();bc.exchange_flowstate_phase2(0,0,0);
    bc.exchange_convective_gradient_phase2(0,0,0);bc.exchange_viscous_gradient_phase2(0,0,0);
    bc.exchange_turbulent_transprops_phase2();bc.exchange_shock_phase2(0,0,0);
    requireFile(bc.mapped_cells.length==0,"RETAINED_DONOR");
}
void mutate(ref JSONValue value,JSONValue[] path,JSONValue replacement) {
    if(path.length==0){value=replacement;return;}
    if(path[0].type==JSONType.string)mutate(value[path[0].str],path[1..$],replacement);
    else mutate(value[to!size_t(path[0].integer)],path[1..$],replacement);
}
int main(string[] args) {
    requireFile(args.length==2,"MODE");GlobalConfig.in_mpi_context=false;GlobalConfig.is_master_task=true;
    initNewtonKrylovSimulation(0,1,1,"60");requireFile(localFluidBlocks.length==1,"BLOCKS");block=localFluidBlocks[0];
    auto meta=parseJSON(readText("binding.json"));auto map=bindings();
    if(args[1]=="export") {
        setupDonors();foreach(c;donors){fill(*(c.fs),c.pos[0],0);donorGradients(c);}
        exportFile("native-snapshot.json",meta,map,donors);
        auto out=File("export.json","w");out.writefln("{\"passed\":true,\"records\":%s,\"scalarsPerRecord\":140,\"CFDAdvances\":0}",map.length);
        return 0;
    }
    requireFile(args[1]=="import","MODE");auto requests=parseJSON(readText("requests.json"));
    auto output=File("import-results.jsonl","w");
    foreach(item;requests["positive"].array){auto path=item["path"].str;
        auto bc=fileBoundary(block,bindex,path,item["sha256"].str,meta,map);exchanges(bc);
        JSONValue[] rows;foreach(i,g;bc.ghost_cells){auto row=bindingJSON(map[i]);row["values"]=encode(packCell(g));rows~=row;}
        JSONValue result;result["metadata"]=meta;result["records"]=JSONValue(rows);
        auto f=File(item["output"].str,"wx");f.write(result.toString());f.close();
        output.writefln("{\"name\":\"%s\",\"passed\":true,\"residentDonorReferences\":%s}",item["name"].str,bc.mapped_cells.length);output.flush();
    }
    auto original=readText("native-snapshot.json");
    foreach(item;requests["negative"].array){auto doc=parseJSON(original);auto expected=meta;
        if(auto p="path" in item.object)mutate(doc,(*p).array,item["replacement"]);
        if(auto p="expectedPath" in item.object)mutate(expected,(*p).array,item["replacement"]);
        string raw=doc.toString();if(auto p="truncate" in item.object)raw=raw[0..raw.length/2];
        auto file=File("negative-current.json","w");file.write(raw);file.close();
        string expectedSHA=digest(raw);if(auto p="wrongSHA" in item.object)expectedSHA="0000000000000000000000000000000000000000000000000000000000000000";
        bool rejected=false;string reason;
        try{auto ignored=importFile("negative-current.json",expectedSHA,expected,map,block.myConfig);}
        catch(Exception e){rejected=true;reason=e.msg;}
        requireFile(rejected,"NOT_REJECTED_"~item["name"].str);
        JSONValue record;record["name"]=item["name"];record["passed"]=true;record["reason"]=reason;
        output.writeln(record.toString());output.flush();
    }
    // 实际配置拒绝，不依赖文件中自报的配置字符串。
    foreach(mode;0..4){bool oldMHD=block.myConfig.MHD,oldElectric=block.myConfig.electric_field_work;
        int oldDimensions=block.myConfig.dimensions;auto oldInterp=block.myConfig.thermo_interpolator;
        if(mode==0)block.myConfig.MHD=true;if(mode==1)block.myConfig.electric_field_work=true;
        if(mode==2)block.myConfig.dimensions=2;if(mode==3)block.myConfig.thermo_interpolator=InterpolateOption.pT;
        bool rejected=false;try{checkConfig(block.myConfig);}catch(Exception e){rejected=true;}
        block.myConfig.MHD=oldMHD;block.myConfig.electric_field_work=oldElectric;
        block.myConfig.dimensions=oldDimensions;block.myConfig.thermo_interpolator=oldInterp;
        requireFile(rejected,"CONFIG_NOT_REJECTED");output.writefln("{\"name\":\"runtime-config-%s\",\"passed\":true}",mode);
    }
    output.close();writeln("RAH_NATIVE_FILE_BRIDGE_PASS");return 0;
}
