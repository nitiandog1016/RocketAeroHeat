// 模块职责：把完整原生 mapped 复制量保存为具名文件，再导入既有冻结边界。
// 只做接口、身份和状态审核；不计算通量/Jacobian，也不授予 CFD 运行资格。
module lmr.bc.ghost_cell_effect.rah_frozen_file;

import std.stdio : File;
import std.file : readText, exists, getSize;
import std.json;
import std.conv : to;
import std.format : format;
import std.math : isFinite, isNaN, fabs;
import std.algorithm : max, sort;
import std.string : toLower;
import std.digest.sha : sha256Of;
import std.digest : toHexString;
import geom;
import gas;
import gas.rah_phase_aware_air;
import nm.number;
import ntypes.complex;
import lmr.globalconfig;
import lmr.coredata;
import lmr.fluidfvcell;
import lmr.fluidblock;
import lmr.flowstate;
import lmr.lsqinterp;
import lmr.bc.ghost_cell_effect.rah_frozen_mapped;

enum schemaV2="RAH.FrozenMappedSnapshot.v2";
enum scalarCount=140;
void requireFile(bool test,string code) {if(!test)throw new Exception("RAH_FILE_"~code);}
string digest(string data) {return toLower(toHexString(sha256Of(data)));}

// 模块一：逐个浮点分量按位编码。保留未用 NaN、负零，不序列化指针或结构体 padding。
union Bits {double d; ulong u;}
string hexDouble(double x) {Bits b;b.d=x;return format("%016x",b.u);}
double unhex(string s) {
    requireFile(s.length==16,"HEX_LENGTH");ulong n=0;
    foreach(c;s){requireFile((c>='0'&&c<='9')||(c>='a'&&c<='f'),"HEX_CHAR");n=(n<<4)|(c<='9'?c-'0':c-'a'+10);}
    Bits b;b.u=n;return b.d;
}
JSONValue encode(number[] values) {
    JSONValue[] a;foreach(q;values)a~=JSONValue(hexDouble(q.re)~hexDouble(q.im));return JSONValue(a);
}
number[] decode(JSONValue a) {
    requireFile(a.type==JSONType.array && a.array.length==scalarCount,"SCALAR_COUNT");
    number[] values;foreach(x;a.array){requireFile(x.type==JSONType.string && x.str.length==32,"SCALAR_ENCODING");
        values~=complex(unhex(x.str[0..16]),unhex(x.str[16..32]));}return values;
}

// 字段顺序独立于内存布局；与 v2 文档中的 geometry/state/reconstruction/viscous 相同。
number[] packCell(FluidFVCell c) {
    number[] a=[c.pos[0].x,c.pos[0].y,c.pos[0].z,c.volume[0],c.areaxy[0],
        c.iLength,c.jLength,c.kLength,c.L_min,c.L_max];
    auto f=c.fs;auto g=&f.gas;
    a~=[g.rho,g.p,g.p_e,g.a,g.T,g.u,g.mu,g.k,g.sigma,g.quality,g.massf[0],g.rho_s[0],
        f.vel.x,f.vel.y,f.vel.z,f.turb[0],f.turb[1],f.mu_t,f.k_t,f.S,f.B.x,f.B.y,f.B.z,f.psi,f.divB];
    auto r=c.gradients;
    static foreach(name;["velx","vely","velz","rho","p","T","u","Bx","By","Bz","psi"])
        mixin("a~=[r."~name~".x,r."~name~".y,r."~name~".z,r."~name~"Phi,r."~name~"Min,r."~name~"Max];");
    foreach(j;0..2)a~=[r.turb[j].x,r.turb[j].y,r.turb[j].z,r.turbPhi[j],r.turbMin[j],r.turbMax[j]];
    a~=[r.rho_s[0].x,r.rho_s[0].y,r.rho_s[0].z,r.rho_sPhi[0],r.rho_sMin[0],r.rho_sMax[0]];
    foreach(row;c.grad.vel)a~=row[];
    a~=c.grad.T[];a~=c.grad.massf[0][];foreach(row;c.grad.turb)a~=row[];
    requireFile(a.length==scalarCount,"PACK_LAYOUT");return a;
}
void unpackCell(FluidFVCell c,number[] a) {
    requireFile(a.length==scalarCount,"SCALAR_COUNT");size_t j=0;
    c.pos[0]=Vector3(a[0],a[1],a[2]);c.volume[0]=a[3];c.areaxy[0]=a[4];
    c.iLength=a[5];c.jLength=a[6];c.kLength=a[7];c.L_min=a[8];c.L_max=a[9];j=10;
    auto f=c.fs;auto g=&f.gas;
    static foreach(name;["rho","p","p_e","a","T","u","mu","k","sigma","quality"])
        mixin("g."~name~"=a[j++];");
    g.massf[0]=a[j++];g.rho_s[0]=a[j++];f.vel=Vector3(a[j],a[j+1],a[j+2]);j+=3;
    f.turb[0]=a[j++];f.turb[1]=a[j++];f.mu_t=a[j++];f.k_t=a[j++];f.S=a[j++];
    f.B=Vector3(a[j],a[j+1],a[j+2]);j+=3;f.psi=a[j++];f.divB=a[j++];
    auto r=c.gradients;
    static foreach(name;["velx","vely","velz","rho","p","T","u","Bx","By","Bz","psi"])
        mixin("r."~name~"=Vector3(a[j],a[j+1],a[j+2]);j+=3;r."~name~"Phi=a[j++];r."~name~"Min=a[j++];r."~name~"Max=a[j++];");
    foreach(i;0..2){r.turb[i]=Vector3(a[j],a[j+1],a[j+2]);j+=3;r.turbPhi[i]=a[j++];r.turbMin[i]=a[j++];r.turbMax[i]=a[j++];}
    r.rho_s[0]=Vector3(a[j],a[j+1],a[j+2]);j+=3;r.rho_sPhi[0]=a[j++];r.rho_sMin[0]=a[j++];r.rho_sMax[0]=a[j++];
    foreach(ref row;c.grad.vel)foreach(ref q;row)q=a[j++];
    foreach(ref q;c.grad.T)q=a[j++];foreach(ref q;c.grad.massf[0])q=a[j++];
    foreach(ref row;c.grad.turb)foreach(ref q;row)q=a[j++];
    requireFile(j==scalarCount,"UNPACK_LAYOUT");
}

// 模块二：只接受已经核对的物理/编译组合，其他配置不能借用未用槽规则。
void checkConfig(LocalConfig c) {
    requireFile(!GlobalConfig.in_mpi_context && c.dimensions==3 && !c.axisymmetric && !c.MHD &&
        !c.electric_field_work && c.n_grid_time_levels==1 && c.n_species==1 && c.n_modes==0 &&
        c.turb_model.nturb==1 && c.turb_model.modelName=="spalart_allmaras_edwards" && c.thermo_interpolator==InterpolateOption.rhou &&
        !c.gmodel.is_plasma && (cast(RAHPhaseAwareAir)c.gmodel)!is null,"CONFIG");
}
void active(number q) {requireFile(isFinite(q.re)&&q.im==0,"ACTIVE_NONREAL_OR_NONFINITE");}
void checkCell(FluidFVCell c) {
    checkConfig(c.myConfig);auto a=packCell(c);
    foreach(q;a)requireFile((isFinite(q.re)||isNaN(q.re)) && (q.im==0||isNaN(q.im)),"UNUSED_INFINITY_OR_COMPLEX");
    foreach(i;[0,1,2,3,5,6,7,8,9])active(a[i]);
    foreach(i;[3,5,6,7,8,9])requireFile(a[i].re>0,"GEOMETRY");
    requireFile(a[8].re<=a[9].re,"GEOMETRY");
    foreach(i;[0,1,3,4,5,6,7,8,9,10,12,13,14,15,17,18,19])active(a[10+i]);
    requireFile(a[20].re==1 && a[18].re==0 && a[19].re==1,"NEUTRAL_AIR");
    foreach(i;[15,17,18])requireFile(a[10+i].re>=0,"SA_TRANSPORT");
    // 仅在临时气体副本上核验 EOS/输运；绝不覆写文件的原生值或补 rho_s/h0。
    auto gasCopy=GasState(c.fs.gas);auto gm=c.myConfig.gmodel;
    gm.update_thermo_from_pT(gasCopy);gm.update_sound_speed(gasCopy);gm.update_trans_coeffs(gasCopy);
    auto wanted=[gasCopy.rho,gasCopy.a,gasCopy.u,gasCopy.mu,gasCopy.k];
    foreach(j,i;[0,3,5,6,7]){double floor=(i==3||i==5)?1.:1e-12;
        requireFile(fabs(a[10+i].re-wanted[j].re)<=5e-10*max(floor,fabs(wanted[j].re)),"THERMO");}
    foreach(i;0..14){auto r=a[35+6*i..41+6*i];bool used=i<4||i==6||i==11;
        // p/T/rho_s 若明确有 limiter，也检查其可用性；关闭 MHD/第二 SA 槽不参与计算。
        if((i==4||i==5||i==13)&&r[3].re!=-1)used=true;
        if(used){foreach(q;r)active(q);requireFile(r[3].re>=0&&r[4].re<=r[5].re,"LIMITER");}}
    foreach(q;a[119..137])active(q);
    foreach(q;a[131..134])requireFile(q.re==0,"SINGLE_SPECIES_GRADIENT");
    requireFile(c.fs.gas.ceaSavedData is null,"CEA_POINTER");
}

// 模块三：文件必须匹配调用方独立绑定的源帧/轮次与母面映射，不信任文件自称身份。
struct FileFaceBinding {
    size_t receiverFaceID, motherCellID;
    size_t[] motherFaceNodes;
    number[] geometry;
}
JSONValue bindingJSON(FileFaceBinding b) {
    JSONValue v;v["receiverFaceID"]=to!long(b.receiverFaceID);v["motherCellID"]=to!long(b.motherCellID);
    long[] ids;foreach(i;b.motherFaceNodes)ids~=to!long(i);v["motherFaceNodes"]=JSONValue(ids);return v;
}
void checkMetadata(JSONValue meta,JSONValue expected) {
    requireFile(meta==expected,"IDENTITY");requireFile(meta["schema"].str==schemaV2,"SCHEMA");
    requireFile(meta["model"].str==RAHPhaseAwareAir.identity && meta["dataSHA256"].str==RAHPhaseAwareAir.dataSHA,"MODEL");
    requireFile(meta["layout"].str=="geometry10-state25-recon14x6-viscous7x3-complexHex64" &&
        meta["configuration"].str=="3D-single-neutral-SA-Edwards-rhou-noMHD-noElectric","LAYOUT");
    foreach(k;["outerIteration","donorIteration","sweepPosition"])requireFile(meta[k].integer>=0,"ITERATION");
    auto t=meta["physicalTimeS"].floating;requireFile(isFinite(t)&&t>=0,"TIME");
}
void exportFile(string path,JSONValue meta,FileFaceBinding[] binding,FluidFVCell[] donors) {
    requireFile(!exists(path),"EXISTS");checkMetadata(meta,meta);requireFile(binding.length==donors.length&&donors.length>0,"COUNT");
    JSONValue[] records;foreach(i,c;donors){checkCell(c);auto a=packCell(c);
        requireFile(encode(a[0..10])==encode(binding[i].geometry),"SOURCE_GEOMETRY");
        auto row=bindingJSON(binding[i]);row["values"]=encode(a);records~=row;}
    JSONValue doc;doc["metadata"]=meta;doc["records"]=JSONValue(records);
    auto f=File(path,"wx");f.write(doc.toString());f.close();
}
FrozenMappedDonor[] importFile(string path,string expectedSHA,JSONValue expected,FileFaceBinding[] binding,LocalConfig config) {
    checkConfig(config);requireFile(getSize(path)<=16*1024*1024,"SIZE");auto raw=readText(path);
    requireFile(digest(raw)==expectedSHA,"SHA");auto doc=parseJSON(raw);
    requireFile(doc.object.length==2 && ("metadata" in doc.object)!is null && ("records" in doc.object)!is null,"TOP_FIELDS");
    checkMetadata(doc["metadata"],expected);auto rows=doc["records"].array;
    requireFile(rows.length==binding.length&&binding.length>0,"COUNT");
    size_t[size_t] byID;foreach(j,row;rows){requireFile(row.object.length==4,"ROW_FIELDS");
        auto id=to!size_t(row["receiverFaceID"].integer);requireFile((id in byID)is null,"DUPLICATE_FACE");byID[id]=j;}
    // 单个 staging cell 只用于既有深拷贝构造；返回对象没有供体 cell/block 的引用。
    FluidCellData data;data.allocate_minimal(1,0,config);data.lsqgradients~=LSQInterpGradients(1,0,1);
    auto staging=new FluidFVCell(config,&data,0);FrozenMappedDonor[] result;
    foreach(b;binding){auto index=b.receiverFaceID in byID;requireFile(index!is null,"FACE_ID");auto row=rows[*index];
        auto saved=row["values"];row.object.remove("values");requireFile(row==bindingJSON(b),"MOTHER_MAPPING");
        auto a=decode(saved);requireFile(encode(a[0..10])==encode(b.geometry),"SOURCE_GEOMETRY");
        unpackCell(staging,a);checkCell(staging);result~=new FrozenMappedDonor(staging,b.motherCellID);}
    return result;
}
RAHFrozenMappedCellCopy fileBoundary(FluidBlock block,int boundary,string path,string expectedSHA,
                                  JSONValue expected,FileFaceBinding[] binding) {
    auto saved=importFile(path,expectedSHA,expected,binding,block.myConfig);
    size_t[] ids;foreach(b;binding)ids~=b.receiverFaceID;
    auto bc=new RAHFrozenMappedCellCopy(block.id,boundary,saved,ids);bc.set_up_cell_mapping();return bc;
}
