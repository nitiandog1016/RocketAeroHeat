"""模块职责：复用控制输入和准确源码，生成无残差/无推进的原生文件桥接 helper。"""
from pathlib import Path
import hashlib,json,shutil,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v215_validation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')

MAIN=r'''
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
'''

def main():
    start=time.monotonic();out=BASE/'prepared01';out.mkdir(exist_ok=False)
    build=json.loads((ROOT/'outputs/v214_validation/build01/report.json').read_text())
    sources={}
    for rel in ['gas/gas_state.d','lmr/flowstate.d','lmr/lsqinterp.d','lmr/flowgradients.d']:
        actual=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source/src')/rel
        assert sha(actual)==build['inputSHA256'][str(actual)];sources[rel]=sha(actual)
    old=ROOT/'outputs/v213_validation/control02';case=BASE/'control01';case.mkdir(exist_ok=False)
    for name in ['donors.json','mapping.npz','mesh.su2','gas-model.lua']:
        if (old/name).exists():shutil.copyfile(old/name,case/name)
    shutil.copytree(old/'lmrsim',case/'lmrsim')
    d=json.loads((case/'donors.json').read_text());z=np.load(str(case/'mapping.npz'))
    for row in d['donors']:row['motherFaceNodes']=sorted(int(z['motherNodeIDs'][i]) for i in row['faceVertices'])
    save(case/'donors.json',d)
    source=ROOT/'native/eilmerFrozenMapped/probe.d'
    prefix=source.read_text(encoding='utf8').split('number[] sample(')[0].replace('module rah_frozen_mapped_probe;','module rah_native_file_probe;')
    (out/'file_probe.d').write_text(prefix+MAIN,encoding='utf8')
    operator=ROOT/'native/eilmerFrozenMapped/lmr/bc/ghost_cell_effect/rah_frozen_file.d'
    meta=dict(schema='RAH.FrozenMappedSnapshot.v2',model='rah-phase-aware-dry-air-eilmer-v1',
        dataSHA256=sha(ROOT/'data/thermo/phase_aware_dry_air_v1.json'),sourceCommit=build['sourceCommit'],
        motherMeshSHA256='22b35b72396c1a7d14195af9bd24d3c6328138484d4009638c1ab48cddec1d6c',
        geometryMappingSHA256=sha(case/'donors.json'),initialFrameSHA256=sha(case/'lmrsim/snapshots/0000/fluid-0000'),
        sourceRecipeSHA256=sha(out/'file_probe.d'),discreteConfigSHA256=sha(case/'lmrsim/config'),
        adapterSHA256=sha(operator),nativeBuildSHA256='PENDING_BUILD',
        donorRegion='v213-control-mother-neighbors',receiverRegion='v213-control-13-columns',
        origin='manufactured-native-donor-state-not-a-flow-solution',outerIteration=0,donorIteration=0,sweepPosition=0,
        physicalTimeS=0.0,configuration='3D-single-neutral-SA-Edwards-rhou-noMHD-noElectric',
        layout='geometry10-state25-recon14x6-viscous7x3-complexHex64')
    save(out/'binding-template.json',meta)
    save(out/'report.json',dict(passed=True,sourceSHA256=sources,helperSHA256=sha(out/'file_probe.d'),
        inputSHA256={p.relative_to(case).as_posix():sha(p) for p in case.rglob('*') if p.is_file()},seconds=time.monotonic()-start))
    print(json.dumps(dict(passed=True,seconds=time.monotonic()-start)))

if __name__=='__main__':main()
