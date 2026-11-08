"""模块职责：复用已验证输入，生成只导出指定位置实际点云的隔离helper，不改求解内核。"""
from pathlib import Path
import json,hashlib,shutil,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v217_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

OBSERVER=r'''
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
'''

def main():
    start=time.monotonic(); out=BASE/'prepared01';out.mkdir(exist_ok=False)
    old=ROOT/'outputs/v214_validation';source=old/'prepared01/nk_probe.d'
    receipt=json.loads((old/'build01/report.json').read_text());expected=receipt['inputSHA256']['/mnt/e/RocketAeroHeat/outputs/v214_validation/prepared01/nk_probe.d']
    assert sha(source)==expected
    prefix=source.read_text(encoding='utf8').split('// 模块四：')[0].replace('module rah_native_nk_probe;','module rah_regional_stencil_probe;')
    target=out/'stencil_probe.d';target.write_bytes((prefix+OBSERVER).encode('utf8'))
    case=BASE/'control01';case.mkdir(exist_ok=False)
    for name in ['phase-air.lua','directions.csv','donors.json','job.lua','grid.lua','mesh.su2']:
        shutil.copyfile(old/'control01'/name,case/name)
    shutil.copytree(old/'control01/lmrsim',case/'lmrsim')
    report=dict(passed=True,oldProbeSHA256=expected,probeSHA256=sha(target),
                inputSHA256={p.relative_to(case).as_posix():sha(p) for p in case.rglob('*') if p.is_file()},
                seconds=time.monotonic()-start,CFDAdvanceRuns=0)
    (out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'))
    print(json.dumps(dict(passed=True,seconds=report['seconds'])))

if __name__=='__main__':main()
