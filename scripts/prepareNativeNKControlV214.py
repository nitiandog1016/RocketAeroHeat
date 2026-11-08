"""模块职责：绑定准确调用源，复用v213输入，生成只调用成熟Jv的隔离控制。"""
from pathlib import Path
import json,shutil,hashlib,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v214_validation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

MAIN=r'''
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
'''

def main():
    start=time.monotonic();prepared=BASE/'prepared01';prepared.mkdir(exist_ok=False)
    b=json.loads((ROOT/'outputs/v213_validation/build02/report.json').read_text())
    src=ROOT/'outputs/v213_validation/source02/source/src/lmr/newtonkrylovsolver.d'
    assert sha(src)==b['inputSHA256']['/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/newtonkrylovsolver.d']
    paths=[src,ROOT/'tools/runtime/eilmer-nk-geometric-v1/overlay/src/lmr/fluidfvcell.d',
           ROOT/'tools/runtime/eilmer-explicit-dense-v1/overlay/src/lmr/fluidblock.d']
    for p in paths[1:]:assert sha(p)==b['inputSHA256'][str(p)]
    source=ROOT/'native/eilmerFrozenMapped/probe.d';assert sha(source)==b['inputSHA256'][str(source)]
    prefix=source.read_text(encoding='utf8').split('number[] sample(')[0]
    prefix=prefix.replace('module rah_frozen_mapped_probe;','module rah_native_nk_probe;')
    prefix=prefix.replace('FluidFVCell[] donors;','FluidFVCell[] donors;\nsize_t[] donorMotherIDs;')
    prefix=prefix.replace('donors~=c;interfaceIDs~=face.id;','donors~=c;interfaceIDs~=face.id;donorMotherIDs~=to!size_t(row["motherCellID"].integer);')
    (prepared/'nk_probe.d').write_text(prefix+MAIN,encoding='utf8')
    old=ROOT/'outputs/v213_validation/control02';case=BASE/'control01';case.mkdir(exist_ok=False)
    for p in old.iterdir():
        if p.name=='operator_results.jsonl':continue
        if p.is_dir():shutil.copytree(p,case/p.name)
        else:shutil.copyfile(p,case/p.name)
    report=dict(passed=True,sourceSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        baseProbeSHA256=sha(source),newProbeSHA256=sha(prepared/'nk_probe.d'),
        caseInputSHA256={p.relative_to(case).as_posix():sha(p) for p in case.rglob('*') if p.is_file()},
        sourceContract='native dispatcher; conserved U[1]; order 1 vs residual 2; no Krylov/time-diagonal qualification',
        seconds=time.monotonic()-start,CFDAdvanceRuns=0)
    (prepared/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(passed=True,seconds=report['seconds'])))

if __name__=='__main__':raise ValueError('v214准备已完成；复用准确输入，禁止重建相同控制。')
