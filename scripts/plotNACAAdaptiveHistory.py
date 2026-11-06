"""模块用途：绘制实际自适应反馈和异常发展史；原固定历史仅参照，不补造异常瞬间场。"""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from su2ThermalRuntime import read,sha
from validateIdealAir import save


def history(path):
    """模块一：读取真实完整有限历史，不插值或按末场代填历史。"""
    with path.open() as f:
        rows=[{k.strip():float(v) for k,v in row.items()} for row in csv.DictReader(f,skipinitialspace=True)]
    assert rows and all(np.isfinite(list(row.values())).all() for row in rows)
    assert [row['Inner_Iter'] for row in rows]==list(range(len(rows)))
    return {key:np.array([row[key] for row in rows]) for key in rows[0]}


def plot(folder,out):
    """模块二：冷重启更新编号对齐原历史，保留局部CFL统计和有符号功率。"""
    out.mkdir(parents=True,exist_ok=False)
    (out/'executed_plot.py').write_bytes(Path(__file__).read_bytes())
    manifest=read(folder/'case_manifest.json');receipt=read(folder/'engine_process.json')
    assert receipt['status']=='terminal' and all(receipt['inputsStable'].values())
    init=manifest['initialization'];source=Path(init['sourceDirectory'])
    original=read(source/'engine_process.json')
    current_path=folder/'history.csv';source_path=source/'history.csv'
    assert sha(current_path)==receipt['outputSHA256']['history.csv']
    assert sha(source_path)==original['outputSHA256']['history.csv']
    data=history(current_path);base=history(source_path);x=data['Inner_Iter']
    source_steps=init['priorSuccessfulUpdates']+x.astype(int)
    common=source_steps<len(base['Inner_Iter']);bx=x[common];bi=source_steps[common]
    plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':10})
    fig,axes=plt.subplots(3,1,figsize=(10,8.5),sharex=True,layout='constrained')
    ax=axes[0]
    ax.fill_between(x,data['Min CFL'],data['Max CFL'],color='#326b91',alpha=.17,label='原生局部CFL范围')
    ax.plot(x,data['Avg CFL'],color='#326b91',label='原生平均CFL')
    ax.axhline(.4,color='.4',ls='--',label='原固定CFL参照0.4')
    ax.set_yscale('log');ax.set_ylabel('下一步CFL（对数轴）');ax.legend(ncol=3,fontsize=8)
    ax=axes[1]
    for key,color,label in [('rms[Rho]','#326b91','密度'),('rms[RhoE]','#ae593e','能量')]:
        ax.plot(x,data[key],color=color,label='自适应 '+label)
        ax.plot(bx,base[key][bi],color=color,ls='--',alpha=.65,label='固定历史参照 '+label)
    ax.axhline(-8,color='.4',lw=.8,ls=':',label='原−8门槛')
    ax.set_ylabel('log10 RMS残差');ax.legend(ncol=3,fontsize=8)
    ax=axes[2]
    ax.plot(x,data['HF'],color='#ae593e',label='自适应实际历史')
    ax.plot(bx,base['HF'][bi],color='.35',ls='--',label='固定历史参照')
    ax.axhline(0,color='.4',lw=.8);ax.set_ylabel('四分之一域净入壁热功率 / W')
    ax.set_xlabel('本段局部迭代编号');ax.legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2);ax.set_xlim(0,max(1,int(x[-1])))
    saved=[row['localIteration'] for row in receipt['checkpoints']['complete']]
    fig.suptitle('NACA TN3344 · 固定原网格与物性 · 原生自适应推进\n'
        '完整场 %s；最后有限历史 %d；退出码 %d；异常瞬间全场未知\n'
        '源保存步 %d，按更新编号对齐；冷重启非连续等价，参照差异不全归因于自适应'%
        (str(saved),int(x[-1]),receipt['exitCode'],init['sourceLocalIteration']),fontsize=11)
    target=out/'adaptive_history.png';fig.savefig(str(target),dpi=160);plt.close(fig)
    inputs=[current_path,source_path,folder/'engine_process.json',folder/'case_manifest.json',Path(__file__).resolve()]
    save(out/'manifest.json',dict(figure=target.name,inputSHA256={str(p):sha(p) for p in inputs},
        completeSavedLocalIterations=saved,lastFiniteHistoryLocalIteration=int(x[-1]),
        referenceIsUninterruptedEquivalent=False,abnormalInstantStateAvailable=False,physicalQualified=False))
    print(str(target))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',required=True);parser.add_argument('--output',required=True)
    a=parser.parse_args();plot(Path(a.folder).resolve(),Path(a.output).resolve())
