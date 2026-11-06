"""模块用途：导出原面空间壁温研究场的科研图；只读保存帧，不产生物理资格。"""
import argparse
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from su2NACALaminarInterface import mesh,load_state
from su2ThermalRuntime import read,sha
from validateIdealAir import save


def plot(folder,audit_folder,out,iteration=None):
    """模块一：原节点字段、给定/流体壁温、热符号和保存步身份随图输出。"""
    out.mkdir(parents=True,exist_ok=False)
    receipt=read(folder/'engine_process.json');assert receipt['status']=='terminal'
    pairs=receipt['checkpoints']['complete'];assert pairs
    pair=next((p for p in pairs if p['localIteration']==iteration),None) if iteration is not None else pairs[-1]
    assert pair is not None
    step=pair['localIteration'];path=folder/pair['volume'];assert sha(path)==pair['volumeSHA256']
    points,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json');a,v=load_state(path)
    # 模块二：图像只消费同一保存场的已通过热审核，不把其他目录的表拼入当前云图。
    audit=read(audit_folder/'report.json')
    assert audit.get('savedFrameInterfacePassed',audit.get('interfacePassed',False))
    assert audit['inputSHA256'][str(path)]==sha(path)
    checked=next(f for f in audit['frames'] if f['localIteration']==step)
    assert checked['interfacePassed']
    wall_csv=audit_folder/checked['wallCSV'];assert sha(wall_csv)==checked['wallCSV_SHA256']
    rows=np.genfromtxt(str(wall_csv),delimiter=',',names=True,dtype=None,encoding='utf8')
    ids=rows['PointID'].astype(int)
    assert len(np.unique(ids))==len(ids)
    assert np.max(abs(np.column_stack([rows[k] for k in ('X_m','Y_m','Z_m')])-points[ids]))<2e-13
    fields=dict(pressure=v['Pressure'],heat=np.full(len(points),np.nan))
    fields['heat'][ids]=rows['qIntoWall_Wm2']/1000
    faces=[]
    for row in markers['wall']:
        f=row[1:]
        for j in range(1,len(f)-1):faces.append([f[0],f[j],f[j+1]])
    faces=np.array(faces,int);xyz=points[faces]*1000
    plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':10})
    fig=plt.figure(figsize=(13,8),layout='constrained')
    titles=[('pressure','壁面压力 / kPa','viridis'),('heat',r'入壁热流 / kW m$^{-2}$','coolwarm')]
    for panel,(key,title,cmap) in enumerate(titles,1):
        ax=fig.add_subplot(2,2,panel,projection='3d')
        value=np.mean(fields[key][faces],axis=1)/(1000 if key=='pressure' else 1)
        if key=='heat':
            limit=max(abs(value.min()),abs(value.max()),1e-9);norm=matplotlib.colors.Normalize(-limit,limit)
        else:norm=matplotlib.colors.Normalize(value.min(),value.max())
        colors=plt.get_cmap(cmap)(norm(value))
        collection=Poly3DCollection(xyz,facecolors=colors,edgecolors=colors,linewidths=.05,antialiased=False,rasterized=True)
        ax.add_collection3d(collection);ax.auto_scale_xyz(points[:,0]*1000,points[ids,1]*1000,points[ids,2]*1000)
        ax.set_xlim(-127,0);ax.set_ylim(0,50.8);ax.set_zlim(0,50.8);ax.set_box_aspect((127,50.8,50.8))
        ax.set_xlabel('X / mm');ax.set_ylabel('Y / mm');ax.set_zlabel('Z / mm');ax.view_init(20,-70)
        ax.set_title(title+'（原四分之一壁面）')
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,orientation='horizontal',shrink=.7,fraction=.035,pad=.03)
    with np.load(str(folder/'mesh_arrays.npz')) as arrays:labels=arrays['nodeIndices']
    station=labels[ids,1];s=np.array(meta['s']);mean_q=[];lo_q=[];hi_q=[];mean_t=[];given=[]
    for i in range(len(s)):
        select=station==i;weight=rows['nativeProjectedArea_m2'][select]
        q=rows['qIntoWall_Wm2'][select]/1000
        mean_q.append(float(np.average(q,weights=weight)));lo_q.append(float(q.min()));hi_q.append(float(q.max()))
        mean_t.append(float(np.average(rows['wallFluidT_K'][select],weights=weight)))
        given.append(float(np.average(rows['givenTw_K'][select],weights=weight)))
    ax=fig.add_subplot(2,2,3);ax.plot(s/.1016,mean_q,color='#255f85',label='原节点投影面积加权均值')
    ax.fill_between(s/.1016,lo_q,hi_q,color='#255f85',alpha=.15,label='当前扇区周向范围')
    ax.axhline(0,color='.4',lw=.7);ax.set_xlabel('壁面弧长 s / D');ax.set_ylabel(r'入壁热流 / kW m$^{-2}$')
    ax.set_title('热流正号入壁；本加热实验预期为负');ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=fig.add_subplot(2,2,4);ax.plot(s/.1016,given,color='#222222',label='给定 Tw')
    ax.plot(s/.1016,mean_t,color='#b44c32',label='方程中的壁面流体 T（扇区均值）')
    ax.set_xlabel('壁面弧长 s / D');ax.set_ylabel('温度 / K');ax.grid(alpha=.2);ax.legend(fontsize=8)
    ax.set_title('给定壁温、流体温度与固体求解温度分开')
    normal=receipt['checkpoints']['normalTerminal'];last=receipt['checkpoints']['lastFiniteHistoryLocalIteration']
    manifest=read(folder/'case_manifest.json');origin=manifest['initialization'];prior=origin['priorSuccessfulUpdates']
    restart_label=('\n冷重启自局部步 %s；累计来源更新 %d（不代表不中断等价）'%
        (origin['sourceLocalIteration'],pair['seriesSuccessfulUpdates'])) if prior else ''
    fig.suptitle(r'NACA TN3344 · M1.97 · $Re_D=2.28\times10^6$ · 层流 / 变比热空气'+'\n'
        '已保存局部步 %d；最后有限历史 %s；正常终态 %s；研究场，尚无实验资格'%(step,last,normal)+restart_label,fontsize=12)
    target=out/'naca_spatial_fields.png';fig.savefig(str(target),dpi=170);plt.close(fig)
    (out/'executed_plot.py').write_bytes(Path(__file__).read_bytes())
    save(out/'manifest.json',dict(figure=target.name,localIteration=step,sourceFieldSHA256=sha(path),
        wallCSV_SHA256=sha(wall_csv),thermalAuditSHA256=sha(audit_folder/'report.json'),sourceMeshSHA256=sha(folder/'mesh.su2'),
        surfaceDisplay='original wall triangles; arithmetic nodal face means for display only',
        curveDisplay='original nodal projected-area weighted mean in each s row; not finite electrode reduction',
        iterationOrigin={k:origin.get(k) for k in ('kind','priorSuccessfulUpdates','sourceLocalIteration','uninterruptedTrajectoryEquivalent')},
        seriesSuccessfulUpdates=pair['seriesSuccessfulUpdates'],advancement=manifest.get('advancement'),
        originalQuarterDomain=True,physicalQualified=False))
    print(str(target))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True)
    parser.add_argument('--audit',required=True);parser.add_argument('--output',required=True);parser.add_argument('--iteration',type=int)
    a=parser.parse_args();plot(Path(a.folder).resolve(),Path(a.audit).resolve(),Path(a.output).resolve(),a.iteration)
