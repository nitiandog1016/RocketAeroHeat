"""独立读取实际原生网格文本并画原箭体诊断；无流场或热图。"""
import json,time,re,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from prepareOpenFoamControl import ROOT,save,sha

def main():
    start=time.monotonic();base=ROOT/'outputs/v208_validation';out=base/'artifact_audit02';out.mkdir(exist_ok=False)
    case=base/'layer_mesh02';d=np.load(case/'exact_mesh.npz');mesh=case/'constant/polyMesh'
    # 模块一：独立解析实际写出的文本，而非仅验证写入前的数组。
    def body(name):
        text=(mesh/name).read_text();match=re.search(r'\n(\d+)\n\(\n(.*)\n\)\n// \*',text,re.S)
        assert match;return int(match.group(1)),match.group(2).splitlines()
    n,rows=body('points');p=np.array([[float(x) for x in row.strip('()').split()] for row in rows]);assert n==len(p) and np.array_equal(p,d['points'])
    for name,key in [('owner','owner'),('neighbour','neighbour')]:
        n,rows=body(name);a=np.array([int(x) for x in rows]);assert n==len(a) and np.array_equal(a,d[key])
    n,rows=body('faces');f=np.full((n,4),-1,np.int64);sizes=np.zeros(n,int)
    assert n==len(rows)
    for i,row in enumerate(rows):
        count,indices=row.split('(',1);a=[int(x) for x in indices.rstrip(')').split()];sizes[i]=int(count);assert len(a)==sizes[i];f[i,:len(a)]=a
    assert np.array_equal(f,d['faces']) and np.array_equal(sizes,d['faceSizes'])
    report=json.loads((base/'layer_audit02/report.json').read_text());assert report['passed'] and not report['geometryQualityAccepted']
    # 模块二：沿原面位置显示真实质量值，灰色部分仍保留整箭几何。
    s=json.loads((ROOT/'outputs/v194_validation/rocket_design03/surface.json').read_text());v=np.array(s['vertices']);wall=np.array(s['facesOneBased'])-1
    data=np.loadtxt(base/'layer_audit02/wall_adjacent_determinant.csv',delimiter=',',skiprows=1)
    ids=data[:,0].astype(int)-1;selected=np.zeros(len(wall),bool);selected[ids]=True
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    fig=plt.figure(figsize=(13,5.6));a=fig.add_subplot(121,projection='3d');b=fig.add_subplot(122)
    cmap=plt.get_cmap('viridis');norm=Normalize(-12,-2)
    a.add_collection3d(Poly3DCollection(v[wall[~selected]],facecolor='#c8cdd1',linewidth=0))
    a.add_collection3d(Poly3DCollection(v[wall[ids]],facecolor=cmap(norm(np.log10(data[:,1]))),linewidth=0))
    a.set(xlim=(v[:,0].min(),v[:,0].max()),ylim=(-.4,.4),zlim=(-.4,.4),xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)')
    a.set_box_aspect((np.ptp(v[:,0]),.8,.8));a.view_init(elev=25,azim=-60)
    a.set_yticks([-.4,0,.4]);a.set_zticks([-.4,0,.4])
    a.set_title('Original rocket: first-layer cell diagnostic',fontsize=10)
    sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap);sm.set_array([])
    color_axis=fig.add_axes([.09,.125,.30,.025])
    cb=fig.colorbar(sm,cax=color_axis,orientation='horizontal');cb.set_label('log10(area-tensor determinant)',fontsize=9)
    layers=report['layers'];x=np.arange(1,33)
    b.bar(x,[r['determinantFailed']/7838*100 for r in layers],color='#b9905b',label='All failing cells')
    b.plot(x,[r['fullyInternalFailed']/7838*100 for r in layers],color='#2b657e',marker='.',label='Fully internal failing cells')
    b.set(xlabel='Original normal layer',ylabel='Cells below unchanged 0.001 gate (%)',ylim=(0,105));b.legend(fontsize=8);b.grid(axis='y',alpha=.2)
    b.set_title('139,027 failed cells; 130,223 fully internal',fontsize=10)
    fig.suptitle('v208 | Actual near-wall mesh compatibility — no pressure / heat solution',fontsize=12)
    fig.text(.5,.025,'32 layers and all 7,838 original wall faces retained. Native mesh acceptance: failed.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.02,right=.98,top=.84,bottom=.25,wspace=.28);fig.savefig(out/'regional_layer_diagnostic.png',dpi=150);plt.close(fig)
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert all(sha(ROOT/n.replace('\\','/'))==h for n,h in protected.items())
    shutil.copyfile(__file__,out/'executed_entry.py')
    save(out/'report.json',dict(passed=True,nativeTextEqualsExactArrays=True,cells=len(d['cells']),faces=n,
        correctedMeshNativeCheckPerformed=False,geometryQualityAccepted=False,protectedFilesUnchanged=29,
        field='geometric area-tensor determinant, not pressure/heat',elapsedSeconds=time.monotonic()-start,
        sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [case/'exact_mesh.npz',base/'layer_audit02/report.json',base/'layer_audit02/wall_adjacent_determinant.csv']},
        figureSHA256=sha(out/'regional_layer_diagnostic.png'),CFDRuns=0,physicalQualification=False))
    print(json.dumps(dict(passed=True,seconds=time.monotonic()-start)))

if __name__=='__main__':main()
