"""模块职责：真实完整循环更新原面研究场接既有WebGL，统一色标、方法和未知温度身份。"""
from pathlib import Path
import csv,json,re,subprocess,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm,SymLogNorm,Normalize
import prepareOuterStencilV249 as run
from terminalCapacityV249 import supervise
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;AUDIT=BASE/'audit03';CLOSED=True
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def export():
    out=BASE/'workbench01';out.mkdir(exist_ok=False);audit=read(AUDIT/'report.json');assert audit['passed']
    g={k:v for k,v in np.load(str(GEOMETRY)).items()};frames={s:np.loadtxt(str(AUDIT/('original_faces_%06d.csv'%s)),delimiter=',',skiprows=1) for s in (0,1)}
    stages=sorted(frames);last=stages[-1]
    model=read(ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json');old={f['key']:f for f in model['fields']};fields=[];ranges={}
    method='SU2 8.6 / 原生Newton-Krylov，FGMRES实际4次乘积（源/阶段推定） / 相态Air线性输运 / Roe隐式Euler / SA-neg WITHFT2'
    scope='原19752面，M5、30 km、α20°、给定Tw300 K；从v248200六量建立的新方法一次控制，未收敛、非热资格。'
    specs=[('pressure','壁面压力','Pa',2,'当前壁节点primitive到原三角面顶点均值'),
      ('heatFlux','入壁对流热流','W/m²',3,'分子k×(当前法向邻点温度−给定Tw)/原法向距离；原面顶点均值，正值入壁'),
      ('wallFluid','壁节点流体温度','K',4,'当前fine primitive；不是固体壁温'),
      ('wallDeviation','壁流体温度减给定Tw','K',None,'壁节点流体温度减300 K；不是已解固体温度或热平衡误差')]
    for key,label,unit,col,detail in specs:
        vals={s:a[:,col] if col is not None else a[:,4]-a[:,5] for s,a in frames.items()};lo=min(float(min(v)) for v in vals.values());hi=max(float(max(v)) for v in vals.values())
        if key=='wallDeviation':lo,hi=-max(abs(lo),abs(hi)),max(abs(lo),abs(hi))
        for s,v in vals.items():
            name=key+'_'+str(s);ranges[name]=[lo,hi]
            fields.append(dict(key=name,label=label+' · NK控制'+str(s),unit=unit,values=v.tolist(),status=[1]*19752,
                methodIndex=[0]*19752,methodNames=[method+'；'+detail],scope=scope,computedStateLabel='未收敛研究数值',
                stateNote='真实NK控制更新%d；不是物理时间。初始化值未经过壁条件收敛；研究覆盖不等于热资格。'%s,
                successfulUpdates=s,comparisonVariable=key,availability='available'))
    for key in ('wallInput','wallSolved'):old[key]['scope']=scope;fields.append(old[key])
    model.update(fields=fields,comparisonColorRanges=ranges,defaultField='heatFlux_'+str(last),
      resultCondition=dict(caption='M5 · α20° · 给定Tw300 K · NK新方法控制；0为v248200 · '+str(stages)+' · 共同色标'),
      caseNote='v249原真实箭体单次原生Newton-Krylov控制。压力/热流均为未收敛研究；给定Tw、壁流体T和未知固体Tw分开。尚无新增合格热面积。',
      researchIdentity=dict(initialization='v248200-new-native-NK-control',physicalQualification=False,productionQualified=False,successfulUpdates=stages,
        sourceFaceCSV_SHA256={str(s):sha(AUDIT/('original_faces_%06d.csv'%s)) for s in frames},auditSHA256=sha(AUDIT/'report.json'),
        workingArrayIsFinalResidual=False,terminalResidual=None,nativeControlAudit=dict(actualFactory=audit['actualFactory'],flowSpaceEvaluations=audit['flowSpaceEvaluations'],initialElevenBitEqual=audit['initialElevenBitEqual'],finalNativeSixCoordinateBitEqual=audit['finalNativeSixCoordinateBitEqual'],originalCapacityPassed=False)))
    payload=json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False);target=out/'rocket_fas_control.json';target.write_text(payload+'\n',encoding='utf8')
    (out/'case-data.js').write_text('window.ROCKET_CASE='+payload+';\n',encoding='utf8')
    html=(ROOT/'web/index.html').read_text(encoding='utf8');html=re.sub(r'((?:src|href)=")((?:viewer|section|export)\.(?:js|css)(?:\?[^\"]*)?)(")',r'\1../../../web/\2\3',html)
    (out/'index.html').write_text(html,encoding='utf8');restored=read(target)
    assert restored['faces']==g['originalFaces'].tolist() and restored['vertices']==g['originalVertices'].tolist()
    for field in restored['fields'][:-2]:
        a=frames[field['successfulUpdates']];col={'pressure':2,'heatFlux':3,'wallFluid':4}.get(field['comparisonVariable']);assert np.array_equal(field['values'],a[:,col] if col is not None else a[:,4]-a[:,5])
    # 模块二：科研图保留真实长宽比例，图例在独立色标区，最近两阶段固定同一范围；工作台保留全阶段共同色标。
    tri=g['originalVertices'][g['originalFaces']];right=np.array([1.,0.,.1]);right/=np.linalg.norm(right);up=np.array([-.025,1.,.25]);up-=right*(up@right);up/=np.linalg.norm(up)
    proj=tri@np.column_stack((right,up));order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(3,2,figsize=(13,7.5));fig.subplots_adjust(left=.05,right=.86,top=.87,bottom=.11,hspace=.32,wspace=.18)
    compare=stages[-2:] if len(stages)>1 else [stages[0],stages[0]]
    for row,(key,col,scale,title,cmap) in enumerate((('pressure',2,1e-3,'Pressure [kPa], log','viridis'),('heatFlux',3,1e-6,'Heat into wall [MW/m2], symlog','coolwarm'),('wallDeviation',None,1,'Wall fluid T - given Tw [K]','coolwarm'))):
        plotted=np.concatenate([(frames[t][:,col] if col is not None else frames[t][:,4]-frames[t][:,5])*scale for t in compare]);lo,hi=float(min(plotted)),float(max(plotted))
        if hi==lo:hi=lo+1
        norm=LogNorm(lo,hi) if key=='pressure' else SymLogNorm(.01,vmin=-max(abs(lo),abs(hi)),vmax=max(abs(lo),abs(hi))) if key=='heatFlux' else Normalize(-max(abs(lo),abs(hi)),max(abs(lo),abs(hi)))
        for colplot,step in enumerate(compare):
            ax=axs[row,colplot];a=frames[step];v=(a[:,col] if col is not None else a[:,4]-a[:,5])*scale
            ax.add_collection(PolyCollection(proj[order],linewidths=0,facecolors=plt.get_cmap(cmap)(norm(v[order]))));ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.set_axis_off()
            ax.set_title(title+' | update '+str(step),loc='left',fontsize=10)
        box=axs[row,1].get_position();cax=fig.add_axes([.9,box.y0+.1*box.height,.012,.8*box.height]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
    fig.suptitle('RocketAeroHeat v249 | Original 19,752 faces | M5 / alpha 20 deg / given Tw 300 K\nNative Newton-Krylov: ONE controlled update from v248200, NOT a qualified thermal solution',fontsize=12)
    fig.text(.05,.035,'Phase-aware Air / Native JFNK + Roe / SA-neg WITHFT2 | Same scale for each pair | Signed heat flux retained\nGiven wall temperature, fluid wall-node temperature and unsolved solid temperature are separate.',fontsize=9)
    fig.savefig(str(out/'native_newton_control.png'),dpi=160);plt.close(fig)
    result=subprocess.run(['node',str(ROOT/'tests/testNonlinearControlV249.js')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
    (out/'frontend_test.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout.decode('utf8',errors='replace')
    save(out/'report.json',dict(passed=True,fields=len(fields),faces=19752,updates=stages,exactResearchValuesChecked=19752*4*len(stages),
      originalAreaM2=float(sum(g['faceArea'])),physicalQualification=False,solidTwSolved=False,sharedWebGLChanged=False,browserInteractionTested=False))


def main():
    if CLOSED:raise SystemExit('v249 native control export closed; follow NEXT_DECISION')
    g=run.setup_guard()
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',GEOMETRY,AUDIT/'report.json',
      ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json',ROOT/'tests/testNonlinearControlV249.js',BASE/'native01/stages.csv',ROOT/'scripts/terminalCapacityV249.py',DOC/'AUDIT_RECOVERY_PLAN.md']+list(AUDIT.glob('original_faces_*.csv'))
    locks += [p for p in (ROOT/'web').iterdir() if p.name in ('index.html','viewer.js','viewer.css','section.js','export.js')]
    return supervise(g,BASE,'export',[sys.executable,str(Path(__file__).resolve()),'task'],45,1024,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v249 native control export closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':export()
    else:raise SystemExit(0 if main() else 1)
