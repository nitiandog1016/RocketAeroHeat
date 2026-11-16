"""模块职责：把保存的三维流管接入只读科研页和可分享科研图，不重新计算轨迹。"""
from pathlib import Path
import json,shutil,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from runSU2RocketPhaseV222 import ROOT,read,save,sha

BASE=ROOT/'outputs/v241_validation';SOURCE=BASE/'streamtubes02'


def export():
    assert read(BASE/'audit02/report.json')['passed'];m=read(SOURCE/'manifest.json');a=read(BASE/'audit02/report.json')
    g=np.load(str(ROOT/'outputs/v230_validation/audit02/geometry.npz'));out=BASE/'workbench01';out.mkdir(exist_ok=False)
    ids=np.loadtxt(str(ROOT/'outputs/v240_validation/run01/original_faces.csv'),delimiter=',',skiprows=1)[:,0].astype(int)
    payload=dict(schema='RAH.StreamtubeWorkbench.v1',status='unqualified-research',physicalQualification=False,edgeAccepted=False,
      newQualifiedAreaM2=0,sourceUpdate=200,sourceSHA256=m['sourceSHA256'],vertices=g['vertices'].tolist(),faces=g['faces'].tolist(),
      originalAreaM2=float(sum(g['area'])),sourceFaceIDs=ids.tolist(),columns=None,tubes=[],
      sourceManifestSHA256=sha(SOURCE/'manifest.json'),auditSHA256=sha(BASE/'audit02/report.json'),
      qStatus='unknown',givenWallTemperatureK=300,solidWallTemperatureStatus='unknown',exactODECurveInsideProof=False)
    count=0
    for info in m['tubes']:
        path=BASE/'audit02'/('tube%02d.csv'%info['id']);columns=path.read_text(encoding='utf8').splitlines()[0].split(',')
        data=np.loadtxt(str(path),delimiter=',',skiprows=1);payload['columns']=columns
        rows=[{k:float(v) for k,v in zip(columns,row)} for row in data]
        for row in rows:
            if row['metricValid']==0:
                for key in ('hBeta','hEta','signedArea','massFluxRatio','HRelativeDrift','entropyDrift'):row[key]=None
        z=np.load(str(SOURCE/info['file']));stops=[r for r in a['boundaryStops'] if r['file']==info['file']]
        summary=' / '.join('%s %d条'%({'wall':'壁面','outlet':'出口','farfield':'远场'}[tag],sum(r['marker']==tag for r in stops)) for tag in sorted(set(r['marker'] for r in stops)))
        payload['tubes'].append(dict(info,rows=rows,rays=[z['ray%d_xyz'%j].tolist() for j in range(9)],
          stopSummary=summary,confirmedStops=stops,CSVSourceSHA256=sha(path)))
        count+=data.size
    save(out/'streamtubes.json',payload)
    (out/'streamtubes-data.js').write_text('window.STREAMTUBES='+json.dumps(payload,ensure_ascii=False,allow_nan=False,separators=(',',':'))+';\n',encoding='utf8')
    shutil.copyfile(ROOT/'web/streamtube-research-v241.html',out/'index.html')
    restored=read(out/'streamtubes.json');assert restored['vertices']==g['vertices'].tolist() and restored['faces']==g['faces'].tolist()
    for t in restored['tubes']:
        data=np.loadtxt(str(BASE/'audit02'/('tube%02d.csv'%t['id'])),delimiter=',',skiprows=1)
        assert np.array_equal(np.array([[r[k] for k in restored['columns']] for r in t['rows']]),data)
    test=subprocess.run(['node',str(ROOT/'tests/testStreamtubeWorkbenchV241.js'),str(out/'streamtubes.json')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(test.stdout);assert test.returncode==0,test.stdout.decode('utf8',errors='replace')
    # 模块二：同一保存数据的科研图；全箭保留、前体放大、H/S和质量变化分开展示。
    fig,axs=plt.subplots(3,2,figsize=(13,10));fig.subplots_adjust(left=.08,right=.96,bottom=.12,top=.88,hspace=.46,wspace=.28)
    colors=['#1f6d8c','#ba6735','#648044','#7d5c98'];xyz=g['vertices'];faces=g['faces'];sel=np.isin(np.arange(len(faces))+1,ids)
    right=np.array([1.,.08,0]);right/=np.linalg.norm(right);up=np.array([0.,-.2,1.]);up-=right*np.dot(up,right);up/=np.linalg.norm(up)
    proj=np.column_stack([right,up]);points=xyz[faces]@proj;order=np.argsort(xyz[faces].mean(axis=1)@np.cross(right,up))
    for col in range(2):
        ax=axs[0,col];which=order if col==0 else order[sel[order]]
        ax.add_collection(PolyCollection(points[which],facecolors=['#d4dfe6' if sel[i] else '#e5e8eb' for i in which],edgecolors='none'))
        for info,color in zip(restored['tubes'],colors):
            ray=np.asarray(info['rays'][0])@proj;ax.plot(ray[:,0],ray[:,1],color=color,lw=1.1,label='Tube '+str(info['id']))
            ax.plot(ray[0,0],ray[0,1],'o',color=color,ms=3)
            pt=np.asarray(info['confirmedStops'][0]['xyz'])@proj;ax.plot(pt[0],pt[1],'x',color=color,ms=4)
        ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.set_xlabel('Projected X [m]');ax.set_ylabel('Projected transverse [m]')
        ax.set_title('Full original surface' if col==0 else 'Forebody detail: seeds and boundary stops',loc='left',fontsize=10)
    for ax,key,label,scale in [(axs[1,0],'HRelativeDrift','Total H drift from inlet [%]',100),
      (axs[1,1],'entropyDrift','Entropy change [J/(kg K)]',1),(axs[2,0],'massFluxRatio','Local mass-flux ratio (not integral)',1),
      (axs[2,1],'hBeta','Projected transverse width ratio',1)]:
        for t,color in zip(restored['tubes'],colors):
            ax.plot([r['sM'] for r in t['rows']],[r[key]*scale for r in t['rows']],color=color,lw=1.2,label='Tube '+str(t['id']))
        ax.set_xlabel('Arc length s [m]');ax.set_ylabel(label);ax.grid(alpha=.18)
    axs[1,0].legend(fontsize=8,ncol=4);fig.suptitle('RocketAeroHeat | M5 / alpha 20 deg | internal 3-D volume streamtubes\nSaved Euler update 200: research only, NO accepted edge or new wall heat',fontsize=12,y=.965)
    fig.text(.08,.035,'Original 19,752 faces / 1.964629 m² retained. Curves show only the common nine-ray prefix.\nNo projection onto wall, no extrapolation through solid. Crosses mark recorded boundary intersections; q and solid Tw remain unknown.',fontsize=9,color='#48545e')
    target=out/'streamtubes.png';fig.savefig(str(target),dpi=160);plt.close(fig)
    save(out/'report.json',dict(passed=True,tubes=4,originalFaces=len(faces),originalAreaM2=payload['originalAreaM2'],
      totalCommonStations=sum(t['commonStations'] for t in restored['tubes']),exactCSVRoundtripNumbers=count,
      actualFrontendPureFunctionsPassed=True,curveUnitsAndResearchIdentityPreserved=True,
      browserInteractionTested=False,sharedWebGLChanged=False,newScientificPage=True,CFDUpdates=0,trajectoryReintegrations=0,
      physicalQualification=False,newQualifiedAreaM2=0,payloadSHA256=sha(out/'streamtubes.json'),scientificImageSHA256=sha(target)))


if __name__=='__main__':export()
