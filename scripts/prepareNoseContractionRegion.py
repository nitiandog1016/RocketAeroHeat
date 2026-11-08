"""准备真实鼻部—收缩区域的封闭流体边界；保留原面，不生成实体端盖。"""
from pathlib import Path
import json,time,shutil
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from phaseAwareViscousAir import PhaseAwareViscousAir

BASE=ROOT/'outputs/v207_validation'
def edges(faces):
    e=np.vstack([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]);u,c=np.unique(np.sort(e,axis=1),axis=0,return_counts=True)
    return e,u,c
def triangulate_annulus(inner,outer,points):
    """两个按极角排序的原环/外环进行拉链连接，保持原切口边。"""
    angles=lambda ids:np.mod(np.arctan2(points[ids,2],points[ids,1]),2*np.pi)
    ai=angles(inner);ao=angles(outer);i=j=0;tri=[];ni=len(inner);no=len(outer)
    while i<ni or j<no:
        ti=ai[(i+1)%ni]+(2*np.pi if i+1>=ni else 0) if i<ni else np.inf
        to=ao[(j+1)%no]+(2*np.pi if j+1>=no else 0) if j<no else np.inf
        if ti<to:tri.append([inner[i%ni],inner[(i+1)%ni],outer[j%no]]);i+=1
        else:tri.append([inner[i%ni],outer[(j+1)%no],outer[j%no]]);j+=1
    tri=np.array(tri,int);xyz=points[tri];bad=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])[:,0]>0
    tri[bad]=tri[bad][:,[0,2,1]]
    return tri

def main():
    start=time.monotonic();out=BASE/'region01';out.mkdir(exist_ok=False)
    if not json.loads((BASE/'thermal_audit03/report.json').read_text())['passed']:raise ValueError('热接口控制未过')
    source=ROOT/'outputs/v194_validation/rocket_design03/surface.json';selection=ROOT/'docs/validation/v189/evidence/selection01/region_geometry.json'
    s=json.loads(source.read_text());r=json.loads(selection.read_text());P=np.array(s['vertices']);F=np.array(s['facesOneBased'],int)-1
    ids=np.array(r['localWall']['faceIDs'],int);wall=F[ids-1];cut=r['cutXM'];e,u,c=edges(wall);rim_edges=u[c==1];rim=np.unique(rim_edges)
    if len(F)!=19752 or len(ids)!=7838 or len(rim)!=46 or len(rim_edges)!=46 or np.any(c>2):raise ValueError('原面或切口拓扑改变')
    if np.max(abs(P[rim,0]-cut))>1e-12:raise ValueError('切口不在原环平面')
    theta=np.mod(np.arctan2(P[rim,2],P[rim,1]),2*np.pi);rim=rim[np.argsort(theta)]
    expected=np.sort(np.column_stack([rim,np.roll(rim,-1)]),axis=1)
    if set(map(tuple,expected))!=set(map(tuple,rim_edges)):raise ValueError('原环不能按极角直接连接')
    xyz=P[F];sf=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])/2;areas=np.linalg.norm(sf,axis=1)
    if abs(areas.sum()-1.9646292565715553)>1e-13 or abs(areas[ids-1].sum()-r['localWall']['areaM2'])>1e-13:raise ValueError('原面面积改变')
    full_volume=float(np.einsum('ij,ij->i',xyz[:,0],sf).sum()/3)
    if full_volume==0:raise ValueError('原壁方向未定义')
    wall_fluid=wall[:,[0,2,1]] if full_volume>0 else wall.copy()
    local_points=np.unique(wall);radius=float(np.linalg.norm(P[local_points,1:],axis=1).max());tip=float(P[local_points,0].max())
    front=tip+2*radius;rf=4*radius;slope=np.tan(np.deg2rad(40));ra=rf+(front-cut)*slope;n=128
    angle=np.arange(n)*2*np.pi/n
    circle=lambda x,rad:np.column_stack([np.full(n,x),rad*np.cos(angle),rad*np.sin(angle)])
    pf=circle(front,rf);pa=circle(cut,ra);all_points=np.vstack([P,pf,pa,[[front,0,0]]]);front_ids=np.arange(len(P),len(P)+n);aft_ids=front_ids+n;center=len(all_points)-1
    inlet=np.column_stack([np.full(n,center),front_ids,np.roll(front_ids,-1)])
    sides=[]
    for i in range(n):
        j=(i+1)%n;sides.extend([[front_ids[i],aft_ids[i],aft_ids[j]],[front_ids[i],aft_ids[j],front_ids[j]]])
    sides=np.array(sides,int);q=all_points[sides];norm=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);centers=q.mean(axis=1)
    flip=np.einsum('ij,ij->i',norm[:,1:],centers[:,1:])<0;sides[flip]=sides[flip][:,[0,2,1]]
    outlet=triangulate_annulus(rim,aft_ids,all_points)
    faces=np.vstack([wall_fluid,inlet,sides,outlet]);patch=np.r_[np.zeros(len(wall),int),np.ones(len(inlet),int),np.full(len(sides),2),np.full(len(outlet),3)]
    directed,unique,count=edges(faces)
    if np.any(count!=2):raise ValueError('流体边界没有闭合两流形')
    signed={}
    for a,b in directed:
        key=(min(a,b),max(a,b));signed[key]=signed.get(key,0)+(1 if a<b else -1)
    if any(signed.values()):raise ValueError('流体边界方向不相容')
    q=all_points[faces];sf=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])/2;ar=np.linalg.norm(sf,axis=1);volume=float(np.einsum('ij,ij->i',q[:,0],sf).sum()/3)
    closure=float(np.linalg.norm(sf.sum(axis=0))/ar.sum())
    if ar.min()<=0 or volume<=0 or closure>1e-13:raise ValueError('面积/体积/向量闭合失败')
    direction=np.array([-np.cos(np.deg2rad(20)),0,np.sin(np.deg2rad(20))]);normal=sf/ar[:,None];mn=5*normal@direction
    if np.any(mn[(patch==1)|(patch==2)]>=-1) or np.any(mn[patch==3]<=1):raise ValueError('均匀来流法向特征不符设计')
    gap=rf+(front-P[local_points,0])*slope-np.linalg.norm(P[local_points,1:],axis=1)
    if gap.min()<=0:raise ValueError('原壁穿过外域')
    used=np.unique(faces);reindex=np.full(len(all_points),-1,int);reindex[used]=np.arange(len(used));mapped=reindex[faces]
    original_vertex=np.where(used<len(P),used+1,0);original_face=np.r_[ids,np.zeros(len(faces)-len(ids),int)]
    np.savez_compressed(out/'boundary.npz',points=all_points[used],faces=mapped,patch=patch,originalVertexIDs=original_vertex,originalFaceIDs=original_face)
    # 模块二：保留实际原面索引及待验证项，不把这张边界面当已通过的体网格。
    contract=dict(schema='RAH.RegionalViscousBoundary.v1',candidate='nose-contraction-of13-phase-aware-v1',
        status='closed-surface-prepared-volume-and-outlet-not-admitted',executionAuthorized=False,
        gas=PhaseAwareViscousAir().viscous_contract(),wallOriginalFaceIDs=ids.tolist(),wallFaceCount=len(ids),wallAreaM2=float(areas[ids-1].sum()),
        fullOriginalFaceCount=len(F),fullOriginalAreaM2=float(areas.sum()),remainingOriginalFaceCount=len(F)-len(ids),
        remainingOriginalAreaM2=float(areas.sum()-areas[ids-1].sum()),cutRingOriginalVertexIDs=(rim+1).tolist(),
        noseXM=tip,cutXM=cut,inletXM=front,inletRadiusM=rf,outletRadiusM=ra,outerHalfAngleDeg=40,outerCircumferentialPoints=n,
        volumeM3=volume,boundaryTriangles=len(faces),boundaryVertices=len(used),vectorAreaClosure=closure,
        uniformFarfieldNormalMachRange=[float(mn[(patch==1)|(patch==2)].min()),float(mn[(patch==1)|(patch==2)].max())],
        uniformOutletNormalMachRange=[float(mn[patch==3].min()),float(mn[patch==3].max())],minimumRadialSeparationM=float(gap.min()),
        artificialPhysicalCap=False,wallBoundary='no-slip; prescribed Tw300K; viscous q from actual temperature gradient',
        farfieldBoundary='fixed physical freestream only if actual inflow characteristics and undisturbed outer-domain tests hold',
        outletBoundary='not frozen: thermally consistent outgoing acoustic condition, no accepted reverse flow; cutoff/extension test required',
        turbulence='laminar verification passed only; actual RANS/transition choice and inlet turbulence not yet qualified',
        nearWall='retain required wall layers and original 19-layer asset; real anisotropic mesh quality/y+ must be independently accepted',
        designComparisons=dict(extendedCutXM=r['extensionTestCutXM'],extendedWallFaces=r['expandedWall']['faceCount'] if 'expandedWall' in r else 9678,outerRadiusFactor=1.5,performed=False),
        missingBeforeCFD=['native variable-cp outlet contract','actual volume/near-wall grid and resource bounds','turbulence/transition scope','public heat benchmark and numerical convergence'],
        sourceSHA256={str(source.relative_to(ROOT)):sha(source),str(selection.relative_to(ROOT)):sha(selection)},
        dataSHA256=sha(out/'boundary.npz'),physicalQualification=False,newQualifiedHeatAreaM2=0)
    save(out/'contract.json',contract);shutil.copyfile(__file__,out/'executed_entry.py')
    save(out/'report.json',dict(passed=True,boundaryClosed=True,originalWallConnectionsPreserved=True,originalWallAreaPreserved=True,
        faceCount=len(faces),vertexCount=len(used),volumeM3=volume,wallAreaM2=contract['wallAreaM2'],farfieldNormalMach=contract['uniformFarfieldNormalMachRange'],
        CFDRuns=0,volumeMeshGenerated=False,physicalQualification=False,elapsedSeconds=time.monotonic()-start))
    if time.monotonic()-start>60 or sum(p.stat().st_size for p in out.iterdir() if p.is_file())>64*1024**2:raise ValueError('区域准备预算超限')
    print(json.dumps(dict(passed=True,originalFaces=len(ids),boundaryTriangles=len(faces),volume=volume,farfieldNormalMach=contract['uniformFarfieldNormalMachRange'],seconds=time.monotonic()-start)))

if __name__=='__main__':main()
