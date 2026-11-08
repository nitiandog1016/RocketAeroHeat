"""按准确复数布局和真实棱柱拓扑核算NK/SA必要内存，不将下界当作峰值。"""
from pathlib import Path
import argparse,hashlib,json,time,shutil
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v211_validation'
MESH_SHA='ae2f2d7224e44e5cf8f11984f4032248d9c2a12d90860fb92fcf2b5d9ff00ec7'
VOLUME_SHA='22b35b72396c1a7d14195af9bd24d3c6328138484d4009638c1ab48cddec1d6c'
PROFILE='rah-phase-air-complex-nk-sa-ilu0-gmres50-v1'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')

def lower_bound(counts,sizes,include_ghosts=True):
    """模块一：逐项计并存载荷，slice/指针不重复计算其指向的dense数组。"""
    n,g,f,v=(counts[k] for k in ['cells','ghosts','faces','vertices'])
    c=n+g if include_ghosts else n
    cv,cf,fv,fi=(counts[k] for k in ['cellVertexLinks','cellFaceLinks','faceVertexLinks','internalFaces'])
    q=6;s=sizes['Scalar'];vec=sizes['Vector3'];sl=sizes['Slice'];m=50
    assert s==16 and vec==48 and sl==16
    rows=[]
    def add(name,count,unit,formula,source,group):
        rows.append(dict(name=name,count=int(count),unitBytes=int(unit),bytes=int(count*unit),
            formula=formula,source=source,group=group))
    # 原生NK设nftl=2；每个U/dUdt数组仍各分配C*q*nftl，不按理想紧凑布局缩算。
    cell=[('stageArrays',8*q*s,'C*8*q*s','fluidblock:274-283;newtonkrylovsolver:673'),
        ('sourceTerms',q*s,'C*q*s','fluidblock:284'),
        ('flowState',sizes['FlowState'],'C*sizeof(FlowState)','fluidblock:286-289'),
        ('gasDynamic',2*s,'C*2*s','gas_state:52-65'),
        ('gradient',sizes['FlowGradients']+3*s,'C*(FlowGradients+3*s)','fluidblock:290;flowgradients:88'),
        ('gradientWorkspace',sizes['WLSQGradWorkspace'],'C*WLSQGradWorkspace','fluidblock:291'),
        ('interpolationWorkspace',sizes['LSQInterpWorkspace'],'C*LSQInterpWorkspace','ufluidblock:329'),
        ('interpolationGradient',sizes['LSQInterpGradients']+vec+3*s,'C*(LSQInterpGradients+vec+3*s)','ufluidblock:332-333;lsqinterp:294-300'),
        ('cellObject',sizes['FluidFVCell'],'C*FluidFVCell','ufluidblock:349,442'),
        ('stageViews',4*sl,'C*2*2*slice','fluidfvcell:200-226'),
        ('udf',q*s,'C*q*s','fluidfvcell:227'),
        ('geometry',2*vec+2*s,'C*(pos+lengths+area+volume)','fluidblock:268-272'),
        ('savedGradient',sizes['FlowGradients']+3*s,'C*(FlowGradients+3*s)','fluidblock:295-296'),
        ('savedInterpolationGradient',sizes['LSQInterpGradients']+vec+3*s,'C*(LSQInterpGradients+vec+3*s)','ufluidblock:335-336'),
        ('savedSource',q*s,'C*q*s','fluidblock:297'),
        ('perturbFlag',1,'C','fluidblock:298'),
        ('denseCellJacobian',q*q*8,'C*q*q*double','fluidblock:294;coredata:48'),
        ('jacobianRowViews',q*sl,'C*q*slice','fluidfvcell:246-250'),
        ('inverseBlock',sizes['Matrix']+q*q*8,'C*(Matrix+q*q*double)','fluidfvcell:251;bbla:34-37'),
        ('haloSlices',2*sl,'C*2*slice','fluidblock:259-260')]
    for name,unit,formula,ref in cell:add(name,c,unit,formula,ref,'cell')
    face=[('faceObject',sizes['FVInterface'],'F*FVInterface','ufluidblock:341-344'),
        ('faceFlow',sizes['FlowState'],'F*FlowState','fluidblock:333-334'),
        ('faceGradient',sizes['FlowGradients']+3*s,'F*(FlowGradients+3*s)','fluidblock:335-336;flowgradients:88'),
        ('faceWorkspace',sizes['WLSQGradWorkspace'],'F*WLSQGradWorkspace','fluidblock:337-338'),
        ('faceGasDynamic',2*s,'F*2*s','gas_state:52-65'),
        ('faceDiffusion',4*s,'F*4*s','fvinterface:136-140'),
        ('faceVectors',7*vec,'F*7*vec','coredata:104-111;fluidblock:318-323;ufluidblock:338-339'),
        ('faceScalarAndIndices',8+s+2+6*8+q*s+2*8,'F*(index+area+2bool+stencil+flux+f2c)','fluidblock:315-328'),
        ('faceInterpolation',9*s,'F*9*s','fluidblock:330;onedinterp:78-87'),
        ('gridFaceObject',sizes['USGFace'],'F*USGFace','usgrid:154-193'),
        ('faceTopReferences',16,'F*(grid+solver pointers)','ufluidblock:343;usgrid:424-432')]
    for name,unit,formula,ref in face:add(name,f,unit,formula,ref,'face')
    add('vertices',v,sizes['FVVertex']+2*vec+sizes['FlowGradients']+3*s+8+vec+2*sl,
        'V*(FVVertex+2vec+FlowGradients+3s+ptr+gridVec+2slice)','fvvertex:26-57;ufluidblock:307-310,409-415','vertex')
    add('gridCellObject',n,sizes['USGCell'],'N*USGCell','usgrid:218-270','internal')
    add('internalTopArrays',n,8+8+8+8+s+2+5*sl+8,
        'N*(gridPtr+allIdx+nfaces+dt+wallDistance+2bool+5slice+solverPtr)','fluidblock:258-267;ufluidblock:328-331,362-363','internal')
    add('cellVertexConnectivity',cv,32,'CV*(gridIdx+solverPtr+denseIdx+vertexIdx)','usgrid:245-269;ufluidblock:364-368,411','connectivity')
    add('cellFaceConnectivity',cf,12+24+vec,'CF*(gridIdx/sign+solver/denseIdx/signs+faceDistance)','usgrid:245-269;ufluidblock:369-405','connectivity')
    add('faceVertexConnectivity',fv,24,'FV*(gridIdx+solverPtr+vertexIdx)','usgrid:167;ufluidblock:357-360,414','connectivity')
    add('interpolationCloud',n+cf,8,'(N+CF)*ptr','ufluidblock:588-606','connectivity')
    add('denseCloudIndices',cf,8,'CF*index','ufluidblock:597,602','connectivity')
    add('viscousCloud',n+cf,16,'(N+CF)*2ptr','ufluidblock:623-634','connectivity')
    # nearest-neighbour ILU0，三角/棱柱邻接不按壁面积缩算。局部CRS原位分解，只计一次。
    nnz_blocks=n+2*fi
    add('cellResidualStencil',nnz_blocks,8,'(N+2Fi)*ptr','fluidfvcell:1534-1582','nk')
    add('faceResidualStencilSubset',cf,8,'CF*ptr (only own faces; conservative subset)','fluidfvcell:1584-1594','nk')
    add('iluValuesAndColumns',nnz_blocks,q*q*16,'(N+2Fi)*q*q*(double+size_t)','fluidblock:1355-1440;jacobian:78-90;smla:319-342','nk')
    add('iluRowOffsets',n*q+1,8,'(N*q+1)*size_t','jacobian:86','nk')
    add('gmresVectors',n*q,8*(9+m+1),'N*q*double*(9+m+1)','fluidblock:145-149,1849-1862','nk')
    add('gmresSmallDense',1,6*sizes['Matrix']+8*(2*(m+1)*m+4*(m+1)**2+6*(m+1)),
        '6 Matrix objects + global/local H/Q/Gamma/g/h payload','newtonkrylovsolver:1592-1603;fluidblock:1859-1862','nk')
    data=json.loads((ROOT/'data/thermo/phase_aware_dry_air_v1.json').read_text());nt=len(data['temperatureK'])
    add('onePhaseAirTable',1,sizes['PhaseAir']+(7*nt-1)*8,
        'one class + ts/cps/hs/ss/slopes/mus/ks (at least one instance)','rah_phase_aware_air:17-37','gas')
    total=sum(r['bytes'] for r in rows)
    return dict(rows=rows,payloadLowerBoundBytes=total,payloadLowerBoundGiB=total/1024**3,
        counts=counts,cellCapacityCounted=c,includeGhosts=include_ghosts,
        krylovDimension=m,variables=q,flowTimeLevels=2,fullPeakBytes=None,fitProven=False,
        runAuthorized=False,physicalQualification=False,
        exclusions=['GC bins/array metadata/spare capacity/page slack', 'all exterior cells/faces/vertices',
            'additional residual face/ghost stencils and boundary objects', 'temporary Jacobian point graph/IO/JSON/mesh-read buffers',
            'runtime/mappings/stacks/model/config copies/allocator peak'])

def assess(profile,mesh_sha,counts,sizes,memory_mib=8192):
    """模块二：只提供必要条件拒绝；余额为正也绝不发放求解准入。"""
    if profile!=PROFILE or mesh_sha!=MESH_SHA:raise ValueError('资源身份不相容')
    expected=dict(cells=250816,ghosts=17148,faces=635614,vertices=130119,internalFaces=618466,
        cellVertexLinks=1504896,cellFaceLinks=1254080,faceVertexLinks=2283802)
    if counts!=expected:raise ValueError('原区域拓扑计数不相容')
    if type(memory_mib) is not int or memory_mib<=0:raise ValueError('预算应为正整数MiB')
    r=lower_bound(counts,sizes);budget=memory_mib*1024**2
    return dict(status='impossible' if r['payloadLowerBoundBytes']>=budget else 'not-excluded-unproven',
        lowerBoundBytes=r['payloadLowerBoundBytes'],budgetBytes=budget,runAuthorized=False)

def extruded_counts(faces,wedges,node_count):
    """模块三：只数既有32层列的组合拓扑，既不造新网格，也不压平切口。"""
    edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
    unique,number=np.unique(edges,axis=0,return_counts=True);assert np.all((number==1)|(number==2))
    nf=len(faces);nv=len(np.unique(faces));b=int(np.sum(number==1));e=len(unique);layers=32
    keys=lambda x:np.sort(x,axis=1)@np.array([1,node_count,node_count*node_count],np.int64)
    regular=(wedges<33*node_count).all(axis=1)&np.all(np.sort(wedges[:,:3]%node_count,axis=1)==np.sort(wedges[:,3:]%node_count,axis=1),axis=1)
    pick=regular & np.isin(keys(wedges[:,:3]%node_count),keys(faces))
    cols,occ=np.unique(keys(wedges[pick,:3]%node_count),return_counts=True)
    assert len(cols)==nf and np.all(occ==layers)
    assert len(np.unique(wedges[pick]))==33*nv
    count=dict(cells=layers*nf,ghosts=2*nf+layers*b,faces=33*nf+layers*e,
        vertices=33*nv,internalFaces=31*nf+layers*int(np.sum(number==2)),
        cellVertexLinks=6*layers*nf,cellFaceLinks=5*layers*nf,faceVertexLinks=3*33*nf+4*layers*e)
    assert count['faces']==count['internalFaces']+count['ghosts']
    return count,dict(wallFaces=nf,columns=nv,boundaryEdges=b,uniqueEdges=e,layers=32,selectedWedges=int(pick.sum()))

def inspect(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    report=dict(passed=False,engineBuilds=0,CFDRuns=0,volumeMeshesGenerated=0)
    try:
        layout=json.loads((BASE/'layout01/report.json').read_text());assert layout['passed']
        sizes=layout['sizes'];mesh=ROOT/'outputs/v208_validation/layer_mesh02/exact_mesh.npz'
        volume=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat'
        assert sha(mesh)==MESH_SHA and sha(volume)==VOLUME_SHA
        z=np.load(mesh);counts=dict(cells=len(z['cells']),ghosts=len(z['faces'])-len(z['neighbour']),
            faces=len(z['faces']),vertices=len(z['points']),internalFaces=len(z['neighbour']),
            cellVertexLinks=z['cells'].size,cellFaceLinks=5*len(z['cells']),faceVertexLinks=int(z['faceSizes'].sum()))
        # 原cell配对唯一：ILU邻接必须是两侧不同且无重复单元对。
        pair=np.sort(np.column_stack((z['owner'][:counts['internalFaces']],z['neighbour'])),axis=1)
        assert np.all(pair[:,0]!=pair[:,1]) and len(np.unique(pair,axis=0))==counts['internalFaces']
        decision=assess(PROFILE,MESH_SHA,counts,sizes)
        d=loadmat(str(volume),variable_names=['sourceVertices','sourceFaces','wedges'])
        v=d['sourceVertices'];f=d['sourceFaces'].astype(np.int64)-1;w=d['wedges'].astype(np.int64)-1
        assert len(f)==19752
        checked,base_topo=extruded_counts(f[z['originalWallFaceIDs']-1],w,len(v));assert checked==counts
        old=json.loads((ROOT/'docs/validation/v189/region_contract.json').read_text())['local']
        ids=np.flatnonzero(np.min(v[f,0],axis=1)>=old['extensionTestCutXM']-1e-10)
        assert len(ids)==old['extensionWallFaceCount']
        extension,ext_topo=extruded_counts(f[ids],w,len(v))
        area=float(np.linalg.norm(np.cross(v[f[ids,1]]-v[f[ids,0]],v[f[ids,2]]-v[f[ids,0]]),axis=1).sum()/2)
        assert abs(area-old['extensionWallAreaM2'])<1e-12
        results=dict(baseLayer=lower_bound(counts,sizes),baseWithoutGhosts=lower_bound(counts,sizes,False),
            extensionLayer=lower_bound(extension,sizes),extensionWithoutGhosts=lower_bound(extension,sizes,False))
        for name,value in results.items():save(out/(name+'.json'),value)
        extension_impossible=results['extensionWithoutGhosts']['payloadLowerBoundBytes']>=8192*1024**2
        report.update(passed=True,profile=PROFILE,sizes=sizes,baseDecision=decision,counts=counts,
            baseTopology=base_topo,extensionTopology=ext_topo,extensionCounts=extension,
            extensionOriginalFaceIDs=(ids+1).tolist(),extensionWallAreaM2=area,
            lowerBoundsBytes={k:r['payloadLowerBoundBytes'] for k,r in results.items()},
            extensionNecessaryPayloadExceeds8GiB=extension_impossible,
            completeOriginalPlanResourceAccepted=False,physicalQualification=False,runAuthorized=False,
            note='Without-ghost cases omit every ghost allocation: permanent internal geometry/state only; outer grid still omitted.',
            inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [mesh,volume,BASE/'layout01/report.json',ROOT/'docs/validation/v211/PLAN.md',Path(__file__)]})
        assert all(sha(ROOT/p)==h for p,h in report['inputSHA256'].items())
    except Exception as exc:
        report['error']=repr(exc);raise
    finally:
        report['seconds']=time.monotonic()-start;save(out/'report.json',report)
    assert report['seconds']<60
    print(json.dumps({k:v for k,v in report.items() if k not in ['extensionOriginalFaceIDs','inputSHA256','sizes']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);inspect(p.parse_args().out)
