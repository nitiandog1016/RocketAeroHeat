"""从同一次组装的边、边界、源、储能和隐式右端复核轴对称整域收支。"""
import argparse
from math import fsum
from pathlib import Path
import numpy as np
from run_ledger import BASE,DOC,PARENT,ROOT,FRAME_HASH,KEYS,numbers,read,save,sha,terminal

STAGES=("start","transport","viscous","source","time","weak","fixed","wall","end","rhs")
OLD=DOC.parent/"v369/evidence"
COLUMNS=109


def total(rows):
    return np.array([fsum(rows[:,j]) for j in range(rows.shape[1])])


def save_csv(path,header,rows):
    with path.open("x",encoding="ascii",newline="\n") as f:
        np.savetxt(f,rows,delimiter=",",fmt="%.17g",header=",".join(header),comments="")


def axial_uniform_moment(folder,xy,i,j,radius,mass_flux,enthalpy):
    """用均匀轴向通量分离中点半径分摊本身的局部正负项；采用原四边形中位对偶面。"""
    with (folder/"mesh.su2").open(encoding="ascii") as f:
        assert f.readline().strip()=="NDIME= 2"
        count=int(f.readline().split("=")[1])
        cells=np.array([[int(v) for v in f.readline().split()] for _ in range(count)])
    rows=len(xy)//COLUMNS
    assert cells.shape==((rows-1)*(COLUMNS-1),5) and (cells[:,0]==9).all()
    corners=cells[:,1:]
    centers=xy[corners].mean(axis=1)
    ids={tuple(sorted((int(a),int(b)))):e for e,(a,b) in enumerate(zip(i,j))}
    normals=np.zeros((len(i),2))
    for a,b in zip(corners.T,np.roll(corners,-1,axis=1).T):
        edge=np.array([ids[(min(left,right),max(left,right))] for left,right in zip(a,b)])
        dr=centers-.5*(xy[a]+xy[b])
        half=np.column_stack((dr[:,1],-dr[:,0]))
        sign=np.sign(np.sum(half*(xy[j[edge]]-xy[i[edge]]),axis=1))
        assert np.all(sign!=0)
        np.add.at(normals,edge,half*sign[:,None])
    nodal=np.zeros_like(xy)
    np.add.at(nodal,i,normals);np.add.at(nodal,j,-normals)
    interior=nodal.reshape(rows,COLUMNS,2)[1:-1,1:-1]
    assert np.max(abs(interior))<1e-12
    uniform=.5*(radius[i]-radius[j])*normals[:,0]*mass_flux
    local=np.zeros(len(xy))
    np.add.at(local,i,uniform);np.add.at(local,j,uniform)
    # 单独核对两个绕轴的边界，均匀轴向输运在各自的闭合曲面进出相抵。
    for band,orientation in ((slice(0,109),1),(slice(-109,None),-1)):
        dr=np.diff(xy[band,1])
        nx=.5*np.r_[dr[0],dr[:-1]+dr[1:],dr[-1]]
        np.testing.assert_allclose(nodal[band,0]*radius[band],orientation*nx*radius[band],atol=1e-12)
        assert abs(fsum(nx*radius[band]))<1e-11
    assert abs(fsum(local)*enthalpy)<1e-5
    return np.column_stack((local,local*enthalpy))


def wall_regions(xy,wall_q,wall_rows,surface):
    """按相邻半边法向在节点对偶法向上的投影分摊接头节点，六区之和等于原生壁面项。"""
    edge=np.diff(xy,axis=0)
    vectors=np.column_stack((edge[:,1],-edge[:,0]))
    half=vectors/2
    dual=np.vstack((half[0],half[1:]+half[:-1],half[-1]))
    area=np.linalg.norm(dual,axis=1)
    unit=dual/area[:,None]
    ratio=2*np.pi*xy[:,1]
    np.testing.assert_allclose(wall_q*area*ratio,wall_rows[:,1],atol=1e-7,rtol=2e-8)
    share_i=np.sum(half*unit[:-1],axis=1)/area[:-1]
    share_j=np.sum(half*unit[1:],axis=1)/area[1:]
    np.testing.assert_allclose(np.r_[share_i,0]+np.r_[0,share_j],1,atol=1e-12)
    segments=share_i*wall_rows[:-1,1]+share_j*wall_rows[1:,1]
    model=read(ROOT/"cases/simple_launcher/geometry.json")
    mesh=np.dtype([("normal","<f4",3),("v","<f4",(3,3)),("attr","<u2")])
    stl=ROOT/"cases/simple_launcher"/model["file"]
    assert sha(stl)==model["fileSHA256"]
    tri=np.fromfile(stl,dtype=mesh,offset=84)["v"]
    mid=model["lengthM"]-(xy[:-1,0]+xy[1:,0])/2
    used=np.zeros(len(mid),int);regions={}
    for r in model["regions"]:
        name=r["name"]
        x=tri[r["firstFace"]-1:r["lastFace"],:,0]
        mask=(abs(mid)<1e-12) if name=="base" else (mid>float(x.min()))&(mid<float(x.max()))
        assert mask.any()
        used+=mask
        regions[name]=dict(nativeWeakWallW=fsum(segments[mask]),originalSurfaceW=surface["regions"][name]["researchIntoWallW"],
                           faces=r["lastFace"]-r["firstFace"]+1)
    np.testing.assert_array_equal(used,np.ones_like(used))
    np.testing.assert_allclose(sum(v["nativeWeakWallW"] for v in regions.values()),fsum(wall_rows[:,1]),atol=1e-8)
    return regions


def export(folder,out):
    folder=folder.resolve();folder.relative_to(BASE)
    surface=read(OLD/"surface.json")
    assert surface["frameSHA256"]==FRAME_HASH
    audit_state(folder,out,PARENT,surface,read(OLD/"runs.json")["half1"]["fieldBalance"],read(OLD/"far.json")["characteristicNodal"])


def audit_state(folder,out,parent,model,old,old_char):
    """在保持 109 个子午列的网格上核对原生项，每份外部对照取自对应网格的步末。"""
    proc=terminal(folder)
    assert all(sha(Path(k))==v for k,v in proc["inputSHA256"].items())
    replay=read(folder/"replay.json")
    assert replay["frameSHA256"]==model["frameSHA256"]==sha(parent/"restart_00001.csv")
    assert replay["previousSHA256"]==sha(parent/"seed_00000.csv") and replay["meshSHA256"]==sha(parent/"mesh.su2")
    assert all(sha(folder/name)==v for name,v in replay["stageFiles"].items())
    def array(name,width):
        a=np.fromfile(folder/(name+".bin"),dtype="<f8").reshape(-1,width)
        assert np.isfinite(a).all()
        return a
    state=array("state",11)
    n=len(state);rows=n//COLUMNS
    assert n==rows*COLUMNS
    frame=numbers(parent/"restart_00001.csv")
    for values,data in ((state[:,3:7],frame),(state[:,7:],numbers(parent/"seed_00000.csv"))):
        np.testing.assert_array_equal(data["PointID"],np.arange(n))
        np.testing.assert_array_equal(values,np.column_stack([data[k] for k in KEYS]))
    np.testing.assert_allclose(state[:,:2],np.column_stack((frame["x"],frame["y"])),atol=3e-15,rtol=1e-14)
    rho,v,p,specific_e,ref_length,ref_time,dt=array("scales",7)[0]
    np.testing.assert_allclose([rho*v*v,specific_e,ref_time*v,ref_length,dt*ref_time],[p,v*v,1,1,6.25e-6],rtol=1e-12)
    scale=np.array([rho*v,rho*specific_e*v])
    radius=2*np.pi*state[:,1]
    edges=array("edges",6);i,j=edges[:,:2].astype(int).T
    assert (edges[:,:2]==np.column_stack((i,j))).all() and min(i)>=0 and min(j)>=0 and max(i)<n and max(j)<n
    assert len(edges)==rows*108+(rows-1)*109 and (i!=j).all()
    pairs=np.sort(np.column_stack((i,j)),axis=1)
    assert len(np.unique(pairs,axis=0))==len(pairs) and np.isin(pairs[:,1]-pairs[:,0],[1,109]).all()
    st={k:array(k,4)[:,[0,3]] for k in STAGES}
    flux=edges[:,4:]
    edge_sum=np.zeros((n,2))
    np.add.at(edge_sum,i,flux);np.add.at(edge_sum,j,-flux)
    raw_far=array("far",10);active=raw_far[:,0]==1
    np.testing.assert_array_equal(np.flatnonzero(active),np.arange(n-109,n))
    assert not raw_far[~active].any()
    far_c=raw_far[:,3:5];far_v=raw_far[:,8:10]
    source=st["source"]-st["viscous"]
    storage=(state[:,[3,6]]-state[:,[7,10]])*state[:,2,None]/dt
    hot=st["wall"]-st["fixed"];axis=st["end"]-st["wall"]
    exact=edge_sum+source+storage+far_c+far_v+hot+axis
    errors={k:float(np.max(abs(val))) for k,val in dict(initial=st["start"],internal=edge_sum-st["transport"],
         separateViscous=st["viscous"]-st["transport"],time=storage-(st["time"]-st["source"]),
         far=far_c+far_v-(st["weak"]-st["time"]),imposed=st["fixed"]-st["weak"],
         nonwallHot=hot[109:],rhsSign=st["end"]+st["rhs"],localClosure=exact+st["rhs"]).items()}
    assert max(errors.values())<1e-12,errors
    def weighted(raw):
        return raw*radius[:,None]*scale
    # 用同一内部边中点半径作为局部守恒基准。两端余下的半径差项相同。
    moment=.5*(radius[i]-radius[j])[:,None]*flux*scale
    local_moment=np.zeros((n,2))
    np.add.at(local_moment,i,moment);np.add.at(local_moment,j,moment)
    np.testing.assert_allclose(total(local_moment),total(weighted(edge_sum)),atol=1e-6,rtol=1e-12)
    terms={"storage":weighted(storage),"farConvective":weighted(far_c),"farViscous":weighted(far_v),
           "wall":weighted(hot),"symmetry":weighted(axis),"internalRadiusMoment":local_moment,
           "axisymmetricSource":weighted(source),"physicalResidual":weighted(-st["rhs"])}
    sums={k:total(z) for k,z in terms.items()}
    native=sums["storage"]+sums["farConvective"]+sums["farViscous"]+sums["wall"]+sums["symmetry"]
    geom=sums["internalRadiusMoment"]+sums["axisymmetricSource"]
    np.testing.assert_allclose(native+geom,sums["physicalResidual"],atol=2e-6,rtol=0)
    prior={"storage":[old["massRateKgS"],old["storedTotalEnergyRateW"]],
           "far":[old["farOutMassKgS"],old["farOutTotalEnergyW"]],"wall":[0,old["wallIntoSolidW"]]}
    deltas={"storage":prior["storage"]-sums["storage"],"far":prior["far"]-sums["farConvective"]-sums["farViscous"],
            "wall":prior["wall"]-sums["wall"],"symmetry":-sums["symmetry"]}
    old_balance=np.array([old["massImbalanceKgS"],old["totalEnergyImbalanceW"]])
    np.testing.assert_allclose(sum(deltas.values())+native,old_balance,atol=1e-6,rtol=0)
    local_geom=local_moment+terms["axisymmetricSource"]
    # r>0 处的轴对称对流源可由守恒量直接复算，r=0 的节点在旋转积分中权重为零。
    density,mx,mr,energy=state[:,3:7].T
    pressure=.4*(energy-.5*(mx*mx+mr*mr)/density)
    axis_adv=2*np.pi*state[:,2,None]*np.column_stack((mr,(energy+pressure)*mr/density))*scale
    axis_adv[state[:,1]==0]=0
    np.testing.assert_allclose(axis_adv[:,0],terms["axisymmetricSource"][:,0],atol=1e-11,rtol=1e-10)
    edge_adv=total((radius[i]-radius[j])[:,None]*edges[:,2:4]*scale)
    geom_adv=edge_adv+total(axis_adv)
    href=old["freeTotalEnthalpyJkg"]
    def excess(row):
        return float(row[1]-href*row[0])
    layer_geom=local_geom.reshape(rows,109,2).sum(axis=1)
    condition=read(ROOT/"outputs/v365_validation/warm2/condition.json")
    uniform=axial_uniform_moment(folder,state[:,:2],i,j,radius,rho*condition["speedMps"],href)
    variation=local_geom-uniform
    corrected_layers=variation.reshape(rows,109,2).sum(axis=1)
    far_data=raw_far[active]
    weight=radius[active]*np.linalg.norm(far_data[:,1:3],axis=1)
    char_mass=far_data[:,5]*far_data[:,7]*weight
    char_native=np.array([fsum(char_mass)*scale[0],fsum(char_mass*far_data[:,6])*scale[1]])
    np.testing.assert_allclose(char_native,[old_char["netMassKgS"],old_char["netAdvectiveW"]],rtol=2e-9,atol=5e-6)
    regions=wall_regions(state[:109,:2],frame["Heat_Flux"][:109],terms["wall"][:109],model)
    out.mkdir(parents=True,exist_ok=False)
    # 发布分层、分列积分及真实远界、热壁的小表；整域逐点阶段和每条内部边保留在 E 盘。
    summary=dict(frameSHA256=model["frameSHA256"],surfaceSHA256=model["surfaceSHA256"],faces=model["faces"],
         mach=model["mach"],altitudeM=model["altitudeM"],wallGivenK=model["wallGivenK"],elapsedUs=model["elapsedFromNumericalSeedUs"],
         units=["kg/s","W"],nativeTerms={k:v.tolist() for k,v in sums.items()},nativeBoundaryPlusStorage=native.tolist(),
         radiusMomentPlusSource=geom.tolist(),savedFieldBalance=old_balance.tolist(),
         savedMinusNative={k:v.tolist() for k,v in deltas.items()},characteristicPhysicalFlux=char_native.tolist(),
         referenceEnthalpyJkg=href,nativeExcessImbalanceW=excess(native),geometryExcessW=excess(geom),
         geometricAdvection=geom_adv.tolist(),geometricDiffusion=(geom-geom_adv).tolist(),
         lastTwoLayersGeometry=total(layer_geom[-2:]).tolist(),lastTwoLayersExcessW=excess(total(layer_geom[-2:])),
         lastTwoLayersAfterUniformSubtraction=total(corrected_layers[-2:]).tolist(),
         regions=regions,checksNativeND=errors,replay=replay,wholeBodyPowerW=None)
    length=condition["lengthM"]
    for label,data in (("mass",variation[:,0]),("energy",variation[:,1]),("excessEnergy",variation[:,1]-href*variation[:,0])):
        idx=int(np.argmax(abs(data)))
        summary[label+"MaxAfterUniformSubtraction"]=dict(point=idx,layer=idx//109,column=idx%109,
             modelXM=float(length-state[idx,0]),radiusM=float(state[idx,1]),signedValue=float(data[idx]))
    for dep,name in ((0,"columns"),(1,"layers")):
        tables=[a.reshape(rows,109,2).sum(axis=dep) for a in (*terms.values(),local_geom,uniform,variation)]
        header=[k+"_"+u for k in (*terms,"combinedGeometry","uniformGeometry","geometryMinusUniform") for u in ("kg_s","W")]
        save_csv(out/(name+".csv"),header,np.column_stack(tables))
    save_csv(out/"far.csv",["point","radiusM","convMassKgS","convEnergyW","viscMassKgS","viscEnergyW"],
             np.column_stack((np.flatnonzero(active),state[active,1],terms["farConvective"][active],terms["farViscous"][active])))
    save_csv(out/"wall.csv",["point","solverX","radiusM","fluxWm2","weakMassKgS","weakPowerW"],
             np.column_stack((np.arange(109),state[:109,:2],frame["Heat_Flux"][:109],terms["wall"][:109])))
    save(out/"balance.json",summary)
    print("原生边界与储能",native,"半径与轴对称项",geom,"保存场与原生的差",deltas,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("folder",type=Path);p.add_argument("out",type=Path)
    args=p.parse_args();export(args.folder,args.out)
