"""用发布的小表重算整域原生收支、分区热壁项及参考焓质量换算。"""
from math import fsum
from pathlib import Path
import json
import numpy as np

DOC=Path(__file__).resolve().parent
DATA=DOC/"evidence"


def read(path):
    return json.loads(path.read_bytes())


def verify():
    s=read(DATA/"balance.json")
    old=read(DOC.parent/"v369/evidence/runs.json")["half1"]
    surface=read(DOC.parent/"v369/evidence/surface.json")
    assert old["passed"] and old["frameSHA256"]==s["frameSHA256"]==surface["frameSHA256"]
    assert s["replay"]["previousSHA256"]==old["sourceSHA256"] and s["replay"]["meshSHA256"]==old["meshSHA256"]
    assert s["faces"]==surface["faces"]==9216 and s["elapsedUs"]==old["nominalEndUs"]==6.25
    assert max(s["checksNativeND"].values())<1e-12 and s["wholeBodyPowerW"] is None
    terms={k:np.array(v) for k,v in s["nativeTerms"].items()}
    def pair(a,k):
        return np.column_stack((a[k+"_kg_s"],a[k+"_W"]))
    for name,count in (("layers",89),("columns",109)):
        tab=np.genfromtxt(DATA/(name+".csv"),delimiter=",",names=True)
        assert len(tab)==count
        for k,v in terms.items():
            z=pair(tab,k)
            np.testing.assert_allclose(fsum(z[:,0]),v[0],atol=1e-10,rtol=0)
            np.testing.assert_allclose(fsum(z[:,1]),v[1],atol=2e-6,rtol=0)
        g=pair(tab,"combinedGeometry")
        np.testing.assert_allclose(g,pair(tab,"internalRadiusMoment")+pair(tab,"axisymmetricSource"),atol=2e-6,rtol=0)
        np.testing.assert_allclose(g-pair(tab,"uniformGeometry"),pair(tab,"geometryMinusUniform"),atol=2e-6,rtol=0)
        np.testing.assert_allclose(pair(tab,"uniformGeometry").sum(axis=0),[0,0],atol=1e-5,rtol=0)
        if name=="layers":
            np.testing.assert_allclose(pair(tab,"geometryMinusUniform")[-2:].sum(axis=0),
                                       s["lastTwoLayersAfterUniformSubtraction"],atol=2e-6,rtol=0)
    far=np.genfromtxt(DATA/"far.csv",delimiter=",",names=True)
    np.testing.assert_array_equal(far["point"],np.arange(9701-109,9701))
    for tag,term in (("conv","farConvective"),("visc","farViscous")):
        np.testing.assert_allclose([fsum(far[tag+"MassKgS"]),fsum(far[tag+"EnergyW"])],terms[term],atol=1e-7,rtol=0)
    native=sum(terms[k] for k in ("storage","farConvective","farViscous","wall","symmetry"))
    np.testing.assert_allclose(native,s["nativeBoundaryPlusStorage"],atol=1e-7)
    geom=terms["internalRadiusMoment"]+terms["axisymmetricSource"]
    np.testing.assert_allclose(geom,s["radiusMomentPlusSource"],atol=1e-6)
    np.testing.assert_allclose(np.add(s["geometricAdvection"],s["geometricDiffusion"]),geom,atol=1e-6)
    np.testing.assert_allclose(native+geom,terms["physicalResidual"],atol=1e-6,rtol=0)
    saved=old["fieldBalance"]
    corrected=saved["totalEnergyImbalanceW"]-sum(d[1] for d in s["savedMinusNative"].values())
    np.testing.assert_allclose(corrected,native[1],rtol=0,atol=1e-7)
    href=saved["freeTotalEnthalpyJkg"]
    np.testing.assert_allclose(native[1]-href*native[0],s["nativeExcessImbalanceW"],atol=1e-7)
    # 用各对偶壁面的真实法向量核对通量符号和区域分配。
    wall=np.genfromtxt(DATA/"wall.csv",delimiter=",",names=True)
    np.testing.assert_array_equal(wall["point"],np.arange(109))
    tangent=np.diff(np.column_stack((wall["solverX"],wall["radiusM"])),axis=0)
    half_normal=np.column_stack((tangent[:,1],-tangent[:,0]))/2
    dual=np.vstack((half_normal[0],half_normal[:-1]+half_normal[1:],half_normal[-1]))
    amps=2*np.pi*wall["radiusM"]*wall["fluxWm2"]
    lengths=np.linalg.norm(dual,axis=1)
    np.testing.assert_allclose(amps*lengths,wall["weakPowerW"],rtol=2e-8,atol=1e-7)
    projected=amps[:,None]*dual/lengths[:,None]
    by_segment=np.einsum("ij,ij->i",half_normal,projected[:-1]+projected[1:])
    geom_model=read(DOC.parents[2]/"cases/simple_launcher/geometry.json")
    x_mid=geom_model["lengthM"]-.5*(wall["solverX"][:-1]+wall["solverX"][1:])
    ranges=read(DOC.parent/"v363/regions.json")["regions"]
    base=abs(x_mid)<1e-12
    for region in ranges:
        key=region["key"]
        select=base if key=="base" else ~base&(x_mid>region["xMinM"])&(x_mid<region["xMaxM"])
        np.testing.assert_allclose(fsum(by_segment[select]),s["regions"][key]["nativeWeakWallW"],atol=2e-6,rtol=0)
        assert s["regions"][key]["originalSurfaceW"]==surface["regions"][key]["researchIntoWallW"]
    np.testing.assert_allclose(fsum(by_segment),terms["wall"][1],atol=1e-5)
    print(f'109 boundary nodes, 89 layers, 109 columns and six wall regions verified; reference-adjusted balance {s["nativeExcessImbalanceW"]:.3f} W')


if __name__=="__main__":
    verify()
