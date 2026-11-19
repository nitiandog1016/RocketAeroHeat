"""以 E 盘保存的同工况原面与两个独立导出核对三维对照的值、身份和显示范围。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from auditSU2SSTState import vtk_tokens
from cm011Tunnel import CP, D, DT, RHO, U, meridian, write_json
from runEmbeddedSU2 import sha


def page(path):
    text=path.read_text(encoding="utf8").strip()
    prefix="window.ROCKET_CASE="
    assert text.startswith(prefix) and text.endswith(";")
    return json.loads(text[len(prefix):-1])


def check(folder,sector):
    old=page(ROOT/"outputs/v272_validation/workbench01/workbench/case-data.js")
    published=folder/"workbench/case-data.js"
    new=page(published)
    assert old["faces"]==new["faces"] and old["vertices"]==new["vertices"]
    assert old["wallCoverage"]==new["wallCoverage"] and old["nativeEnergyBalance"]==new["nativeEnergyBalance"]
    assert new["researchIdentity"]["productionQualified"] is False
    info=json.loads((folder/"sector_audit/summary.json").read_bytes())
    exposed=new["cm011SectorSummary"]
    assert [s["name"] for s in exposed["byRegion"]]==[s["name"] for s in info["byRegion"]]
    for k in ("savedIndex","lastHistoryIndex","normalTerminal","timedOut","numericalWallAccepted","bodySHA256"):
        assert exposed[k]==info[k],k
    if "nativeEnergyBalance" in info:
        receipt=info["nativeEnergyBalance"]
        posted=exposed["nativeEnergyBalance"]
        for key in ("sourceRestartSHA256","wallMinusFarFraction","nativeWallReleaseW","nativeFarExcessOutW","belowHalfPercent"):
            assert posted[key]==receipt[key]
    rows=np.genfromtxt(folder/"sector_audit/wall.csv",delimiter=",",names=True)
    fields={f["key"]:f for f in new["fields"]}
    assert all(fields[f["key"]]==f for f in old["fields"])
    matrix=(("sectorGaugePressure","gaugePressurePa",1.),("sectorHeatFlux","heatFluxWm2",.001),
            ("sectorWallNode","wallNodeK",1.),("sectorFirstNode","firstNodeK",1.))
    assert len(fields)==len(old["fields"])+len(matrix)
    assert np.array_equal(rows["faceID"],np.arange(len(new["faces"])))
    # 独立读原生表面帧。按坐标找到母线，直接逐顶点调用一维角度插值以核对折叠和原编号。
    surf,native=vtk_tokens(sector/f'surface_{info["savedIndex"]:06d}.vtk')
    rays=np.deg2rad([0.,3.,6.])
    line,_=meridian(2)
    coordinates=np.array(new["vertices"])/D
    distances,meridian_ids=cKDTree(line).query(np.column_stack((coordinates[:,0],np.hypot(coordinates[:,1],coordinates[:,2]))))
    assert max(distances)<1e-11
    az=np.arctan2(coordinates[:,2],coordinates[:,1])
    folded=abs(az-2*rays[-1]*np.rint(az/(2*rays[-1])))
    tree=cKDTree(surf)
    source=[]
    for phi in rays:
        distances,ids=tree.query(np.column_stack((line[:,0],line[:,1]*np.cos(phi),line[:,1]*np.sin(phi))))
        assert max(distances)<2e-6
        source.append(ids)
    for variable,column,scale in (("Heat_Flux","heatFluxWm2",RHO*CP*U*DT),("Pressure","gaugePressurePa",RHO*U**2)):
        values=native[variable][:,0][np.array(source)]*scale
        points=np.array([np.interp(a,rays,values[:,j]) for a,j in zip(folded,meridian_ids)])
        field=points[np.array(new["faces"])].mean(axis=1)
        np.testing.assert_allclose(rows[column],field,rtol=1e-11,atol=1e-13)
    for key,column,scale in matrix:
        f=fields[key]
        np.testing.assert_allclose(f["values"],rows[column]*scale,rtol=1e-12,atol=1e-15)
        assert all(s==1 for s in f["status"])
        assert all("三维" in f["methodNames"][i] for i in f["methodIndex"])
        assert str(info["savedIndex"]) in f["stateNote"]
        assert ("未达稳态验收" in f["computedStateLabel"]) == (not info["numericalWallAccepted"])
        if "nativeEnergyBalance" in info:
            assert f'{100*receipt["wallMinusFarFraction"]:.4f}%' in f["stateNote"]
    for k1,k2 in (("heatFlux","sectorHeatFlux"),("pressureGauge","sectorGaugePressure")):
        limits=new["comparisonColorRanges"][k1]
        assert limits==new["comparisonColorRanges"][k2]
        vals=fields[k1]["values"]+fields[k2]["values"]
        np.testing.assert_allclose(limits,[min(vals),max(vals)],rtol=1e-13,atol=1e-15)
    area=np.array(new["faceArea"])
    power=np.dot(area,1000*np.array(fields["sectorHeatFlux"]["values"]))
    np.testing.assert_allclose(power,info["fullTrianglePowerW"],rtol=1e-12,atol=1e-14)
    return dict(faces=len(area),referenceFields=len(old["fields"]),newFields=[e[0] for e in matrix],
                savedIndex=info["savedIndex"],originalFieldDataEqual=True,independentNativeFaceMapping=True,coverageEqual=True,
                sectorWallAccepted=info["numericalWallAccepted"],independentPowerW=power,
                webPayloadSHA256=sha(published),sectorRowsSHA256=sha(folder/"sector_audit/wall.csv"),
                browserInteractionTested=False)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbench_run",type=Path)
    parser.add_argument("sector_run",type=Path)
    parser.add_argument("new_report",type=Path)
    args=parser.parse_args()
    write_json(args.new_report,check(args.workbench_run,args.sector_run))
