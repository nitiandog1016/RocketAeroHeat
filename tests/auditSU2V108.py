"""独立读取v108原面结果，核对原始VTK与两种热负荷；不运行求解器。"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from auditEmbeddedSU2FlatPlate import read_wall

ROOT=Path(__file__).resolve().parents[1]


def main():
    folder=ROOT/'outputs/v108_validation';native=folder/'normal01/native_input'
    paths=[folder/'normal01/normal_volume.mat',folder/'import01/mapped_startup.mat',native/'surface_flow.vtk',
           native/'engine_process.json',native/'history.csv',native/'engine_stdout.log',Path(__file__).resolve()]
    saved=loadmat(str(paths[0]),squeeze_me=True,struct_as_record=False)
    result=loadmat(str(paths[1]),squeeze_me=True,struct_as_record=False)['result'];p=saved['projection'];v=saved['volume']
    xyz,wall=read_wall(paths[2]);ids=np.asarray(p.wallVertexIDs,dtype=int)-1
    distance,index=cKDTree(v.vertices[ids]).query(xyz);assert len(np.unique(index))==len(ids) and max(distance)<1e-6
    fields=np.zeros((len(v.vertices),3));fields[ids[index]]=np.column_stack([wall[k] for k in ('Pressure','Heat_Flux','Temperature')])
    expected=p.weights@fields;actual=np.column_stack((result.pressureRaw,result.heatFluxRaw,result.fluidNodeTemperatureMeanRaw))
    error=float(np.max(abs(expected-actual)/np.maximum(1,abs(expected))));assert error<1e-12
    area=np.linalg.norm(np.cross(result.V[result.F[:,1]-1]-result.V[result.F[:,0]-1],
                                result.V[result.F[:,2]-1]-result.V[result.F[:,0]-1]),axis=1)/2
    assert np.allclose(area,result.faceArea,rtol=2e-12,atol=1e-15)
    assert np.array_equal(result.V,saved['surface'].V) and np.array_equal(result.F,saved['surface'].F)
    physical=float(area@result.heatFluxRaw);load=float(area@result.conservativeHeatLoadDensityRaw)
    history=np.genfromtxt(str(paths[4]),delimiter=',',names=True);native_hf=float(history['HF'][-1])
    assert abs(load/native_hf-1)<2e-6 and not result.productionQualified and not result.eligibleForSurfaceCoverage
    log=paths[5].read_text(encoding='utf8')
    orth=float(re.search(r'Orthogonality Angle \(deg.\)\|\s*([\d.eE+-]+)',log)[1])
    cv=float(re.search(r'CV Volume\|\s*([\d.eE+-]+)',log)[1]);assert cv>0
    report=dict(originalFaces=len(area),mappingMaximumNormalizedError=error,originalArea=float(area.sum()),
                physicalAreaHeatPower=physical,nativeMappedHeatLoadPower=load,engineHeatPower=native_hf,
                engineHeatPowerRelativeError=abs(load/native_hf-1),
                fullDomainNativeOrthogonalityMinDeg=orth,nativeDualControlVolumeMinimum=cv,
                fullDomainMetricIsNotNearWallAlignment=True,
                lastEnergyResidual=float(history['rmsRhoE'][-1]),physicsQualified=False,newQualifiedArea=0,
                inputSHA256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    dest=ROOT/'docs/validation/v108/independent_mapping_audit.json';assert not dest.exists()
    dest.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'))
    print(json.dumps({k:v for k,v in report.items() if k!='inputSHA256'},indent=2))


if __name__=='__main__':main()
