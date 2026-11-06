"""模块用途。只读核对与v166共有完整场，补齐逐字比较证据，不启动求解器。"""
import json
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save


def main():
    base=ROOT/'outputs/v167_validation';out=base/'common_frame_comparison.json';assert not out.exists()
    assert read(base/'review01/report.json')['passedChecks']
    comparison=[]
    for step in (0,100):
        for pattern in ('restart_flow_%06d.csv','surface_flow_%06d.vtk'):
            name=pattern%step;old=ROOT/'outputs/v166_validation/real02'/name;new=base/'real01'/name
            equal=sha(old)==sha(new);assert equal
            comparison.append(dict(localIteration=step,baselinePath=str(old),candidatePath=str(new),
                baselineSHA256=sha(old),candidateSHA256=sha(new),allBytesEqual=equal))
    save(out,dict(passed=True,pairs=comparison,nativeRuns=0,physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(dict(passed=True,completeFieldAndSurfacePairs=len(comparison))))


if __name__=='__main__':main()
