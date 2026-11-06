"""模块用途：保存已目视核对的原文Cr归约曲线，理论辅助量不当作热流测量。"""
from pathlib import Path
import csv
import shutil

import numpy as np

from naca3344Model import DATA
from prepareNACA3344Benchmark import digitize
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save


def main():
    out=ROOT/'outputs/v169_validation/recovery01';out.mkdir(parents=True,exist_ok=False)
    (out/'executed_prepare.py').write_bytes(Path(__file__).read_bytes())
    source=ROOT/'outputs/v168_validation/sources/normalization-47.png'
    pdf=ROOT/'outputs/v168_validation/sources/NACA_TN_3344.pdf'
    assert sha(pdf)==read(DATA/'source_manifest.json')['sha256']
    # 三点标定和第四点独立核验；实线是作者对实验的恢复因子归约，不是热流实测。
    spec=dict(pixelCalibration=[[414,641,0,.82],[2339,653,1.3,.82],[414,494,0,.84]],
        independentAxisCheck=[2338,504,1.3,.84],
        pixels=[[414,597],[562,599.5],[710,600.5],[858,603.5],[1006,607],[1154,612.5],
            [1302,623.5],[1448,629.5],[1523,587.5],[1597,515],[1672,487.5],
            [1746,478],[1894,474],[2042,476]])
    values,error=digitize(spec)
    assert max(error)<.001 and np.all((values[:,1]>.8)&(values[:,1]<.86))
    folder=DATA/'normalization';folder.mkdir(exist_ok=False)
    shutil.copyfile(str(source),str(folder/'figure3b.png'))
    save(folder/'digitization.json',dict(spec,source='NACA TN3344 PDF page47/printed46 Figure3(b), solid Present method',
        coefficientIsDirectHeatMeasurement=False,readingUncertaintyCr=.0005,readingUncertaintyS_over_D=.004))
    with (folder/'recovery_factor.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.writer(f);writer.writerow(['s_over_D','Cr','reading_uncertainty_Cr','source_pixel_x','source_pixel_y'])
        writer.writerows([[*v,.0005,*px] for v,px in zip(values,spec['pixels'])])
    manifest=dict(sourcePDF_SHA256=sha(pdf),sourcePDFPage=47,printedPage=46,
        sourceURL='https://ntrs.nasa.gov/citations/19930084172',
        meaning='authors equation12 recovery-factor curve used in equation10 reduction; not measured heat flux',
        equation9='Cr0=sqrt(Pr((T0+Te)/2)); source sum-temperature bound 2000 degree Rankine',
        equation12='Cr=Cr0*(1-0.023*mbar); transformed mbar is not raw local Mach gradient',
        chosenImplementation='interpolate published equation12 curve; do not rederive mbar from a raw gradient',
        localEdge='must be supplied by qualified local solution; freestream is not interchangeable',
        conductivityDifference='solver uses Lemmon/Sutherland/Pr0.72; paper used reference26 tables; equivalence not established',
        actualExperimentalSegmentBounds=None,
        unresolvedSegmentMapping='Figure1c interlaced meridians, 12deg nose electrode intervals and 27/64 inch cylinder intervals; individual heat points not yet assigned to exact bounds',
        exactExperimentalComparisonQualified=False,
        files={p.name:sha(p) for p in folder.iterdir() if p.is_file()})
    save(folder/'manifest.json',manifest)
    save(out/'report.json',dict(passed=True,points=len(values),independentGridError=error.tolist(),
        sourceSHA256={str(source):sha(source),str(pdf):sha(pdf)},manifestSHA256=sha(folder/'manifest.json'),
        experimentalQualification=False,newQualifiedArea_m2=0))
    print(dict(points=len(values),independentGridError=error.tolist()))


if __name__=='__main__':main()
