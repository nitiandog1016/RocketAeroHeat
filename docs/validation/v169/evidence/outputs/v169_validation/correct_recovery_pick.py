"""目视叠图后修正与测量圆圈重叠的一点；准确旧表保留，不重跑原生任务。"""
import csv,json,sys,shutil
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from prepareNACA3344Benchmark import digitize
from su2ThermalRuntime import sha,read
from validateIdealAir import save
data=ROOT/'tests/data/naca_tn3344/normalization';out=ROOT/'outputs/v169_validation/recovery02'
out.mkdir(exist_ok=False)
for name in ('manifest.json','digitization.json','recovery_factor.csv'):
    shutil.copyfile(str(data/name),str(out/('before_'+name)))
spec=read(data/'digitization.json');assert spec['pixels'][7]==[1448,629.5]
spec['pixels'][7]=[1430,635]
values,error=digitize(spec)
save(data/'digitization.json',spec)
with (data/'recovery_factor.csv').open('w',newline='',encoding='utf8') as f:
    writer=csv.writer(f);writer.writerow(['s_over_D','Cr','reading_uncertainty_Cr','source_pixel_x','source_pixel_y'])
    writer.writerows([[*v,.0005,*px] for v,px in zip(values,spec['pixels'])])
manifest=read(data/'manifest.json')
manifest['files']={name:sha(data/name) for name in manifest['files']}
save(data/'manifest.json',manifest)
save(out/'report.json',dict(passed=True,beforePixel=[1448,629.5],afterPixel=[1430,635],
    reason='point overlapped measured circle; use adjacent unambiguous solid line',
    independentGridError=error.tolist(),manifestSHA256=sha(data/'manifest.json'),CFDReruns=0))
