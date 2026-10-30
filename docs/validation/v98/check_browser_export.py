"""模块用途：把真实浏览器下载的13行CSV和两个探针，独立对照冻结MAT。"""
from pathlib import Path
import csv
import json
import hashlib
import sys
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'tests'))
from auditPressureDrivenSurfaceFlow import rows

p=loadmat(ROOT/'tests/data/shock_layer_transport/parabolic_shock_layer_evidence.mat',simplify_cells=True)
t=loadmat(ROOT/'tests/data/shock_layer_transport/translated_conic_layer_evidence.mat',simplify_cells=True)
cases=[('browser_held_out_column24.csv',rows(t['heldOut']['columns'])[23]),
       ('browser_NASA7_column1.csv',rows(p['cases'])[1]['column'])]
checked=0
for name,col in cases:
 with (HERE/'raw'/name).open(encoding='utf-8-sig',newline='') as stream:records=list(csv.DictReader(stream))
 assert len(records)==len(col['eta'])
 for i,record in enumerate(records):
  for key,value in [('eta',col['eta'][i]),('normal_distance_m',col['normalDistanceM'][i]),
                    *[(axis+'_m',col['position'][i,j]) for j,axis in enumerate('xyz')],
                    ('layer_pressure_Pa',col['state']['pressure'][i]),('layer_temperature_K',col['state']['temperature'][i]),
                    ('layer_density_kg_m3',col['state']['density'][i]),('layer_speed_m_s',col['state']['speed'][i]),
                    ('entropy_J_kg_K',col['state']['entropy'][i]),('normal_jacobian_ratio',col['jacobianRatio'][i])]:
   assert record[key]!='' and float(record[key])==value,(name,i,key)
  for key in ('heat_flux_W_m2','wall_temperature_input_K','wall_temperature_solved_K','edge_temperature_K','recovery_temperature_K'):
   assert record[key]=='',key
  assert record['heat_flux_status']==record['wall_temperature_status']=='unknown'
  assert record['pressure_method']==col['state']['pressureMethod']
  assert record['temperature_method']==col['state']['temperatureMethod']
  checked+=1
receipt=json.loads((HERE/'raw/browser_final.json').read_text(encoding='utf8'))
for key,col,index in [('pickedPoint',cases[0][1],1),('profilePoint',cases[1][1],5)]:
 probe=receipt[key];values=dict(probe['values']);assert int(probe['station'])==index
 for label,expected in [('层内压力 / kPa',col['state']['pressure'][index]*.001),
                        ('层内温度 / K',col['state']['temperature'][index]),('层内速度 / m/s',col['state']['speed'][index]),
                        ('层内密度 / kg/m³',col['state']['density'][index]),('同流线来源熵 / J/(kg·K)',col['state']['entropy'][index])]:
  # 屏幕显示7位有效数字；CSV则按上方逐值完全相等核对，不借显示舍入放宽CSV。
  assert float(values[label])==float(format(float(expected),'.7g')),(key,label)
 assert values['热流 / W/m²']=='未求解'
 assert values['外缘Te / 恢复Tr / 给定Tw / 求解Tw']=='未求 / 未求 / 未设置 / 未求'
report=dict(status='passed-real-browser-CSV-and-rounded-probe-to-frozen-MAT-audit',CSVrows=checked,
            numericCSVValuesExactlyEqual=True,displaySignificantDigits=7,unknownValuesPreserved=True,
            downloads=[{'path':str((HERE/'raw'/name).relative_to(ROOT)).replace('\\','/'),
                        'sha256':hashlib.sha256((HERE/'raw'/name).read_bytes()).hexdigest()} for name,_ in cases],
            newProductionSolvedAreaM2=0)
target=HERE/'raw/browser_export_audit.json';assert not target.exists()
target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))
