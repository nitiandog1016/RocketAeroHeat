"""严格读取原生真实壁面p/T/q；使用面编号，热流转为流入壁面为正。"""
from pathlib import Path
import re,json
import numpy as np
from prepareOpenFoamControl import sha
from openFoamFields import internal_field
from openFoamConservation import array_entry

def dictionary_body(text,name):
    matches=re.findall(r'\b'+re.escape(name)+r'\s*\{([^{}]*)\}',text,re.S)
    if len(matches)!=1:raise ValueError('缺失或重复字典 '+name)
    return matches[0]

def read_wall_heat(case,time_name,cell_count):
    """模块一：只导出实际wall patch；该读取结果本身不授热资格。"""
    case=Path(case);folder=case/time_name
    if not np.isfinite(float(time_name)) or float(time_name)<=0:raise ValueError('必须是实际非零保存时刻')
    contract=json.loads((case/'thermo_contract.json').read_text())
    if contract.get('transport')!='rah-lemmon-dilute-icoTabulated-of13-v1' or contract.get('caloric',{}).get('model')!='rah-phase-aware-dry-air-of13-v1':raise ValueError('不是已注册热物性组合')
    meshes=case/'constant/polyMesh';text=(meshes/'boundary').read_text()
    patches=[]
    for name,body in re.findall(r'(\w+)\s*\{([^{}]*)\}',text,re.S):
        if re.search(r'\btype\s+wall\s*;',body):
            patches.append((name,int(re.search(r'\bnFaces\s+(\d+);',body).group(1)),int(re.search(r'\bstartFace\s+(\d+);',body).group(1))))
    if not patches:raise ValueError('没有物理wall patch')
    owner_text=(meshes/'owner').read_text();m=re.search(r'\n(\d+)\s*\(\s*([\s\d]+)\)',owner_text)
    if not m:raise ValueError('无法读取原生owner')
    owner=np.array([int(v) for v in m.group(2).split()])
    if len(owner)!=int(m.group(1)) or np.any(owner>=cell_count):raise ValueError('owner数量或范围错误')
    fields={n:(folder/n).read_text() for n in ['p','T','wallHeatFlux']}
    for n,s in fields.items():
        if not re.search(r'format\s+ascii;',s) or not s.rstrip().endswith('// ************************************************************************* //'):raise ValueError('场未完整写出')
        if not re.search(r'object\s+'+n+r';',s):raise ValueError('场身份不符')
        if n=='wallHeatFlux' and not re.search(r'dimensions\s+\[1\s+0\s+-3\s+0\s+0\s+0\s+0\]',s):raise ValueError('热流单位错误')
    pi=internal_field(folder/'p',cell_count);rows=[];seen=set()
    for name,count,start in patches:
        ids=np.arange(start,start+count)
        if start<0 or start+count>len(owner) or seen.intersection(ids.tolist()):raise ValueError('wall范围或重叠错误')
        seen.update(ids.tolist());bodies={n:dictionary_body(s[s.index('boundaryField'):],name) for n,s in fields.items()}
        p=array_entry(bodies['p'],'value',count);T=array_entry(bodies['T'],'value',count);raw=array_entry(bodies['wallHeatFlux'],'value',count)
        if p is None and re.search(r'type\s+zeroGradient;',bodies['p']):p=pi[owner[ids]]
        if p is None or T is None or raw is None:raise ValueError('没有真实壁面值 '+name)
        if np.any(p<=0) or np.any(T<=0):raise ValueError('壁面物理值非法')
        for j in range(count):rows.append(dict(patch=name,nativeFaceIndex=int(ids[j]),ownerCellIndex=int(owner[ids[j]]),
            pressurePa=float(p[j]),givenWallTemperatureK=float(T[j]),heatFluxIntoWallWm2=float(-raw[j]),
            nativeWallHeatFluxWm2=float(raw[j]),method='OpenFOAM13 shockFluid + Fourier + Lemmon icoTabulated',
            nativeSignConversion=-1,solidWallTemperatureK=None,edgeTemperatureK=None,physicalQualification=False))
    return dict(timeName=time_name,timeSeconds=float(time_name),rows=rows,contract=contract,
        fieldSHA256={n:sha(folder/n) for n in fields},meshSHA256={n:sha(meshes/n) for n in ['owner','boundary']},
        originalRocketMapped=False,physicalQualification=False)
