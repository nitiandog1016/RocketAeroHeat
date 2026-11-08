"""严格读取原生真实壁面p/T/q；使用面编号，热流转为流入壁面为正。"""
from pathlib import Path
import re,json
import numpy as np
from prepareOpenFoamControl import sha
from openFoamFields import internal_field
from openFoamConservation import array_entry
from phaseAwareViscousAir import PhaseAwareViscousAir

FOOTER='// ************************************************************************* //'

def checked_header(text,name,kind,location=None,dimensions=None):
    """模块一：同时核对对象、单位和保存时刻，拒绝混帧或未写完文件。"""
    if not text.rstrip().endswith(FOOTER):raise ValueError('场或网格未完整写出')
    body=dictionary_body(text,'FoamFile')
    for key,value in [('format','ascii'),('object',name),('class',kind)]:
        if re.findall(r'\b'+key+r'\s+([^;]+);',body)!=[value]:raise ValueError('文件身份不符 '+key)
    if location is not None:
        values=re.findall(r'\blocation\s+"([^"]+)";',body)
        if values!=[location]:raise ValueError('保存时刻不符')
    if dimensions is not None:
        values=re.findall(r'\bdimensions\s+\[([^\]]+)\]\s*;',text)
        if len(values)!=1 or values[0].split()!=dimensions.split():raise ValueError('场单位错误')

def checked_contract(case):
    """模块二：本读取器只支持已测层流Fourier/给定壁温组合，不泛化湍流身份。"""
    gas=PhaseAwareViscousAir()
    contract=json.loads((case/'thermo_contract.json').read_text())
    if contract!=gas.viscous_contract():raise ValueError('不是已注册热物性合同')
    paths=[case/'constant'/n for n in ['physicalProperties','momentumTransport','thermophysicalTransport']]
    props,momentum,transport=[p.read_text() for p in paths]
    expected=dict(type='hePsiThermo',mixture='pureMixture',transport='icoTabulated',thermo='hIcoTabulated',
                  equationOfState='perfectGas',specie='specie',energy='rahPhaseAwareInternalEnergy')
    body=dictionary_body(props,'thermoType')
    if dict(re.findall(r'(\w+)\s+(\w+)\s*;',body))!=expected:raise ValueError('实际物性类型不符')
    compact=lambda s:re.sub(r'\s+','',s)
    mixture=props[props.index('mixture {'):]
    if compact(mixture)!=compact('mixture {'+gas.dictionary()+'}'):raise ValueError('实际物性表不符')
    if re.findall(r'\bsimulationType\s+(\w+);',momentum)!=['laminar']:raise ValueError('未验证的动量输运模型')
    if re.findall(r'\bsimulationType\s+(\w+);',transport)!=['laminar'] or not re.search(r'\bmodel\s+Fourier;',dictionary_body(transport,'laminar')):
        raise ValueError('未验证的热输运模型')
    return contract,{p.name:sha(p) for p in paths},gas

def dictionary_body(text,name):
    matches=re.findall(r'\b'+re.escape(name)+r'\s*\{([^{}]*)\}',text,re.S)
    if len(matches)!=1:raise ValueError('缺失或重复字典 '+name)
    return matches[0]

def read_wall_heat(case,time_name,cell_count):
    """模块三：只导出实际wall patch；该读取结果本身不授热资格。"""
    case=Path(case);folder=case/time_name
    if not isinstance(cell_count,(int,np.integer)) or isinstance(cell_count,bool) or cell_count<=0:raise ValueError('单元数错误')
    if not re.fullmatch(r'[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?',time_name) or not np.isfinite(float(time_name)) or float(time_name)<=0:raise ValueError('必须是实际非零保存时刻')
    contract,input_sha,gas=checked_contract(case)
    meshes=case/'constant/polyMesh';text=(meshes/'boundary').read_text()
    checked_header(text,'boundary','polyBoundaryMesh')
    patches=[]
    for name,body in re.findall(r'(\w+)\s*\{([^{}]*)\}',text,re.S):
        if re.search(r'\btype\s+wall\s*;',body):
            patches.append((name,int(re.search(r'\bnFaces\s+(\d+);',body).group(1)),int(re.search(r'\bstartFace\s+(\d+);',body).group(1))))
    if not patches:raise ValueError('没有物理wall patch')
    if len({p[0] for p in patches})!=len(patches):raise ValueError('重复wall名字')
    owner_text=(meshes/'owner').read_text();checked_header(owner_text,'owner','labelList')
    m=re.search(r'\n(\d+)\s*\(\s*([\s\d]+)\)',owner_text)
    if not m:raise ValueError('无法读取原生owner')
    owner=np.array([int(v) for v in m.group(2).split()])
    if len(owner)!=int(m.group(1)) or np.any(owner>=cell_count):raise ValueError('owner数量或范围错误')
    fields={n:(folder/n).read_text() for n in ['p','T','wallHeatFlux']}
    units={'p':'1 -1 -2 0 0 0 0','T':'0 0 0 1 0 0 0','wallHeatFlux':'1 0 -3 0 0 0 0'}
    for n,s in fields.items():
        checked_header(s,n,'volScalarField',time_name,units[n])
    pi=internal_field(folder/'p',cell_count);rows=[];seen=set()
    for name,count,start in patches:
        ids=np.arange(start,start+count)
        if count<=0 or start<0 or start+count>len(owner) or seen.intersection(ids.tolist()):raise ValueError('wall范围或重叠错误')
        seen.update(ids.tolist());bodies={n:dictionary_body(s[s.index('boundaryField'):],name) for n,s in fields.items()}
        if re.findall(r'\btype\s+(\w+);',bodies['T'])!=['fixedValue']:raise ValueError('不是给定壁温边界')
        p=array_entry(bodies['p'],'value',count);T=array_entry(bodies['T'],'value',count);raw=array_entry(bodies['wallHeatFlux'],'value',count)
        if p is None and re.search(r'type\s+zeroGradient;',bodies['p']):p=pi[owner[ids]]
        if p is None or T is None or raw is None:raise ValueError('没有真实壁面值 '+name)
        if np.any(p<=0) or np.any(T<=0):raise ValueError('壁面物理值非法')
        gas.evaluate(p,T)
        for j in range(count):rows.append(dict(patch=name,nativeFaceIndex=int(ids[j]),ownerCellIndex=int(owner[ids[j]]),
            pressurePa=float(p[j]),givenWallTemperatureK=float(T[j]),heatFluxIntoWallWm2=float(-raw[j]),
            nativeWallHeatFluxWm2=float(raw[j]),method='OpenFOAM13 shockFluid + Fourier + Lemmon icoTabulated',
            nativeSignConversion=-1,solidWallTemperatureK=None,edgeTemperatureK=None,physicalQualification=False))
    return dict(timeName=time_name,timeSeconds=float(time_name),rows=rows,contract=contract,
        fieldSHA256={n:sha(folder/n) for n in fields},meshSHA256={n:sha(meshes/n) for n in ['owner','boundary']},
        inputSHA256=input_sha,contractSHA256=sha(case/'thermo_contract.json'),
        originalRocketMapped=False,physicalQualification=False)
