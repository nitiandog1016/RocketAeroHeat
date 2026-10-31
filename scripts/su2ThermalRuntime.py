"""原生空间壁温运行时选择和实际施加值审核；未知/不匹配输入不能降级。"""
from pathlib import Path
import csv,hashlib,json,math,re,bisect

ROOT=Path(__file__).resolve().parents[1]
FLAVOR='rah-axial-wall-v1'
FROZEN_FLAVOR='rah-frozen-energy-v1'
COOLPROP_FLAVOR='rah-coolprop-v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))

def resolve(folder,config):
    """只选择固定官方包或明确受控构建；绑定清单、温度表和运行时身份。"""
    manifest=folder/'case_manifest.json';locked={};profile=None
    data=read(manifest) if manifest.exists() else {}
    flavor=data.get('requiredRuntimeFlavor','official')
    assert flavor in ('official',FLAVOR,FROZEN_FLAVOR,COOLPROP_FLAVOR),'未知SU2运行时类型'
    if flavor!=FROZEN_FLAVOR:assert 'frozenEnergyMode' not in data,'运行时类型与冻结能量选项不匹配'
    runtime_dir=ROOT/'tools/runtime/su2'/('8.6.0' if flavor=='official' else '8.6.0-'+flavor)
    runtime_file=runtime_dir/'runtime.json'
    assert runtime_file.is_file(),'空间壁温受控运行时尚未构建，不能退回官方标量壁温'
    runtime=read(runtime_file)
    if flavor==FLAVOR:
        assert runtime['flavor']==FLAVOR and data['thermalBoundary']['mode']=='axialTemperature'
        assert re.search(r'^\s*MARKER_PYTHON_CUSTOM\s*=\s*\(\s*wall\s*\)',config,re.M)
        assert re.search(r'^\s*MARKER_ISOTHERMAL\s*=\s*\(\s*wall\s*,',config,re.M)
        assert not re.search(r'^\s*MARKER_HEATFLUX\s*=',config,re.M)
        assert re.search(r'^\s*MGLEVEL\s*=\s*0\s*$',config,re.M)
        assert re.search(r'^\s*SYSTEM_MEASUREMENTS\s*=\s*SI\s*$',config,re.M)
        assert not re.search(r'^\s*TIME_DOMAIN\s*=\s*YES',config,re.M)
        profile=folder/'wall_temperature.csv';assert profile.is_file()
        assert sha(profile)==data['thermalProfileSHA256'],'温度表与输入清单不符'
        locked.update({str(manifest):sha(manifest),str(profile):sha(profile)})
    else:
        assert not (folder/'wall_temperature.csv').exists(),'有温度表却未声明受控运行时'
        if flavor==FROZEN_FLAVOR:
            assert runtime['flavor']==FROZEN_FLAVOR
            validate_frozen_manifest(data,config)
            locked[str(manifest)]=sha(manifest)
        elif flavor==COOLPROP_FLAVOR:
            assert runtime['flavor']==COOLPROP_FLAVOR
            validate_coolprop_manifest(data,config)
            locked[str(manifest)]=sha(manifest)
        else:
            assert 'frozenEnergyMode' not in data,'官方运行时不能默默忽略冻结能量选项'
    exe=Path(runtime['executable']).resolve();exe.relative_to(runtime_dir.resolve())
    for item in runtime['files']:
        p=(runtime_dir/item['path']).resolve();p.relative_to(runtime_dir.resolve())
        assert sha(p)==item['sha256'],'运行时文件身份变化'
    locked[str(runtime_file)]=sha(runtime_file)
    return exe,runtime,flavor,profile,locked

def validate_coolprop_manifest(data,config):
    """首阶段只接受已准备的固定Air/单区层流物性任务；旧rhoE或湍流列不能混用。"""
    settings={}
    for line in config.splitlines():
        line=line.split('%',1)[0].strip()
        if '=' in line:
            key,value=map(str.strip,line.split('=',1));assert key not in settings;settings[key]=value
    required={'SOLVER':'NAVIER_STOKES','KIND_TURB_MODEL':'NONE','FLUID_MODEL':'COOLPROP','FLUID_NAME':'Air',
              'SYSTEM_MEASUREMENTS':'SI','REF_DIMENSIONALIZATION':'DIMENSIONAL','MGLEVEL':'0','RESTART_SOL':'NO',
              'VISCOSITY_MODEL':'SUTHERLAND','CONDUCTIVITY_MODEL':'CONSTANT_PRANDTL'}
    for k,v in required.items():assert settings.get(k)==v,'当前CoolProp适配不接受 '+k
    assert settings.get('CONV_NUM_METHOD_FLOW') in ('ROE','HLLC')
    assert settings.get('TIME_DOMAIN','NO')=='NO' and 'MARKER_RIEMANN' in settings
    for k in ('MARKER_FAR','MARKER_INLET','MARKER_OUTLET','CONFIG_LIST','NZONES','MARKER_PYTHON_CUSTOM'):
        assert k not in settings,'需要已验证的一般流体边界：'+k
    assert ('MARKER_ISOTHERMAL' in settings)!=('MARKER_HEATFLUX' in settings)
    assert data.get('thermodynamicModel')=='CoolProp-7.2.0-HEOS-Air'
    assert float(settings['PRANDTL_LAM'])==.72,'当前审核与显示仅支持显式Pr0.72'
    return settings

def validate_frozen_manifest(data,config):
    """实验恢复只接受显式选择的固定单区稳态NS；绝不从环境继承或默默降级。"""
    assert data.get('frozenEnergyMode') in ('diagnostic','apply'),'须明确选择诊断或恢复'
    assert data.get('thermalBoundary',{}).get('mode') in ('adiabatic','isothermal')
    required={'SOLVER':'NAVIER_STOKES','KIND_TURB_MODEL':'NONE','MATH_PROBLEM':'DIRECT',
              'SYSTEM_MEASUREMENTS':'SI','REF_DIMENSIONALIZATION':'DIMENSIONAL','FLUID_MODEL':'IDEAL_GAS',
              'MGLEVEL':'0','TIME_DISCRE_FLOW':'EULER_IMPLICIT'}
    settings={}
    for line in config.splitlines():
        line=line.split('%',1)[0].strip()
        if '=' in line:
            key,value=map(str.strip,line.split('=',1));assert key not in settings,'重复配置键';settings[key]=value
    for key,value in required.items():assert settings.get(key)==value,'实验恢复不支持配置 '+key
    for key in ('CONFIG_LIST','MARKER_PYTHON_CUSTOM','NZONES'):assert key not in settings
    assert settings.get('TIME_DOMAIN','NO')=='NO'
    assert ('MARKER_HEATFLUX' in settings)!=('MARKER_ISOTHERMAL' in settings)
    if data['thermalBoundary']['mode']=='adiabatic':
        assert re.fullmatch(r'\(\s*wall\s*,\s*0(?:\.0*)?\s*\)',settings.get('MARKER_HEATFLUX',''))
    else:
        match=re.fullmatch(r'\(\s*wall\s*,\s*([\d.eE+-]+)\s*\)',settings.get('MARKER_ISOTHERMAL',''))
        assert match and math.isfinite(float(match.group(1))) and float(match.group(1))>0
    return data['frozenEnergyMode']

def audit_applied(folder,profile,mesh):
    """从原网格节点/墙面集合和表重新计算，审核C++实际写入的逐节点值。"""
    with profile.open(encoding='utf-8') as f:
        assert f.readline().strip()=='RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K'
        pairs=[list(map(float,row)) for row in csv.reader(f)]
    assert len(pairs)>=2 and all(len(row)==2 for row in pairs)
    xx=[q[0] for q in pairs];tt=[q[1] for q in pairs]
    assert all(math.isfinite(x) for q in pairs for x in q) and min(tt)>0
    assert all(b>a for a,b in zip(xx,xx[1:]))
    def expected(x):
        assert xx[0]<=x<=xx[-1]
        j=bisect.bisect_left(xx,x)
        if j==0:return tt[0]
        w=(x-xx[j-1])/(xx[j]-xx[j-1]);return (1-w)*tt[j-1]+w*tt[j]
    # SU2全局节点编号来自本程序显式写出的最后一列，不能按重新排序后的行号猜测。
    nodes={};wall=set();lines=mesh.read_text(encoding='utf-8').splitlines();i=0;ndim=None
    while i<len(lines):
        line=lines[i].strip();i+=1
        if line.startswith('NDIME='):ndim=int(line.split('=')[1])
        elif line.startswith('NPOIN='):
            count=int(line.split('=')[1])
            for row in lines[i:i+count]:
                v=row.split();nodes[int(v[ndim])]=tuple(map(float,v[:ndim]))
            i+=count
        elif line.startswith('MARKER_TAG=') and line.split('=')[1].strip()=='wall':
            count=int(lines[i].split('=')[1]);i+=1
            for row in lines[i:i+count]:wall.update(map(int,row.split()[1:]))
            i+=count
    assert ndim in (2,3) and wall
    file=folder/'rah_wall_temperature_applied.csv'
    with file.open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
    ids=[int(r['globalNodeID']) for r in rows];assert len(set(ids))==len(rows) and set(ids)==wall
    errors=[]
    for row,j in zip(rows,ids):
        xyz=tuple(float(row[n]) for n in ('X_m','Y_m','Z_m'))
        assert all(abs(a-b)<=1e-12*max(1,abs(b)) for a,b in zip(xyz,nodes[j]))
        value=float(row['givenTemperature_K']);want=expected(nodes[j][0]);errors.append(abs(value-want))
        assert math.isfinite(value) and abs(value-want)<1e-11*max(1,abs(want))
    return dict(passed=True,wallNodes=len(rows),maximumTemperatureDifferenceK=max(errors),
                profileX_m=xx,profileTemperature_K=tt,
                appliedFileSHA256=sha(file),profileSHA256=sha(profile),
                temperatureRangeK=[min(float(r['givenTemperature_K']) for r in rows),
                                   max(float(r['givenTemperature_K']) for r in rows)],
                physicsQualified=False)
