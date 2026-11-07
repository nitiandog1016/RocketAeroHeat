"""读取原生ASCII保存场并审计物性；不用显示色标推断求解资格。"""
from pathlib import Path
import re,json
import numpy as np
from openFoamAirThermo import evaluate


def internal_field(path,count,vector=False):
    """模块一：严格读取均匀/非均匀内部场，拒绝数目不符、非有限及尾部缺失。"""
    text=Path(path).read_text(encoding='utf-8')
    if not re.search(r'format\s+ascii\s*;',text):raise ValueError('只接受已登记ASCII')
    m=re.search(r'internalField\s+uniform\s+([^;]+);',text)
    columns=3 if vector else 1
    if m:
        values=np.fromstring(m.group(1).replace('(',' ').replace(')',' '),sep=' ')
        if values.size!=columns:raise ValueError('均匀场维度错误')
        a=np.tile(values,(count,1))
    else:
        m=re.search(r'internalField\s+nonuniform\s+List<(scalar|vector)>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
        if not m or int(m.group(2))!=count or m.group(1)!=('vector' if vector else 'scalar'):
            raise ValueError('保存场不完整或数量不匹配')
        values=np.fromstring(m.group(3).replace('(',' ').replace(')',' '),sep=' ')
        if values.size!=count*columns:raise ValueError('保存场数字数量错误')
        a=values.reshape(count,columns)
    if not np.isfinite(a).all():raise ValueError('保存场包含非有限值')
    return a if vector else a[:,0]


def audit_state(fields):
    """模块二：热内能允许负参考值；正性靠p/rho/T与同一NASA7物性核对。"""
    p,T,rho,e=[fields[k] for k in ['p','T','rho','e']]
    if not all(np.isfinite(fields[k]).all() for k in fields):raise ValueError('非有限流场')
    if np.any(rho<=0):raise ValueError('非正密度')
    gas=evaluate(p,T)
    erho=float(np.max(abs(rho/gas['densityKgm3']-1)))
    # 原生温度逆算用1e-4相对种子容差；事前严格2e-8能量缩放门槛独立于显示精度。
    ee=float(np.max(abs(e-gas['internalEnergyJkg'])/np.maximum(1,np.abs(gas['internalEnergyJkg']))))
    if erho>2e-10 or ee>2e-8:raise ValueError('保存场EOS/内能不相容')
    return dict(TminK=float(T.min()),TmaxK=float(T.max()),rhoMinKgm3=float(rho.min()),pMinPa=float(p.min()),
        densityEOSMaxRelative=erho,internalEnergyMaxScaled=ee)


def read_state(directory,count):
    return {name:internal_field(Path(directory)/name,count,name=='U') for name in ['p','T','rho','U','e']}


class SavedGuard:
    """模块三：只读取完成的保存目录；日志警告导致拒绝，不声称瞬时查询前保护。"""
    def __init__(self,case,count,report):
        self.case=Path(case);self.count=count;self.report=Path(report);self.checked={}
    def __call__(self,log,final):
        text=Path(log).read_text(errors='replace')
        for token in ['out of temperature range','Negative initial temperature','Maximum number of iterations exceeded','FOAM FATAL']:
            if token in text:return dict(reason='native-warning-or-fatal',token=token)
        for step in range(100,601,100):
            if str(step) in self.checked:continue
            d=self.case/str(step)
            if not d.is_dir():continue
            files=[d/n for n in ['p','T','rho','U','e']]
            complete=all(p.is_file() and p.read_bytes().rstrip().endswith(b'// ************************************************************************* //') for p in files)
            if not complete and not final:continue
            try:self.checked[str(step)]=audit_state(read_state(d,self.count))
            except (ValueError,OSError) as error:return dict(reason='saved-state-rejected',step=step,error=str(error))
            self.report.write_text(json.dumps(dict(frames=self.checked,physicalQualification=False),indent=2)+'\n')
        return None
