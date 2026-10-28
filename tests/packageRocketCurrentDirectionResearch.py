# -*- coding: utf-8 -*-
"""冻结当前方向研究源码、完整替换文件及离线预览；不声称隔离MATLAB工程重算。"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
target=VAL/'current_direction_research_package_verification.json';assert not target.exists()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()


manifest=ROOT/'docs/validation/v83/progress_current_state_direction.json'
m=json.loads(manifest.read_text(encoding='utf8'))
assert m['production']['goalStatus']=='active' and not m['production']['version83IsolatedCompleteEngineeringDelivery']
names=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT).split(b'\0')
source=sorted({ROOT/Path(n.decode('utf8')) for n in names if n and (ROOT/Path(n.decode('utf8'))).is_file()})
replace=[ROOT/p for p in m['currentSourcesSHA256']]+[manifest,Path(__file__)]
preview=VAL/'workbench35'
readme='''RocketAeroHeat v83 当前真实状态方向研究源码（2026-10-26）

性质：研究源码快照/完整替换文件；本包没有隔离MATLAB工程重算，不是正式v83工程交付。
最新正式冻结交付仍为v82。正式合格热流覆盖66.5626327424%，本轮新增0。

源码快照保留当前非忽略源码、测试、模型和文档。main.m为工程入口。
每个替换文件都是完整文件，沿工程相对路径替换；替换包需要原工程及相应依赖。
保留MATLAB后端、X主轴、鼻尖+X、零攻角来流-X和SI单位。

research-preview/index.html为真实krylov22的离线研究工作台预览。
新储库问题该分支累计2接受/9次主NS；旧特征问题13接受/258次不拼入新历史。
最新一次RMS略降、最大局部残差升高，尚未稳态，所有场量保留算法和研究身份。
Tw=300K是输入，Te/Tr/求解Tw未知。科研PNG和原MAT/Node数据检查不是浏览器/WebGL视觉验收。
Q21与正则化方向Q22都是Q20的子分支，Q22结果较好但步长仍小；不合并两分支计数。

数值研究重算依赖本机outputs中的大型MAT检查点及公开参考缓存，本源码包不包含这些依赖。
ZIP CRC和逐文件SHA通过不表示独立重算或区域物理验证完成。
详见docs/CURRENT_STATE_DIRECTION_V83.md与docs/validation/v83/progress_current_state_direction.json。
'''
packages=[]
for filename,files,prefix in [('RocketAeroHeat_v83_current_state_research_source.zip',source,'RocketAeroHeat/'),
                              ('RocketAeroHeat_v83_current_state_complete_replacements.zip',sorted(set(replace)),'')]:
    path=VAL/filename;assert not path.exists();entries={}
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for item in files:
            item.resolve().relative_to(ROOT.resolve());name=prefix+item.relative_to(ROOT).as_posix()
            assert name not in entries
            data=item.read_bytes();z.writestr(name,data);entries[name]=hashlib.sha256(data).hexdigest()
        for item in sorted(p for p in preview.rglob('*') if p.is_file()):
            name=prefix+'research-preview/'+item.relative_to(preview).as_posix();assert name not in entries
            data=item.read_bytes();z.writestr(name,data);entries[name]=hashlib.sha256(data).hexdigest()
        name=prefix+'RESEARCH_PACKAGE_README.txt';data=readme.encode('utf8');z.writestr(name,data);entries[name]=hashlib.sha256(data).hexdigest()
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None and len(z.namelist())==len(entries)
        assert all(hashlib.sha256(z.read(name)).hexdigest()==digest for name,digest in entries.items())
    packages.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path),
        entries=len(entries),everyEntrySHA256=entries,CRCAndEveryEntryPassed=True,isolatedMATLABRecomputation=False))
out=dict(status='passed-current-state-research-archives-CRC-and-every-entry-SHA',recordedUTC=datetime.now(timezone.utc).isoformat(),
    processID=os.getpid(),progress=manifest.relative_to(ROOT).as_posix(),progressSHA256=sha(manifest),packages=packages,
    formalVersion83DeliveryClaimed=False,largeMATCheckpointsIncluded=False,isolatedMATLABRecomputation=False,
    newProductionSolvedAreaM2=0,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],packages=[{k:v for k,v in p.items() if k!='everyEntrySHA256'} for p in packages])))
