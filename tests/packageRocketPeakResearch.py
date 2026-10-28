# -*- coding: utf-8 -*-
"""冻结局部峰值研究的完整源码/替换文件与离线预览，不冒充完整隔离工程交付。"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,subprocess,zipfile

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
verification=VAL/'peak_research_package_verification.json';assert not verification.exists()
manifest=ROOT/'docs/validation/v83/progress_peak_weighted_direction.json';m=json.loads(manifest.read_text(encoding='utf8'))
assert m['production']['goalStatus']=='active' and not m['production']['version83IsolatedCompleteEngineeringDelivery']


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


names=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT).split(b'\0')
source=sorted({ROOT/Path(n.decode('utf8')) for n in names if n and (ROOT/Path(n.decode('utf8'))).is_file()})
prior=json.loads((ROOT/'docs/validation/v83/progress_current_state_direction.json').read_text(encoding='utf8'))
replace=sorted({ROOT/n for n in set(prior['currentSourcesSHA256'])|set(m['currentSourcesSHA256'])}|{manifest,Path(__file__).resolve()})
selected=VAL/('workbench38' if m['selectedEndpointCase']=='krylov24' else 'workbench35')
preview=[(selected,'selected-preview'),(VAL/'workbench36','coordinate-study-preview')]
readme=f'''RocketAeroHeat v83 局部峰值约束研究（2026-10-26）

研究源码快照/完整替换文件；本包没有隔离MATLAB完整工程重算，不是正式v83交付。
最新正式交付仍为v82，正式合格热流覆盖66.5626327424%，本轮新增0。

源码包保留当前非忽略的源码、测试、模型和文档。main.m为工程入口。
替换文件均为完整文件，需要既有RocketAeroHeat工程。MATLAB后端/X主轴/鼻尖+X/来流-X/SI保持。
selected-preview/index.html展示当前选用端点{m['selectedEndpointCase']}。
coordinate-study-preview/index.html展示同Q20、同方向的Q23候选坐标对照。
Q22和Q23是并列子分支，不能合并累计接受和NS调用。所有研究值仍未具区域物理资格。
Tw=300K是输入，Te/Tr/求解Tw未知。数据/科研PNG检查不是浏览器/WebGL视觉验收。

数值研究依赖本机outputs中的大型MAT检查点和公开参考缓存，本包不包含这些大型依赖。
ZIP CRC和逐文件SHA通过不代表隔离物理解重算或实验验证。
详见docs/PEAK_WEIGHTED_DIRECTION_V83.md与docs/validation/v83/progress_peak_weighted_direction.json。
'''
packages=[]
for name,files,prefix in [('RocketAeroHeat_v83_peak_research_source.zip',source,'RocketAeroHeat/'),
                          ('RocketAeroHeat_v83_peak_complete_replacements.zip',replace,'')]:
    path=VAL/name;assert not path.exists();entries={}
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in files:
            f.resolve().relative_to(ROOT.resolve());entry=prefix+f.relative_to(ROOT).as_posix();assert entry not in entries
            data=f.read_bytes();z.writestr(entry,data);entries[entry]=hashlib.sha256(data).hexdigest()
        for folder,label in preview:
            for f in sorted(p for p in folder.rglob('*') if p.is_file()):
                entry=prefix+label+'/'+f.relative_to(folder).as_posix();assert entry not in entries
                data=f.read_bytes();z.writestr(entry,data);entries[entry]=hashlib.sha256(data).hexdigest()
        entry=prefix+'RESEARCH_PACKAGE_README.txt';data=readme.encode('utf8');z.writestr(entry,data);entries[entry]=hashlib.sha256(data).hexdigest()
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None and len(z.namelist())==len(entries)
        assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in entries.items())
    packages.append(dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),bytes=path.stat().st_size,
        entries=len(entries),allEntrySHA256=entries,CRCAndEveryEntryPassed=True,isolatedMATLABRecomputation=False))
r=dict(status='passed-peak-research-source-and-complete-file-archives-CRC-and-all-SHA',recordedUTC=datetime.now(timezone.utc).isoformat(),
       progressSHA256=sha(manifest),packages=packages,formalVersion83DeliveryClaimed=False,newProductionSolvedAreaM2=0,scriptSHA256=sha(__file__))
verification.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=r['status'],packages=[{k:v for k,v in p.items() if k!='allEntrySHA256'} for p in packages])))
