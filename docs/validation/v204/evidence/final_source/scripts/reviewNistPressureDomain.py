"""模块：绑定现有CoolProp与一手资料，区分参考压力、模型域和混合物身份。"""
from pathlib import Path
import hashlib,json,sys,time,shutil

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'tools/runtime/su2-thermal-build/source/subprojects/CoolProp'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    start=time.monotonic();out=ROOT/'outputs/v204_validation/source02';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    sys.path.insert(0,str(ROOT/'tools/runtime/coolprop-7.2.0'))
    import CoolProp
    from CoolProp.CoolProp import get_global_param_string
    assert CoolProp.__version__=='7.2.0'
    assert get_global_param_string('gitrevision')=='98b3523d5daa98454618d381d2ae53f7471d216b'
    fluids={}
    files=[]
    for name in ['Nitrogen','Oxygen','Argon','Air']:
        p=SRC/'dev/fluids'/f'{name}.json';d=json.loads(p.read_text());e=d['EOS'][0]
        fluids[name]={k:e.get(k) for k in ['BibTeX_EOS','T_max','p_max','gas_constant','molar_mass']};files.append(p)
    p=SRC/'dev/mixtures/mixture_binary_pairs.json';pairs=json.loads(p.read_text())
    pairs=[d for d in pairs if {d['CAS1'],d['CAS2']} <= {'7727-37-9','7782-44-7','7440-37-1'}];assert len(pairs)==3;files.append(p)
    files += [SRC/'CoolPropBibTeXLibrary.bib', SRC/'LICENSE',SRC/'src/Backends/Helmholtz/FlashRoutines.cpp']
    copies=[]
    for p in files:
        q=out/'fixed_source'/p.relative_to(SRC);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        copies.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    from pypdf import PdfReader
    paper=ROOT/'outputs/v153_validation/source_review/lemmon2000.pdf'
    reader=PdfReader(str(paper));pages=[]
    for i,page in enumerate(reader.pages):
        txt=page.extract_text() or ''
        if any(s in txt.lower() for s in ['extrapolat','2000 k','1000 k']):
            pages.append(dict(pdfPage=i+1,text=txt))
    # 原文仅留本地阅读；公开证据只保存来源位置、SHA与本项目自己的结论。
    (out/'lemmon_selected_pages.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2),encoding='utf-8')
    r=dict(passed=True,python=sys.version,version=CoolProp.__version__,revision=get_global_param_string('gitrevision'),
        binarySHA256={p.name:sha(p) for p in (ROOT/'tools/runtime/coolprop-7.2.0/CoolProp').glob('*.pyd')},
        fluids=fluids,binaryPairs=pairs,files=copies,paper=dict(path=paper.relative_to(ROOT).as_posix(),sha256=sha(paper),pages=len(reader.pages)),
        conclusions=['100000Pa standard entropy reference is not an accuracy ceiling',
        'Lemmon2000 standard Air validity to2000K is distinct from variable-mixture validity to1000K',
        'CoolProp7.2 actual mixture pairs cite Gernert2013, not identical to Lemmon2000 mixture equation',
        'Unimposed PT uses Gernert stability; returned iphase_liquid is hard-coded for any stable single-phase mixture',
        'No imposed-phase result alone proves phase stability; metadata Tmax alone does not prove accuracy'],
        elapsedSeconds=time.monotonic()-start,CFDRuns=0,builds=0,physicalQualification=False)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(passed=True,seconds=r['elapsedSeconds'],fluids=fluids,pages=[x['pdfPage'] for x in pages])))

if __name__=='__main__':main()
