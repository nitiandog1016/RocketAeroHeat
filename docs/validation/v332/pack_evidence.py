"""归档本轮小型运行回执、边界核对和测站图；大场只登记哈希及本机位置。"""
import hashlib
import io
import json
from zipfile import ZipFile, ZIP_DEFLATED
import numpy as np
from sector_grid import HERE,ROOT,read,sha


def pack(base,ui,out):
    files=[base/"grid.json",base/"prepared01/prepared.json",base/"grid01_process.json",base/"prepare01_process.json",
           base/"audit.json",base/"post_process.json",base/"audit_process.json",ui/"ui.json",ui/"export_check.json",ui/"debugger_process.json",
           *list((base/"build01").glob("*.json"))]
    work=base/"run01";files.extend(work.glob("*.json"));large={}
    for stage in work.iterdir():
        if not stage.is_dir():continue
        for path in stage.iterdir():
            if path.name in ("record.json","history.csv","stages.csv","far_pairs.csv","case.cfg","flow_far.csv","debugger.log","debugger_process.json"):
                files.append(path)
            if path.name.startswith(("restart_","surface_")) or path.suffix==".bin":
                large[path.relative_to(base).as_posix()]=sha(path)
    # 加入本机同态的壁节点原始行，让只取 Git 版本的读者也能独立复核两套功率。
    from thermistorWireSurface import residual
    proof={}
    for tag,rel in (("early","releaseEarly"),("late","releaseLate")):
        n=np.fromfile(work/tag/"nodes.bin","<f8").reshape(-1,16);n=n[n[:,0].argsort()]
        wet=np.flatnonzero(n[:,4]);far=np.loadtxt(work/tag/"flow_far.csv",delimiter=",");ids=far[:,0].astype(int)
        proof[tag+"_wall"]=np.c_[wet,residual(work/rel/"restart.csv")[wet],n[wet,9:11]]
        proof[tag+"_far"]=np.c_[ids,n[ids,1],far[:,1],n[ids,11]]
    with ZipFile(out,"x",compression=ZIP_DEFLATED,compresslevel=9) as z:
        packed={}
        for file in files:
            name=file.relative_to(base).as_posix();raw=file.read_bytes();z.writestr(name,raw);packed[name]=hashlib.sha256(raw).hexdigest()
        for name,values in proof.items():
            tmp=io.StringIO();np.savetxt(tmp,values,delimiter=",",fmt=["%d","%.17g","%.17g","%.17g"])
            raw=tmp.getvalue().encode("ascii");z.writestr(name+".csv",raw);packed[name+".csv"]=hashlib.sha256(raw).hexdigest()
        z.writestr("manifest.json",json.dumps(dict(included=packed,largeLocal=large),indent=2))
    from verify_result import check_receipts
    check_receipts(out,read(work/"thermal.json"))


if __name__=="__main__":
    import sys
    base=ROOT/"outputs/v332_validation";pack(base,base/sys.argv[1],HERE/"evidence/receipts.zip")
