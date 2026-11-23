"""从已校验的论文渲染图读图 12 局部曲线，保留像素宽度和轴标定。"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha
from native_check import appendix


def measure(pdf,case,out):
    assert sha(pdf)=="41e9098cdfafea78a5a613a1e81c165d1cb910c6075f106b9fede8fb78b94976"
    out.mkdir(parents=True,exist_ok=False)
    subprocess.run(["pdftoppm","-f","9","-l","9","-r","320","-png",str(pdf),str(out/"page")],check=True)
    page=Image.open(out/"page-09.png").convert("L")
    w,h=page.size
    assert (w,h)==(2645,3743)
    box=(int(.48*w),int(.77*h),int(.9*w),int(.925*h))
    crop=page.crop(box);ink=np.array(crop)<120
    re_l=appendix(read(case))["reynoldsL"]
    # 实测点在十至百的对数十年内。局部垂直方向分别用零线及 0.02 水平刻度标定。
    x=299+(482-299)*np.log10(re_l/10)
    readings=[]
    for i in range(round(x)-2,round(x)+3):
        def stroke(lo,hi):
            y=np.flatnonzero(ink[lo:hi,i])+lo
            assert len(y)>0 and np.all(np.diff(y)==1)
            return [int(y.min()),int(y.max())]
        baseline,grid,curve=stroke(331,343),stroke(276,286),stroke(298,320)
        def value(z,g,c):return .02*(z-c)/(z-g)
        # 逐列枚举墨迹外缘并加入半像素量化，避免把线宽当作模型或仪器不确定度。
        bounds=[value(z,g,c) for z in (baseline[0]-.5,baseline[1]+.5) for g in (grid[0]-.5,grid[1]+.5)
                for c in (curve[0]-.5,curve[1]+.5)]
        readings.append(dict(x=i,baseline=baseline,grid=grid,curve=curve,
                             center=value(np.mean(baseline),np.mean(grid),np.mean(curve)),interval=[min(bounds),max(bounds)]))
    data=dict(pdfSHA256=sha(pdf),page=818,figure=12,dpi=320,cropBox=list(box),pagePixels=[w,h],xAtRe10=299,xAtRe100=482,
              ratioAtGrid=.02,inkThreshold=120,xAtPoint=x,reynoldsL=re_l,columns=readings,
              value=float(np.median([r["center"] for r in readings])),range=[min(r["interval"][0] for r in readings),max(r["interval"][1] for r in readings)])
    (out/"reading.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("pdf","case","out"):p.add_argument(k,type=Path)
    a=p.parse_args();measure(a.pdf.resolve(),a.case.resolve(),a.out.resolve())
