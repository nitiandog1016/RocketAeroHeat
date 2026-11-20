"""核对两张一手图的坐标与单位。图读范围、稳定性拟合线和实验精度分开保存。"""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
ORIGINAL="41e9098cdfafea78a5a613a1e81c165d1cb910c6075f106b9fede8fb78b94976"
CROP="c4c81c99f21ab0b42037d09aaadead1aa28b80e387ab0e6dc238c4e23f16daed"
FIG5="9ca4bad73465ff6d2e5a534fb3e5d5ba064a3b57a5a847636629592983e491d5"
SVG_URL="https://arxiv.org/html/2504.17925v1/Rec_vs_AR_new-wide.svg"


def reading(center, zero, two_hundred):
    """图7的刻度定位沿用 v282，保留图读范围的全部端点组合。"""
    def scale(c,z,h):
        return (c-z)*200/(h-z)
    ends=[scale(center+a,zero+b,two_hundred+c)
          for a,b,c in product((-2,2),(-.75,.75),(-.75,.75))]
    return scale(center,zero,two_hundred),[min(ends),max(ends)]


def check(svg):
    source=ROOT/"outputs/v282_validation/source01"
    for path,expected in ((source/"tsubouchi1959.pdf",ORIGINAL),(source/"figure7.png",CROP),(svg,FIG5)):
        assert hashlib.sha256(path.read_bytes()).hexdigest()==expected, path
    dt,dt_range=reading(117.5,92.5,588)
    h_raw,h_range=reading(239,363.5,264.5)
    h=h_raw*4184/3600
    fig=ET.fromstring(svg.read_bytes())
    # 图5的绿色实线对应圆形长椭球标记。横轴每单位 62，纵轴每 500 对应 80。
    curves=[p for p in fig.iter() if p.get("stroke")=="#00b300" and len(p.get("d",""))>2000]
    assert len(curves)==1
    xy=np.asarray([float(t) for t in re.findall(r"[-+]?[\d.]+(?:[eE][-+]?\d+)?",curves[0].get("d"))]).reshape(-1,2)
    aspect=(1-.53**2)**(-.5)
    def critical(a):
        return float((271-np.interp(73+62*(a-1),xy[:,0],xy[:,1]))*500/80)
    # 用两处论文表2直接印出的数值验证曲线身份和坐标轴。
    assert abs(critical(1)-212.67)<1 and abs(critical(4)-999.66)<5
    return dict(heatExperiment=dict(doi="10.1299/kikai1938.25.810",figure=7,printedPage=815,
                pdfSHA256=ORIGINAL,cropSHA256=CROP,pointPixels=[117.5,239],speedMps=6,
                deltaK=dt,deltaReadRangeK=dt_range,hKcalM2HourK=h_raw,hWm2K=h,
                graphReadRangeWm2K=[x*4184/3600 for x in h_range],
                experimentalPointUncertainty=None,measuredQuantity="whole-bead mean h"),
                stability=dict(doi="10.1017/jfm.2025.10402",authorVersionFigure=5,
                svgURL=SVG_URL,svgSHA256=FIG5,axialAspectRatio=aspect,
                criticalReFigureCurve=critical(aspect),curveCheckAt1=critical(1),curveCheckAt4=critical(4),
                interpretation="regular m=1 primary bifurcation, curve reading, not Hopf or a frequency"))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("svg",type=Path)
    parser.add_argument("result",type=Path)
    args=parser.parse_args()
    result=check(args.svg)
    with args.result.open("x",encoding="utf8") as f:
        json.dump(result,f,indent=2)
        f.write("\n")
    print(json.dumps(result))
