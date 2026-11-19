"""用项目的 CoolProp 7.2.0 Air 独立核对给定近室温工况的干空气物性。"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[3]
PACKAGE=ROOT/"tools/runtime/coolprop-7.2.0"
sys.path.insert(0,str(PACKAGE))
import CoolProp
from CoolProp.CoolProp import PropsSI


def check():
    assert CoolProp.__version__=="7.2.0"
    peaks={};worst={}
    # 多留 1 K 包含新工况恢复温升对参考温度的影响。
    for t in range(280,312):
        mu=1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)
        for p in range(90000,110001,1000):
            values={"D":p/(287.05*t),"C":1005.,"V":mu,"L":mu*1005/.7}
            for key,value in values.items():
                ref=PropsSI(key,"T",t,"P",p,"Air")
                error=abs(value/ref-1)
                if error>peaks.get(key,-1):
                    peaks[key]=error;worst[key]=[t,p,value,ref]
    return dict(version=CoolProp.__version__,gridT=[280,311,1],gridPa=[90000,110000,1000],
                maxRelative=peaks,worst=worst,extensionSHA256={f.name:hashlib.sha256(f.read_bytes()).hexdigest()
                for f in (PACKAGE/"CoolProp").glob("*.pyd")})


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_output",type=Path)
    result=check()
    with parser.parse_args().new_output.open("x",encoding="utf8") as stream:
        json.dump(result,stream,indent=2)
    assert max(result["maxRelative"].values())<.03
