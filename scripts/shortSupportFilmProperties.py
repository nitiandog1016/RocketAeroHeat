"""把每个等温候选的膜温、空气物性及不可压原生配置绑定到同一试算。"""
import hashlib
import math
from shortSphere65 import CASE, condition, configuration, kelvin, read

# v340 的原生配置身份，统一 LF。更换物性时仅流动热导、被动标量扩散和停止输出时刻可变化。
BASE_CONFIG={(200,False):"45d962c158cdac36f3c9a85c5a2c99cee345653d792ab77e7f565f35500e5260",
             (20,False):"8794149ec83c5038a900c5fcdb1147a851d3f7f5e9cb1474b39875cc61d95566",
             (1,False):"1fac0dedb8a5098b2dbb8361fe517075dc496b71c57084002f66cff3bca4145a",
             (1,True):"f0a24e3078dd2d32f08ca68a663440428cccc7a25f2596398fa58b528968bcf5"}
OLD_DIFF="0.0016366381190378977"


def film_state(case,excess):
    """只使用来流、试件和候选超温；可直接传入没有表温和七站测量列的字典。"""
    if not math.isfinite(excess) or not 20<=excess<=40:
        raise ValueError("本试件的登记试算超温为 20–40 K。")
    air=float(kelvin(case["airF"]));wall=air+excess
    c=condition(case,wall_k=wall)
    return c,dict(model="air-conductivity-at-candidate-film",airK=air,wallK=wall,filmK=(air+wall)/2,
                  kWmK=c["kFilm"],pe=c["pe"],diffusionND=1/c["pe"])


def native_text(spec,steps,release):
    c,identity=film_state(read(CASE),spec["uniformExcessK"])
    assert spec["condition"]==c and spec["filmTransport"]==identity
    text=configuration(c,steps,release,"BOUNDED_SCALAR")
    baseline=text
    for key in ("THERMAL_CONDUCTIVITY_CONSTANT","DIFFUSIVITY_CONSTANT"):
        line=f"{key}= {identity['diffusionND']}\n"
        assert text.count(line)==1
        baseline=baseline.replace(line,f"{key}= {OLD_DIFF}\n")
    assert hashlib.sha256(baseline.encode("ascii")).hexdigest()==BASE_CONFIG[steps,release]
    # 保存固定长度的历史；科学接受仍由六分量末窗、无限幅方程和热量账决定。
    assert text.count("CONV_STARTITER= 50\n")==1
    return text.replace("CONV_STARTITER= 50\n",f"CONV_STARTITER= {steps+2}\n")
