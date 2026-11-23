"""由发布小数据与归档轨迹绘制越界节点数和首次局部残差，不读取大场。"""
import io
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from verify_flux import DATA, HERE, read


if __name__ == "__main__":
    d = read(DATA/"supportedFlux65.json")
    with ZipFile(HERE/"evidence/receipts.zip") as archive:
        trace = np.loadtxt(io.BytesIO(archive.read("outputs/v326_validation/detail01/extrema.csv")), delimiter=",")
    plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.spines.top": False, "axes.spines.right": False})
    fig, (timeline, bars) = plt.subplots(1, 2, figsize=(10.8, 3.8), layout="constrained")
    timeline.plot(trace[:, 0]+1, trace[:, 4], color="#246487", label="质量通量对照 · 每次更新")
    frames = d["frames"]
    timeline.scatter([f["index"]+1 for f in frames], [f["baselineBelow"] for f in frames],
                     marker="s", facecolor="white", edgecolor="#3a4145", label="原迎风格式 · 已保存帧")
    step = d["first"]["update"]
    timeline.scatter([step], [len(d["first"]["belowIDs"])], color="#aa5029")
    timeline.annotate("首次越界 · 6 个节点", (step, len(d["first"]["belowIDs"])), (20, 12), textcoords="offset points")
    timeline.set(xlabel="稳态更新次数", ylabel=r"低于 $-10^{-12}$ 的流体节点数", title="全域温差标量")
    timeline.legend(loc="lower right", fontsize=8, frameon=False)
    keys = ["advection", "massDivergence", "twoPoint", "gradientCorrection", "spatial"]
    values = np.array([d["row"][k] for k in keys])*1e14
    bars.barh(["对流", "质量散度修正", "两点扩散", "扩散梯度修正", "装配空间残差"], values,
              color=["#567888", "#567888", "#567888", "#ab5634", "#36464c"])
    bars.invert_yaxis();bars.axvline(0, color="#333333", linewidth=.7)
    bars.set(xlabel=r"空间残差贡献（原生无量纲量，$\times 10^{-14}$）", title="第 8 次更新 · 流体节点 78423")
    fig.savefig(HERE/"evidence/flux.png", dpi=160)
    plt.close(fig)
