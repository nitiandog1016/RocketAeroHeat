"""Draw the original RocketAeroHeat icon; no external artwork or network input."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "RocketAeroHeat-desktop-icon-v1"
import matplotlib.pyplot as plt
from matplotlib.path import Path as VectorPath
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch, Polygon
from PIL import Image


def create_icon():
    output = Path(__file__).resolve().parents[1] / "assets" / "icons"
    output.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(4, 4), dpi=256)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set(xlim=(0, 256), ylim=(256, 0), aspect="equal")
    ax.axis("off")

    def curve(points, color, width):
        path = VectorPath(points, [VectorPath.MOVETO] + [VectorPath.CURVE4] * 3)
        ax.add_patch(PathPatch(path, fill=False, edgecolor=color,
                               linewidth=width, capstyle="round"))

    ax.add_patch(FancyBboxPatch((8, 8), 240, 240, boxstyle="round,pad=0,rounding_size=50",
                              facecolor="#102c46", edgecolor="none"))
    # Two swept airflow contours remain recognizable at desktop sizes.
    curve([(31, 151), (34, 94), (84, 39), (145, 31)], "#21627c", 4)
    curve([(106, 225), (170, 221), (221, 167), (225, 112)], "#28768a", 4)
    # Flame and core, behind the angled rocket.
    flame = VectorPath([(95, 142), (48, 147), (54, 178), (37, 217),
                        (78, 202), (112, 204), (119, 163), (95, 142)],
                       [VectorPath.MOVETO] + [VectorPath.CURVE4] * 6 + [VectorPath.CLOSEPOLY])
    ax.add_patch(PathPatch(flame, facecolor="#ff7947", edgecolor="none"))
    core = VectorPath([(99, 157), (73, 160), (81, 181), (64, 192),
                       (89, 185), (111, 186), (111, 165), (99, 157)],
                      [VectorPath.MOVETO] + [VectorPath.CURVE4] * 6 + [VectorPath.CLOSEPOLY])
    ax.add_patch(PathPatch(core, facecolor="#ffd477", edgecolor="none"))
    ax.add_patch(Polygon([(105, 106), (68, 124), (62, 158), (105, 149)],
                         closed=True, facecolor="#56c6d8", edgecolor="none"))
    ax.add_patch(Polygon([(148, 150), (135, 192), (168, 185), (185, 147)],
                         closed=True, facecolor="#56c6d8", edgecolor="none"))
    body = VectorPath([(95, 140), (112, 91), (158, 53), (207, 48),
                       (205, 99), (170, 146), (120, 166), (95, 140)],
                      [VectorPath.MOVETO] + [VectorPath.CURVE4] * 6 + [VectorPath.CLOSEPOLY])
    ax.add_patch(PathPatch(body, facecolor="#f4f9fb", edgecolor="none"))
    ax.add_patch(Polygon([(99, 142), (120, 163), (129, 158), (104, 133)],
                         closed=True, facecolor="#b3d8e2", edgecolor="none"))
    ax.add_patch(Circle((158, 100), 21, facecolor="#102c46", edgecolor="none"))
    ax.add_patch(Circle((158, 100), 14, facecolor="#5dd8e6", edgecolor="none"))
    ax.add_patch(Circle((153, 95), 5, facecolor="#d9fcff", edgecolor="none"))
    fig.savefig(output / "rocket-aero-heat.svg", transparent=True,
                metadata={"Date": None, "Creator": "RocketAeroHeat"})
    svg = output / "rocket-aero-heat.svg"
    clean_svg = "\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()) + "\n"
    with svg.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(clean_svg)
    fig.savefig(output / "rocket-aero-heat.png", transparent=True)
    plt.close(fig)
    with Image.open(output / "rocket-aero-heat.png") as source:
        source.resize((256, 256), Image.Resampling.LANCZOS).save(
            output / "rocket-aero-heat.ico", format="ICO",
            sizes=[(size, size) for size in (16, 24, 32, 48, 64, 128, 256)])


if __name__ == "__main__":
    create_icon()
