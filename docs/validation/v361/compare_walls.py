"""只用已发布的两份球壁表，重算压力与速度误差并生成角度对照。"""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def compare(plot_path=None):
    series = []
    for name in ("sphere", "sphere_fine"):
        data = np.genfromtxt(HERE / f"{name}_wall.csv", delimiter=",", names=True)
        theta = np.arctan2(data["y"], data["x"])
        i = np.argsort(theta)
        theta = theta[i]
        sine, cosine = np.sin(theta), np.cos(theta)
        velocity = -sine * data["Velocity_x"][i] + cosine * data["Velocity_y"][i]
        pressure = 2 * data["Pressure"][i]
        theory_u = 1.5 * sine
        theory_p = 1 - theory_u**2
        band = -np.diff(cosine)
        weights = (np.r_[0, band] + np.r_[band, 0]) / 2
        assert np.all(np.diff(theta) > 0) and abs(weights.sum()-2) < 1e-12
        ref = json.loads((HERE / f"{name}.json").read_bytes())
        errors = [np.sqrt(np.average((pressure-theory_p)**2, weights=weights)),
                  np.max(np.abs(pressure-theory_p)), np.sqrt(np.average((velocity-theory_u)**2, weights=weights))]
        np.testing.assert_allclose(errors, [ref["pressureCoefficientRMS"], ref["pressureCoefficientMax"],
                                            ref["tangentialVelocityRMS"]], rtol=1e-12)
        series.append((np.rad2deg(theta), pressure-theory_p, velocity-theory_u))
    np.testing.assert_allclose(series[0][0], series[1][0][::2], atol=1e-10)
    if plot_path is not None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7.0, 5.6), constrained_layout=True)
        for (angle, cp_error, speed_error), label in zip(series, ("64 × 36", "128 × 72")):
            axes[0].plot(angle, cp_error, label=label)
            axes[1].plot(angle, speed_error)
        axes[0].set(title="Inviscid sphere: numerical minus analytical wall values",
                    ylabel=r"Pressure error, $\Delta C_p$")
        axes[0].legend(frameon=False)
        axes[1].set(ylabel=r"Tangential speed error / $U_\infty$", xlabel="Angle from front stagnation point (deg)")
        for ax in axes:
            ax.axhline(0, color="#666666", linewidth=0.6)
            ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
            ax.grid(axis="y", alpha=0.2)
        fig.savefig(plot_path, dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", type=Path)
    compare(parser.parse_args().figure)
    print("Both published wall tables agree with the full-field readback.")
