"""生成问题三/四的实际轨迹与控制律曲线（与问题一风格一致，用于论文）。

所有路径基于 __file__ 相对定位，可在任意目录运行。
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = str(ROOT / "论文")
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

from q34_direct_collocation import DirectCollocationSolver  # noqa: E402


def plot_profiles(result, question, t2, t3, fname, prop_used):
    """h / V / gamma / m 四联图（与问题一一致）。"""
    t = result.t
    h = (np.hypot(result.x[:, 0], result.x[:, 1]) - 6371e3) / 1e3
    v = np.hypot(result.x[:, 2], result.x[:, 3])
    gamma = np.degrees(np.arctan2(
        result.x[:, 2] * result.x[:, 0] + result.x[:, 3] * result.x[:, 1],
        result.x[:, 2] * (-result.x[:, 1]) + result.x[:, 3] * result.x[:, 0],
    ))
    m = result.x[:, 4] / 1e3

    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5))
    axes[0, 0].plot(t / 60, h, lw=1.8, color="#1F4E79")
    axes[0, 0].axhline(400, ls="--", color="#b23c3c", lw=1.2, label="目标轨道 400 km")
    axes[0, 0].axvline(t2 / 60, ls=":", color="k", lw=0.8, label="二子级点火")
    axes[0, 0].set_xlabel("时间 t / min", fontsize=11); axes[0, 0].set_ylabel("高度 h / km", fontsize=11)
    axes[0, 0].legend(fontsize=8); axes[0, 0].grid(alpha=0.3)
    axes[0, 0].set_title("a 高度", fontsize=12)

    axes[0, 1].plot(t / 60, v / 1e3, lw=1.8, color="#b23c3c")
    axes[0, 1].axhline(7.6726, ls="--", color="gray", lw=1.2, label="圆轨道速度 7.6726 km/s")
    axes[0, 1].axvline(t2 / 60, ls=":", color="k", lw=0.8)
    axes[0, 1].set_xlabel("时间 t / min", fontsize=11); axes[0, 1].set_ylabel("惯性速度 V / km/s", fontsize=11)
    axes[0, 1].legend(fontsize=8); axes[0, 1].grid(alpha=0.3)
    axes[0, 1].set_title("b 惯性速度", fontsize=12)

    axes[1, 0].plot(t / 60, gamma, lw=1.8, color="#2E6B4E")
    axes[1, 0].axhline(0, ls="--", color="gray", lw=1.2)
    axes[1, 0].axvline(t2 / 60, ls=":", color="k", lw=0.8)
    axes[1, 0].set_xlabel("时间 t / min", fontsize=11); axes[1, 0].set_ylabel("飞行路径角 γ / deg", fontsize=11)
    axes[1, 0].grid(alpha=0.3)
    axes[1, 0].set_title("c 飞行路径角", fontsize=12)

    axes[1, 1].plot(t / 60, m, lw=1.8, color="#7a4d8f")
    axes[1, 1].axvline(t2 / 60, ls=":", color="k", lw=0.8)
    axes[1, 1].set_xlabel("时间 t / min", fontsize=11); axes[1, 1].set_ylabel("总质量 m / t", fontsize=11)
    axes[1, 1].grid(alpha=0.3)
    axes[1, 1].set_title("d 总质量", fontsize=12)

    fig.suptitle("问题%s 最优轨迹（滑行 %.3f s + 燃烧 %.3f s，耗药 %.1f kg）"
                 % (question, t2 - 149.06, t3 - t2, prop_used),
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT + "/" + fname, dpi=150)
    plt.close(fig)
    print("已生成", fname)


def plot_controls(result, question, t2, t3, phi_deg, sigma, fname):
    """Q3：仅俯仰角控制律（σ 恒为满推力，非优化量，不画）；
    Q4：俯仰角 + 节流比。"""
    tau_nodes = np.linspace(0, 1, len(phi_deg))
    t_burn = np.linspace(t2, t3, 300)
    tau = (t_burn - t2) / (t3 - t2)
    phi_curve = np.interp(tau, tau_nodes, phi_deg)

    if question == 3:
        fig, ax = plt.subplots(figsize=(7.5, 4))
        ax.plot(t_burn / 60, phi_curve, lw=2.2, color="#1F4E79")
        ax.set_xlabel("时间 t / min", fontsize=12)
        ax.set_ylabel("俯仰角 φ / deg", fontsize=12)
        ax.grid(alpha=0.3)
        ax.set_title("a 二子级俯仰角控制律", fontsize=13)
        ax.annotate("极小值 1.75°", xy=(7.5, 1.75), xytext=(5.5, 3.2),
                    fontsize=10, color="#1F4E79",
                    arrowprops=dict(arrowstyle="->", color="#1F4E79", lw=1))
        fig.tight_layout()
        fig.savefig(OUT + "/" + fname, dpi=150)
        plt.close(fig)
        print("已生成", fname)
        return

    sig_curve = np.interp(tau, tau_nodes, sigma)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(t_burn / 60, phi_curve, lw=2.0, color="#1F4E79")
    axes[0].set_xlabel("时间 t / min", fontsize=11)
    axes[0].set_ylabel("俯仰角 φ / deg", fontsize=11)
    axes[0].grid(alpha=0.3)
    axes[0].set_title("a 二子级俯仰角控制律", fontsize=12)

    axes[1].plot(t_burn / 60, sig_curve, lw=2.0, color="#b23c3c")
    axes[1].axhline(1.0, ls="--", color="gray", lw=1, label="额定推力 1.0")
    axes[1].axhline(0.6, ls="--", color="gray", lw=1, label="下限 0.6")
    axes[1].set_xlabel("时间 t / min", fontsize=11)
    axes[1].set_ylabel("节流比 σ", fontsize=11)
    axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3)
    axes[1].set_title("b 节流比控制律", fontsize=12)
    fig.suptitle("问题%s 最优控制律" % question, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT + "/" + fname, dpi=150)
    plt.close(fig)
    print("已生成", fname)


if __name__ == "__main__":
    # 运行 Q3 配点（20/40）并画图
    solver3 = DirectCollocationSolver(question=3, n_coast=20, n_burn=40)
    sol3 = solver3.solve()
    res3, val3 = solver3.validate(sol3)
    t2_3 = solver3.rk.t_burn1 + sol3.coast_duration
    t3_3 = t2_3 + sol3.burn_duration
    plot_profiles(res3, 3, t2_3, t3_3, "fig_q3_trajectory.pdf", sol3.propellant_used)
    plot_controls(res3, 3, t2_3, t3_3, sol3.phi_deg, sol3.sigma, "fig_q3_controls.pdf")

    # 运行 Q4 配点（20/40）并画图
    solver4 = DirectCollocationSolver(question=4, n_coast=20, n_burn=40)
    sol4 = solver4.solve()
    res4, val4 = solver4.validate(sol4)
    t2_4 = solver4.rk.t_burn1 + sol4.coast_duration
    t3_4 = t2_4 + sol4.burn_duration
    plot_profiles(res4, 4, t2_4, t3_4, "fig_q4_trajectory.pdf", sol4.propellant_used)
    plot_controls(res4, 4, t2_4, t3_4, sol4.phi_deg, sol4.sigma, "fig_q4_controls.pdf")

    print("全部图表生成完成")
