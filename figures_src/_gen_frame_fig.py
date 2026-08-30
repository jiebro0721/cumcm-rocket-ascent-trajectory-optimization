"""生成论文 §5.1 用图：坐标系与方向规定示意图（教学用，最终版）。

左图：三维地心惯性系（x 轴指向发射点经线、y 轴东、z 轴自转轴、发射点位置、
      初始速度方向、自转右手方向）
右图：赤道平面俯视图（发射点、飞行中火箭、rhat/that 局部基、phi/gamma 角）
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

OUT = r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2\论文"


def arrow(ax, p0, p1, color, lw=1.8, ms=13):
    ax.annotate("", xy=p1, xytext=p0,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                mutation_scale=ms))


def fig_frame():
    fig = plt.figure(figsize=(13.5, 6.4))

    # ================= 左图：三维坐标系 =================
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.set_box_aspect((1, 1, 0.85))
    ax.view_init(elev=22, azim=-58)

    Re = 1.0
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 40)
    xs = Re * np.outer(np.cos(u), np.sin(v))
    ys = Re * np.outer(np.sin(u), np.sin(v))
    zs = Re * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="#dbe7f3", alpha=0.42, linewidth=0)

    ax.plot([0, 1.55], [0, 0], [0, 0], color="#1F4E79", lw=2.2)
    ax.plot([0, 0], [0, 1.55], [0, 0], color="#b23c3c", lw=2.2)
    ax.plot([0, 0], [0, 0], [0, 1.30], color="#2E6B4E", lw=2.2)
    ax.text(1.60, -0.06, -0.10, r"$x$ 轴：指向发射点经线", color="#1F4E79", fontsize=13)
    ax.text(-0.10, 1.60, -0.06, r"$y$ 轴：当地东向", color="#b23c3c", fontsize=13)
    ax.text(-0.08, -0.08, 1.36, r"$z$ 轴：地球自转轴", color="#2E6B4E", fontsize=13)

    ax.text(0.06, 0.02, 0.12, r"$O$ 地心", color="#333", fontsize=14)

    ax.plot([Re], [0], [0], "o", ms=9, color="#e07b39")
    ax.text(Re + 0.04, 0.14, 0.02, r"发射点 $(R_E,0,0)$", color="#e07b39", fontsize=13)

    ax.plot([Re, Re], [0, 0.66], [0, 0], color="#b23c3c", lw=2.4, ls="--")
    ax.text(Re + 0.04, 0.74, 0.0, r"$v_0=\omega_E R_E$  东向初速",
            color="#b23c3c", fontsize=13)

    th = np.linspace(0.0, 1.55, 40)
    ax.plot(0.98 * np.cos(th), 0.98 * np.sin(th), 0.18 * np.zeros_like(th),
            color="#555", lw=1.6, alpha=0.8)
    ax.plot([0.98 * np.cos(1.55)], [0.98 * np.sin(1.55)], [0.18],
            marker=(3, 0, -75), color="#555", ms=9)
    ax.text(0.52, 0.98, 0.16, r"自转方向 $\omega_E$", color="#555", fontsize=13)

    ax.set_xlim(-0.2, 1.9)
    ax.set_ylim(-0.2, 1.9)
    ax.set_zlim(-1.0, 1.5)
    ax.set_axis_off()
    ax.set_title("a  地心惯性坐标系", fontsize=17)

    # ================= 右图：赤道平面 =================
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect("equal")
    ax2.set_xlim(-1.55, 2.45)
    ax2.set_ylim(-1.55, 2.15)
    ax2.axis("off")

    ax2.add_patch(Circle((0, 0), 0.9, fc="#e8edf3", ec="#7a8ea8", lw=1.4))
    ax2.plot(0, 0, "o", ms=5, color="#333")
    ax2.text(-0.18, -0.28, r"$O$ 地心", fontsize=14, color="#333")

    arrow(ax2, (-1.4, 0), (2.3, 0), "#1F4E79", 1.8)
    arrow(ax2, (0, -1.4), (0, 1.9), "#b23c3c", 1.8)
    ax2.text(2.38, -0.12, r"$x$（发射点经线方向）", color="#1F4E79", fontsize=13)
    ax2.text(0.06, 1.97, r"$y$（东向）", color="#b23c3c", fontsize=13)

    # 发射点与初速（标签放右下）
    ax2.plot(0.9, 0, "o", ms=9, color="#e07b39")
    arrow(ax2, (0.9, 0), (0.9, 0.72), "#b23c3c", 2.0, ms=13)
    ax2.text(1.06, 0.80, r"$v_0=\omega_E R_E$", color="#b23c3c", fontsize=14)
    ax2.text(1.06, 0.60, r"发射点 $(R_E,0)$", color="#e07b39", fontsize=13)

    # 飞行中的火箭（标签放左上）
    th = np.deg2rad(52)
    rc = 1.35
    rx, ry = rc * np.cos(th), rc * np.sin(th)
    ax2.plot(rx, ry, "o", ms=10, color="#1F4E79")
    ax2.text(rx - 0.36, ry + 0.48, r"火箭 $(x,y)$", color="#1F4E79",
             fontsize=14, ha="center")

    # 位置向量与局部基
    arrow(ax2, (0, 0), (rx, ry), "#333", 1.6, ms=11)
    ax2.text(0.52 * np.cos(th) - 0.32, 0.52 * np.sin(th) - 0.22,
             r"$r$", color="#333", fontsize=15)

    ux, uy = np.cos(th), np.sin(th)
    tx, ty = -np.sin(th), np.cos(th)
    arrow(ax2, (rx, ry), (rx + 0.50 * ux, ry + 0.50 * uy), "#1F4E79", 1.6, ms=11)
    ax2.text(rx + 0.44 * ux - 0.13, ry + 0.62 * uy, r"$\hat{r}$", color="#1F4E79", fontsize=15)
    arrow(ax2, (rx, ry), (rx + 0.50 * tx, ry + 0.50 * ty), "#1F4E79", 1.6, ms=11)
    ax2.text(rx + 0.70 * tx - 0.06, ry + 0.50 * ty, r"$\hat{t}$", color="#1F4E79", fontsize=15)

    # 当地水平面（过火箭的虚线）
    hl = 1.0
    ax2.plot([rx - hl * tx, rx + hl * tx], [ry - hl * ty, ry + hl * ty],
             ls="--", color="#a5a5a5", lw=1.0)

    # 速度 v 与推力 T
    gv, pv = np.deg2rad(18), np.deg2rad(30)
    vx, vy = np.cos(th + gv), np.sin(th + gv)
    fx, fy = np.cos(th + pv), np.sin(th + pv)
    arrow(ax2, (rx, ry), (rx + 0.80 * vx, ry + 0.80 * vy), "#2E6B4E", 2.0, ms=13)
    ax2.text(rx + 0.94 * vx - 0.04, ry + 0.76 * vy - 0.15, r"$v$", color="#2E6B4E", fontsize=15)
    arrow(ax2, (rx, ry), (rx + 0.85 * fx, ry + 0.85 * fy), "#b23c3c", 2.3, ms=14)
    ax2.text(rx + 0.98 * fx - 0.10, ry + 0.90 * fy, r"$T$", color="#b23c3c", fontsize=16)

    # 角 γ、φ
    a1 = np.linspace(0, gv, 30)
    ax2.plot(rx + 0.32 * np.cos(th + a1), ry + 0.32 * np.sin(th + a1),
             color="#2E6B4E", lw=1.2)
    ax2.text(rx + 0.50 * np.cos(th + gv * 0.45) - 0.20,
             ry + 0.50 * np.sin(th + gv * 0.45) - 0.20,
             r"$\gamma$", color="#2E6B4E", fontsize=15)
    a2 = np.linspace(0, pv, 30)
    ax2.plot(rx + 0.58 * np.cos(th + a2), ry + 0.58 * np.sin(th + a2),
             color="#b23c3c", lw=1.2)
    ax2.text(rx + 0.80 * np.cos(th + pv * 0.45) + 0.02,
             ry + 0.80 * np.sin(th + pv * 0.45) + 0.02,
             r"$\varphi$", color="#b23c3c", fontsize=15)

    # 底部说明（符号定义，字号 11 以保证清晰）
    ax2.text(-1.4, -1.42,
             r"$\hat{r}=(x/r,\ y/r)$ 径向，$\hat{t}=(-y/r,\ x/r)$ 东向，"
             r"$\hat{u}=\cos\varphi\,\hat{t}+\sin\varphi\,\hat{r}$",
             fontsize=11.5, color="#555")
    ax2.text(-1.4, -1.56, r"$v_r=v\cdot\hat{r}$，$v_t=v\cdot\hat{t}$；"
             r"$\gamma$: 速度与水平面夹角，$\varphi$: 推力与水平面夹角",
             fontsize=11.5, color="#555")

    ax2.set_title("b  赤道平面内：局部基与角的规定", fontsize=17)

    fig.tight_layout()
    fig.savefig(OUT + r"\fig_frame.pdf")
    plt.close(fig)
    print("fig_frame 完成")


if __name__ == "__main__":
    fig_frame()
