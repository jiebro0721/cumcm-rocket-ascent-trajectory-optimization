"""生成论文示意图 v2：更大画布、更大字号、布局合理（对齐 A127 图标准）。"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["mathtext.fontset"] = "stix"

OUT = r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2\论文"

# 统一字号
FS_TITLE = 16
FS_LABEL = 14
FS_TEXT = 13
FS_BOX = 12.5


def _set_default(ax):
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(FS_TEXT)


# ===========================================================================
# 图 1：受力分析与三阶段控制
# ===========================================================================
def fig_forces_coords():
    fig, ax1 = plt.subplots(figsize=(7.8, 7.0))

    # 受力图（放大、仅符号——不含完整数学公式）
    ax1.set_xlim(-2.35, 2.35)
    ax1.set_ylim(-2.35, 2.35)
    ax1.set_aspect("equal")
    ax1.axis("off")

    earth = plt.Circle((0, 0), 0.80, fc="#e8edf3", ec="#7a8ea8", lw=1.4)
    ax1.add_patch(earth)
    ax1.plot(0, 0, "o", ms=5, color="#555")
    ax1.text(0.045, -0.09, r"$O$", fontsize=15, ha="left", va="top")

    th = np.deg2rad(42)
    rc = 1.55
    rx, ry = rc * np.cos(th), rc * np.sin(th)
    ax1.plot(rx, ry, "o", ms=11, color="#1F4E79")

    ax1.annotate("", xy=(rx, ry), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="-", color="#333", lw=1.8))
    ax1.text(0.45 * np.cos(th * 0.55), 0.45 * np.sin(th * 0.55) - 0.08,
             r"$r$", fontsize=16, color="#333")

    ux, uy = np.cos(th), np.sin(th)
    tx, ty = -np.sin(th), np.cos(th)
    ax1.annotate("", xy=(rx + 0.70 * ux, ry + 0.70 * uy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#1F4E79", lw=1.5, mutation_scale=13))
    ax1.text(rx + 0.82 * ux, ry + 0.60 * uy, r"$\hat{r}$", fontsize=13, color="#1F4E79")
    ax1.annotate("", xy=(rx + 0.70 * tx, ry + 0.70 * ty), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#1F4E79", lw=1.5, mutation_scale=13))
    ax1.text(rx + 0.14 * tx, ry + 0.88 * ty, r"$\hat{t}$", fontsize=13,
             color="#1F4E79", ha="center")

    hlen = 1.7
    ax1.plot([rx - hlen * tx, rx + hlen * tx], [ry - hlen * ty, ry + hlen * ty],
             ls="--", color="#bbb", lw=0.9)

    g = np.deg2rad(20)
    vx, vy = np.cos(th + g), np.sin(th + g)
    ax1.annotate("", xy=(rx + 1.05 * vx, ry + 1.05 * vy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#2E6B4E", lw=2.4, mutation_scale=16))
    ax1.text(rx + 1.18 * vx, ry + 1.18 * vy, r"$v$", fontsize=18, color="#2E6B4E")

    p = np.deg2rad(42)
    fx, fy = np.cos(th + p), np.sin(th + p)
    ax1.annotate("", xy=(rx + 1.10 * fx, ry + 1.10 * fy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#b23c3c", lw=2.8, mutation_scale=17))
    ax1.text(rx + 1.26 * fx, ry + 1.26 * fy, r"$T$", fontsize=19, color="#b23c3c")

    ax1.annotate("", xy=(rx - 0.55 * vx, ry - 0.55 * vy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#8a6d3b", lw=2.2, mutation_scale=15))
    ax1.text(rx - 0.70 * vx - 0.06, ry - 0.70 * vy, r"$D$", fontsize=16, color="#8a6d3b")

    ax1.annotate("", xy=(rx - 0.55 * ux, ry - 0.55 * uy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#555", lw=1.8, mutation_scale=14))
    ax1.text(rx - 0.70 * ux, ry - 0.70 * uy, r"$g$", fontsize=14, color="#555")

    arc = np.linspace(0, g, 40)
    ar = 0.42
    ax1.plot(rx + ar * np.cos(th + arc), ry + ar * np.sin(th + arc),
             color="#2E6B4E", lw=1.3)
    ax1.text(rx + 0.58 * np.cos(th + g * 0.5) - 0.06, ry + 0.58 * np.sin(th + g * 0.5) - 0.12,
             r"$\gamma$", fontsize=15, color="#2E6B4E")
    arc2 = np.linspace(0, p, 40)
    ar2 = 0.68
    ax1.plot(rx + ar2 * np.cos(th + arc2), ry + ar2 * np.sin(th + arc2),
             color="#b23c3c", lw=1.3)
    ax1.text(rx + 0.86 * np.cos(th + p * 0.5) - 0.08, ry + 0.86 * np.sin(th + p * 0.5) + 0.04,
             r"$\varphi$", fontsize=15, color="#b23c3c")

    fig.savefig(OUT + r"\fig_forces_coords.pdf")
    plt.close(fig)
    print("fig_forces_coords 完成")


# ===========================================================================
# 图 2：飞行剖面
# ===========================================================================
def fig_flight_profile():
    fig, ax = plt.subplots(figsize=(10.8, 6.4))

    t = np.linspace(0, 1, 600)
    h = np.zeros_like(t)
    tb1, tb2, tb3 = 0.235, 0.485, 0.90
    m1 = t <= 0.035
    h[m1] = 42 * (t[m1] / 0.035) ** 1.2
    m2 = (t > 0.035) & (t <= tb1)
    s = (t[m2] - 0.035) / (tb1 - 0.035)
    h[m2] = 42 + (208 - 42) * (s ** 0.75)
    m3 = (t > tb1) & (t <= tb2)
    s = (t[m3] - tb1) / (tb2 - tb1)
    h[m3] = 208 + (272 - 208) * s ** 0.9
    m4 = t > tb2
    s = (t[m4] - tb2) / (1 - tb2)
    h[m4] = 272 + (400 - 272) * s ** 0.82

    ax.plot(t, h, color="#1F4E79", lw=2.6)
    ax.axhline(400, ls="--", color="#b23c3c", lw=1.4)

    # 阶段分隔竖线 + 区间双箭头标注（直观对应，间距拉开）
    for x in [0.035, tb1, tb2]:
        ax.axvline(x, ls=":", color="#888", lw=1.1)

    def seg_arrow(x0, x1, y, label, sub, color):
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=2.0))
        ax.text((x0 + x1) / 2, y + 12, label, fontsize=FS_TEXT + 1,
                ha="center", color=color, fontweight="bold")
        ax.text((x0 + x1) / 2, y - 18, sub, fontsize=FS_TEXT - 1,
                ha="center", color="#555")

    # 区间箭头逐层抬高，避免相互遮挡；曲线下方留出 t1/t2 标注区
    seg_arrow(0.005, 0.030, 40, "垂直起飞段", "0–10 s", "#111")
    seg_arrow(0.05, 0.225, 92, "程序转弯段", "10 s – $t_1$", "#111")
    seg_arrow(0.245, 0.475, 148, "无动力滑行段", "$t_1$ – $t_2$", "#8a6d3b")
    seg_arrow(0.495, 0.995, 190, "入轨修正段", "$t_2$ – $t_3$", "#b23c3c")

    # 关键时刻标注（t1/t2 放曲线上方，t3 放曲线下方避免与目标轨道文字重叠）
    for x, label, y, ha, dy in [
        (tb1, "$t_1$ 一级关机", 208, "center", 26),
        (tb2, "$t_2$ 二级点火", 272, "center", 26),
        (tb3, "$t_3$ 关机入轨", 400, "center", -30),
    ]:
        ax.plot(x, y, "o", ms=7, color="#1F4E79")
        ax.text(x, y + dy, label, fontsize=FS_TEXT, ha=ha, color="#333")
    ax.text(0.965, 415, "目标轨道 $H=400$ km", ha="right", fontsize=FS_TEXT,
            color="#b23c3c")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 445)
    ax.set_xlabel("归一化时间 $t$", fontsize=FS_LABEL)
    ax.set_ylabel("高度 $h$ / km", fontsize=FS_LABEL)
    ax.grid(alpha=0.3)
    _set_default(ax)
    ax.set_title("飞行剖面与阶段划分", fontsize=FS_TITLE)
    fig.savefig(OUT + r"\fig_flight_profile.pdf")
    plt.close(fig)
    print("fig_flight_profile 完成")


# ===========================================================================
# 图 3：入轨条件
# ===========================================================================
def fig_orbit_condition():
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-2.20, 2.20)
    ax.set_ylim(-1.85, 1.85)

    earth = plt.Circle((0, 0), 1.0, fc="#e8edf3", ec="#7a8ea8", lw=1.4)
    ax.add_patch(earth)
    orbit = plt.Circle((0, 0), 1.27, fill=False, ec="#b23c3c", lw=1.6, ls="--")
    ax.add_patch(orbit)

    th = np.deg2rad(46)
    rx, ry = 1.27 * np.cos(th), 1.27 * np.sin(th)
    ax.plot(rx, ry, "o", ms=9, color="#1F4E79")

    ax.annotate("", xy=(rx, ry), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-", color="#333", lw=1.8))
    ax.text(0.42 * np.cos(th * 0.5) - 0.12, 0.42 * np.sin(th * 0.5) - 0.12,
            r"$r$", fontsize=16, color="#333")

    tv = th - np.deg2rad(70)
    ax.annotate("", xy=(rx + 0.62 * np.cos(tv), ry + 0.62 * np.sin(tv)),
                xytext=(rx, ry),
                arrowprops=dict(arrowstyle="-|>", color="#2E6B4E", lw=2.4, mutation_scale=16))
    ax.text(rx + 0.80 * np.cos(tv) - 0.04, ry + 0.80 * np.sin(tv) + 0.06,
            r"$v$", fontsize=17, color="#2E6B4E")
    ax.text(rx + 0.10, ry + 0.22, r"关机点 $t_3$", fontsize=FS_TEXT,
            color="#1F4E79", ha="left")

    # 入轨条件（纯文字描述——不含完整数学公式）
    ax.text(-2.02, 1.42, "入轨条件", fontsize=FS_TEXT + 2, color="#111",
            ha="left", fontweight="bold")
    items = [
        "轨心距等于目标轨道半径 6771 km",
        "速度大小等于该半径处的圆轨道速度",
        "速度方向沿当地切向，即径向速度为零",
    ]
    for i, it in enumerate(items):
        ax.text(-2.02, 1.14 - 0.26 * i, it, fontsize=FS_TEXT + 1,
                color="#333", ha="left")
    ax.text(-2.02, 0.16, "速度在惯性系中度量，取顺行方向",
            fontsize=FS_TEXT - 1, color="#666", ha="left")

    fig.savefig(OUT + r"\fig_orbit_condition.pdf")
    plt.close(fig)
    print("fig_orbit_condition 完成")


# ===========================================================================
# 图 4：俯仰角律
# ===========================================================================
def fig_pitch_laws():
    fig, ax = plt.subplots(figsize=(10.0, 5.8))
    tau = np.linspace(0, 1, 400)
    phi0 = 26.03
    phi_q1 = phi0 - 25.5 * tau ** 1.35
    phi_q2 = phi0 - 22.03 * tau
    q3_tau = [0.0, 0.25, 0.50, 0.75, 1.0]
    q3_phi = [30.033, 11.142, 11.425, 7.715, 6.189]
    ax.plot(tau, phi_q1, ls="--", color="#999", lw=2.0, label="问题1 恒攻角")
    ax.plot(tau, phi_q2, color="#1F4E79", lw=2.4, label="问题2 恒定速率")
    ax.plot(tau, np.interp(tau, q3_tau, q3_phi), color="#b23c3c", lw=2.6,
            label="问题3 分段线性最优")
    ax.plot(q3_tau, q3_phi, "o", ms=7, color="#b23c3c")
    ax.set_xlabel(r"归一化时间 $\tau=(t-t_2)/(t_3-t_2)$", fontsize=FS_LABEL)
    ax.set_ylabel(r"俯仰角 $\varphi$ / deg", fontsize=FS_LABEL)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 36)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=FS_TEXT, loc="lower left")
    _set_default(ax)
    ax.set_title("三种控制策略下的俯仰角律", fontsize=FS_TITLE)
    ax.annotate(r"$\varphi_0$ 点火时刻速度方向", xy=(0.01, phi0),
                xytext=(0.12, 32.5), fontsize=FS_TEXT, color="#333",
                arrowprops=dict(arrowstyle="->", color="#333", lw=1.0))
    fig.savefig(OUT + r"\fig_pitch_laws.pdf")
    plt.close(fig)
    print("fig_pitch_laws 完成")


# ===========================================================================
# 图 5：问题四求解框架
# ===========================================================================
def fig_q4_framework():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.6))
    for ax in (ax1, ax2):
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    ax1.set_title("优化问题", fontsize=FS_TITLE, fontweight="bold", color="#111")
    lines_l = [
        (r"$\min\ J=m(t_2^{+})-m(t_3)$", "#333", 0.95),
        (r"二子级推进剂消耗", "#666", 0.86),
        (r"s.t.  $|r|=R_E+H$", "#333", 0.74),
        (r"$\qquad |v|=\sqrt{\mu/r}$", "#333", 0.64),
        (r"$\qquad r\cdot v=0$", "#333", 0.54),
        (r"$\qquad \sigma(t)\in[0.6,1.0]$", "#333", 0.44),
        (r"决策变量", "#1F4E79", 0.30),
        (r"$t_c$，滑行时间", "#333", 0.20),
        (r"$\varphi$，5 节点分段线性", "#333", 0.10),
        (r"$\sigma$，4 节点分段线性", "#333", 0.00),
        (r"$t_3$，关机时刻，自由", "#333", -0.10),
    ]
    for s, c, y in lines_l:
        ax1.text(0.05, y, s, fontsize=FS_TEXT, color=c, va="center")

    ax2.set_title("求解流程", fontsize=FS_TITLE, fontweight="bold", color="#b23c3c")
    steps = [
        "① 点火前轨迹按 $t_c$ 缓存",
        "② DOP853 积分二级段",
        "③ 计算目标与入轨残差",
        "④ SLSQP，入轨残差为等式约束",
        "⑤ 88 组初值并行，16 进程",
    ]
    y = 0.92
    for i, s in enumerate(steps):
        ax2.add_patch(plt.Rectangle((0.05, y - 0.085), 0.90, 0.085,
                                    fc="white", ec="#1F4E79" if i < 3 else "#b23c3c",
                                    lw=1.5, transform=ax2.transAxes))
        ax2.text(0.50, y - 0.042, s, fontsize=FS_TEXT, ha="center", va="center")
        if i < 4:
            ax2.annotate("", xy=(0.50, y - 0.085), xytext=(0.50, y - 0.0),
                         arrowprops=dict(arrowstyle="-|>", color="#333", lw=1.3))
        y -= 0.155
    ax2.text(0.50, y + 0.02, "完整 Simulator 复验，残差小于 $10^{-5}$ 才接受",
             fontsize=FS_TEXT - 1, color="#666", ha="center")

    fig.savefig(OUT + r"\fig_q4_framework.pdf")
    plt.close(fig)
    print("fig_q4_framework 完成")


if __name__ == "__main__":
    fig_forces_coords()
    fig_flight_profile()
    fig_orbit_condition()
    fig_pitch_laws()
    fig_q4_framework()
