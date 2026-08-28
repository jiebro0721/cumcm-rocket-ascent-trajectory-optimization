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
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 6.2))

    # ---- (a) 受力图 ----
    ax1.set_xlim(-1.7, 1.95)
    ax1.set_ylim(-1.85, 1.80)
    ax1.set_aspect("equal")
    ax1.axis("off")
    ax1.set_title("a 横向坐标系与受力", fontsize=FS_TITLE)

    earth = plt.Circle((0, 0), 1.0, fc="#e8edf3", ec="#7a8ea8", lw=1.4)
    ax1.add_patch(earth)
    ax1.plot(0, 0, "o", ms=5, color="#555")
    ax1.text(0.04, -0.08, r"$O$", fontsize=14, ha="left", va="top")
    ax1.text(0.0, -1.28, "地球", ha="center", fontsize=FS_TEXT, color="#555")

    th = np.deg2rad(38)
    rx, ry = 1.30 * np.cos(th), 1.30 * np.sin(th)
    ax1.plot(rx, ry, "o", ms=9, color="#1F4E79")
    ax1.text(rx - 0.08, ry - 0.26, "火箭", ha="center", fontsize=FS_TEXT, color="#1F4E79")

    ax1.annotate("", xy=(rx, ry), xytext=(0, 0),
                 arrowprops=dict(arrowstyle="-", color="#333", lw=1.8))
    ax1.text(0.42 * np.cos(th * 0.55), 0.42 * np.sin(th * 0.55),
             r"$r$", fontsize=15, color="#333")

    ux, uy = np.cos(th), np.sin(th)
    tx, ty = -np.sin(th), np.cos(th)
    ax1.annotate("", xy=(rx + 0.55 * ux, ry + 0.55 * uy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#1F4E79", lw=1.5, mutation_scale=13))
    ax1.text(rx + 0.66 * ux, ry + 0.66 * uy, r"$\hat{r}$", fontsize=14, color="#1F4E79")
    ax1.annotate("", xy=(rx + 0.55 * tx, ry + 0.55 * ty), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#1F4E79", lw=1.5, mutation_scale=13))
    ax1.text(rx + 0.12 * tx, ry + 0.70 * ty, r"$\hat{t}$", fontsize=14,
             color="#1F4E79", ha="center")

    hlen = 1.3
    ax1.plot([rx - hlen * tx, rx + hlen * tx], [ry - hlen * ty, ry + hlen * ty],
             ls="--", color="#999", lw=1.0)
    hx, hy = rx - hlen * tx, ry - hlen * ty
    ax1.text(hx + 0.30 * tx, hy + 0.30 * ty, "当地水平面", fontsize=12.5,
             color="#999", ha="center")

    g = np.deg2rad(22)
    vx, vy = np.cos(th + g), np.sin(th + g)
    ax1.annotate("", xy=(rx + 0.85 * vx, ry + 0.85 * vy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#2E6B4E", lw=2.2, mutation_scale=15))
    ax1.text(rx + 0.92 * vx, ry + 0.92 * vy, r"$v$", fontsize=16, color="#2E6B4E")

    p = np.deg2rad(35)
    fx, fy = np.cos(th + p), np.sin(th + p)
    ax1.annotate("", xy=(rx + 0.88 * fx, ry + 0.88 * fy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#b23c3c", lw=2.5, mutation_scale=16))
    ax1.text(rx + 0.98 * fx, ry + 0.98 * fy, r"$T$", fontsize=17, color="#b23c3c")

    ax1.annotate("", xy=(rx - 0.45 * vx, ry - 0.45 * vy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#8a6d3b", lw=2.0, mutation_scale=14))
    ax1.text(rx - 0.56 * vx - 0.06, ry - 0.56 * vy, r"$D$", fontsize=14, color="#8a6d3b")

    ax1.annotate("", xy=(rx - 0.42 * ux, ry - 0.42 * uy), xytext=(rx, ry),
                 arrowprops=dict(arrowstyle="-|>", color="#555", lw=1.6, mutation_scale=13))
    ax1.text(rx - 0.54 * ux, ry - 0.54 * uy, r"$g$", fontsize=13, color="#555")

    arc = np.linspace(0, g, 40)
    ar = 0.32
    ax1.plot(rx + ar * np.cos(th + arc), ry + ar * np.sin(th + arc),
             color="#2E6B4E", lw=1.2)
    ax1.text(rx + 0.44 * np.cos(th + g * 0.5) - 0.10, ry + 0.44 * np.sin(th + g * 0.5) - 0.11,
             r"$\gamma$", fontsize=14, color="#2E6B4E")
    arc2 = np.linspace(0, p, 40)
    ar2 = 0.52
    ax1.plot(rx + ar2 * np.cos(th + arc2), ry + ar2 * np.sin(th + arc2),
             color="#b23c3c", lw=1.2)
    ax1.text(rx + 0.64 * np.cos(th + p * 0.5) - 0.11, ry + 0.64 * np.sin(th + p * 0.5),
             r"$\varphi$", fontsize=14, color="#b23c3c")

    ax1.text(0.0, -1.62, r"$\gamma$：速度与水平面夹角，飞行路径角", fontsize=12.5,
             color="#555", ha="center")
    ax1.text(0.0, -1.76, r"$\varphi$：推力与水平面夹角，俯仰角", fontsize=12.5,
             color="#555", ha="center")

    # ---- (b) 三阶段控制输入 ----
    ax2.axis("off")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_title("b 三个阶段的控制输入", fontsize=FS_TITLE)
    boxes = [
        ("阶段Ⅰ 垂直起飞与程序转弯", "#1F4E79",
         [r"$T=F_{\max 1}$ 恒定",
          r"$\varphi(t)=90^\circ-0.4^\circ/s\cdot(t-10\,\mathrm{s})$",
          r"终止条件：一子级推进剂耗尽 $t_1$"]),
        ("阶段Ⅱ 无动力滑行", "#8a6d3b",
         [r"$T=0$，重力转向，无阻力",
          r"终止条件：滑行时间 $t_c$ 到期"]),
        ("阶段Ⅲ 二子级入轨修正", "#b23c3c",
         [r"问题1：$\varphi=\gamma$，$T=F_{\max 2}$",
          r"问题2：$\varphi=\varphi_0+k\cdot(t-t_2)$",
          r"问题4：$T=\sigma(t)\,F_{\max 2}$，$\sigma\in[0.6,1]$"]),
    ]
    ypos = [0.885, 0.565, 0.18]
    height = 0.25
    for (title, color, lines), y in zip(boxes, ypos):
        ax2.add_patch(plt.Rectangle((0.03, y - height), 0.94, height,
                                    fc="#f7f9fc", ec="#ccc", lw=1.0))
        ax2.text(0.5, y - 0.03, title, fontsize=FS_BOX + 1.5, ha="center",
                 color=color, fontweight="bold")
        ax2.plot([0.10, 0.90], [y - 0.085, y - 0.085], color=color, lw=1.5)
        for i, line in enumerate(lines):
            ax2.text(0.10, y - 0.145 - 0.062 * i, line, fontsize=FS_BOX,
                     va="center", color="#333")

    fig.savefig(OUT + r"\fig_forces_coords.pdf")
    plt.close(fig)
    print("fig_forces_coords 完成")


# ===========================================================================
# 图 2：飞行剖面
# ===========================================================================
def fig_flight_profile():
    fig, ax = plt.subplots(figsize=(10.5, 6.0))

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

    ax.plot(t, h, color="#1F4E79", lw=2.4)
    ax.axhline(400, ls="--", color="#b23c3c", lw=1.4)
    ax.text(0.93, 418, "目标轨道 $H=400$ km", ha="right", fontsize=FS_TEXT, color="#b23c3c")

    for x, label, y, ha, dy in [
        (tb1, r"$t_1$ 一级关机", 208, "right", 14),
        (tb2, r"$t_2$ 二级点火", 272, "right", 14),
        (tb3, r"$t_3$ 关机入轨", 400, "left", 14),
    ]:
        ax.axvline(x, ls=":", color="#999", lw=1.0)
        ax.text(x - 0.015 if ha == "right" else x + 0.015, y + dy, label,
                fontsize=FS_TEXT, ha=ha, color="#333")

    ax.text(0.19, 68, "垂直起飞段", fontsize=FS_TEXT, ha="center", color="#111")
    ax.text(0.19, 38, "0–10 s", fontsize=FS_TEXT - 1, ha="center", color="#666")
    ax.text(0.39, 108, "程序转弯段", fontsize=FS_TEXT, ha="center", color="#111")
    ax.text(0.39, 78, "10 s – $t_1$", fontsize=FS_TEXT - 1, ha="center", color="#666")
    ax.text(0.66, 243, "无动力滑行段", fontsize=FS_TEXT, ha="center", color="#111")
    ax.text(0.66, 213, "$t_1$ – $t_2$", fontsize=FS_TEXT - 1, ha="center", color="#666")
    ax.text(0.78, 335, "入轨修正段", fontsize=FS_TEXT, ha="center", color="#b23c3c")
    ax.text(0.78, 305, "$t_2$ – $t_3$", fontsize=FS_TEXT - 1, ha="center", color="#666")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 470)
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
    fig, ax = plt.subplots(figsize=(10.8, 6.4))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-2.55, 2.45)
    ax.set_ylim(-1.80, 1.90)

    earth = plt.Circle((0, 0), 1.0, fc="#e8edf3", ec="#7a8ea8", lw=1.4)
    ax.add_patch(earth)
    orbit = plt.Circle((0, 0), 1.27, fill=False, ec="#b23c3c", lw=1.6, ls="--")
    ax.add_patch(orbit)
    ax.text(-2.30, 1.22, "目标圆轨道", fontsize=FS_TEXT, color="#b23c3c", ha="left")
    ax.text(-2.30, 1.06, r"$|r|=R_E+H$", fontsize=FS_TEXT - 1, color="#b23c3c", ha="left")

    th = np.deg2rad(46)
    rx, ry = 1.27 * np.cos(th), 1.27 * np.sin(th)
    ax.plot(rx, ry, "o", ms=9, color="#1F4E79")
    ax.text(rx + 0.08, ry + 0.18, "关机点 $t_3$", fontsize=FS_TEXT, color="#1F4E79",
            ha="left")

    ax.annotate("", xy=(rx, ry), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-", color="#333", lw=1.8))
    ax.text(0.36 * np.cos(th * 0.5) - 0.11, 0.36 * np.sin(th * 0.5) - 0.11,
            r"$r$", fontsize=15, color="#333")

    tv = th - np.deg2rad(70)
    ax.annotate("", xy=(rx + 0.62 * np.cos(tv), ry + 0.62 * np.sin(tv)),
                xytext=(rx, ry),
                arrowprops=dict(arrowstyle="-|>", color="#2E6B4E", lw=2.2, mutation_scale=15))
    ax.text(rx + 0.78 * np.cos(tv) - 0.05, ry + 0.78 * np.sin(tv),
            r"$v$", fontsize=16, color="#2E6B4E")

    ax.text(-2.30, 0.62, "入轨条件", fontsize=FS_TEXT + 1, color="#111",
            ha="left", fontweight="bold")
    items = [
        r"$|r(t_3)|=R_E+H=6771\ \mathrm{km}$",
        r"$|v(t_3)|=\sqrt{\mu/r}=7672.6\ \mathrm{m/s}$",
        r"$r(t_3)\cdot v(t_3)=0$",
    ]
    for i, it in enumerate(items):
        ax.text(-2.30, 0.34 - 0.20 * i, it, fontsize=FS_TEXT, color="#333", ha="left")
    ax.text(-2.30, -0.42, r"第三条等价于 $\gamma=0$", fontsize=FS_TEXT - 1,
            color="#666", ha="left")
    ax.text(-2.30, -0.60, r"速度在惯性系中度量", fontsize=FS_TEXT - 1,
            color="#666", ha="left")

    ax.text(-2.30, -1.52, r"地球自转初速 $v_0=\omega_E R_E\approx 464.6$ m/s，东向",
            fontsize=FS_TEXT - 1, color="#555", ha="left")
    ax.text(-2.30, -1.70, r"阻力按相对速度 $v_{rel}=v-\omega_E\times r$ 计算",
            fontsize=FS_TEXT - 1, color="#555", ha="left")
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
