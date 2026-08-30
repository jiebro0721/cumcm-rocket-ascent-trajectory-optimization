"""生成论文图6：Hermite--Simpson 直接配点法示意图（贴合本项目）。

仿照经典配点法示意：
  - 真实轨迹（细灰线 = 由高精度积分器得到）
  - 配点近似（分段彩色曲线 = 每段上的 Hermite 三次插值）
  - 标注：knot point（配点节点，空心圈）、mid-point（中点，实心点）、
    segment（分段区间，下方括号）
  - 横轴为燃烧段归一化时间 tau，纵轴为示意状态（取火箭高度 h 量级）
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

OUT = r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2\论文"

# ---- 真实轨迹（模拟 Q3 最优解的高度曲线量级）----
N = 1001
tau = np.linspace(0.0, 1.0, N)
h_true = 250 + 160 * np.sin(np.pi * tau * 0.78) ** 1.15  # 示意高度 / km


def hermite_segment(t0, t1, x0, x1, f0, f1):
    """三次 Hermite 插值端点值+导数。"""
    ts = np.linspace(t0, t1, 60)
    s = (ts - t0) / (t1 - t0)
    h00 = 2 * s**3 - 3 * s**2 + 1
    h10 = s**3 - 2 * s**2 + s
    h01 = -2 * s**3 + 3 * s**2
    h11 = s**3 - s**2
    dt = t1 - t0
    return ts, h00 * x0 + h10 * dt * f0 + h01 * x1 + h11 * dt * f1


def true(t):
    return 250 + 160 * np.sin(np.pi * t * 0.78) ** 1.15


def dtrue(t):
    h = 1e-6
    return (true(t + h) - true(t - h)) / (2 * h)


fig, ax = plt.subplots(figsize=(9.6, 5.4))

# 真实轨迹
ax.plot(tau, h_true, color="#333", lw=1.3, label="真实运动轨迹")

# 配点近似：4 个 segment，每个用 Hermite 三次插值（端点值+端点导数）
nseg = 4
tk = np.linspace(0.0, 1.0, nseg + 1)
colors = ["#1F4E79", "#b23c3c", "#8a6d3b", "#2E6B4E"]
for i in range(nseg):
    t0, t1 = tk[i], tk[i + 1]
    x0, x1 = true(t0), true(t1)
    f0, f1 = dtrue(t0), dtrue(t1)
    ts, xs = hermite_segment(t0, t1, x0, x1, f0, f1)
    ax.plot(ts, xs, lw=2.6, color=colors[i])

# 节点 knot point（空心圈）
for t in tk:
    ax.plot(t, true(t), "o", ms=13, mfc="white", mec="#111", mew=2.0, zorder=5)
# 中点 mid-point（实心小点）
for t in (tk[:-1] + tk[1:]) / 2:
    ax.plot(t, true(t), "o", ms=6, color="#111", zorder=5)

# 分段括号（segment）
ybox = 215
for i in range(nseg):
    ax.plot([tk[i] + 0.012, tk[i] + 0.012, tk[i + 1] - 0.012, tk[i + 1] - 0.012],
            [ybox, ybox - 10, ybox - 10, ybox],
            color="#555", lw=1.2)
    ax.text((tk[i] + tk[i + 1]) / 2, ybox - 22, f"segment $[t_k, t_{{k+1}}]$",
            ha="center", fontsize=10, color="#555")

# 文字标注（箭头）
ax.annotate("knot point 配点节点", xy=(tk[0], true(tk[0])),
            xytext=(0.05, 340), fontsize=11, color="#111",
            arrowprops=dict(arrowstyle="->", color="#555", lw=1.2))
ax.annotate("mid-point 中点", xy=((tk[0] + tk[1]) / 2 + 0.01, true((tk[0] + tk[1]) / 2)),
            xytext=(0.22, 165), fontsize=11, color="#111",
            arrowprops=dict(arrowstyle="->", color="#555", lw=1.2))
ax.annotate("Hermite--Simpson 近似\n（在节点与中点匹配动力学）",
            xy=(0.62, true(0.62)), xytext=(0.55, 355), fontsize=11,
            color="#1F4E79",
            arrowprops=dict(arrowstyle="->", color="#1F4E79", lw=1.4))

# 坐标轴
ax.set_xlabel(r"燃烧段归一化时间 $\tau=(t-t_2)/t_b$", fontsize=12)
ax.set_ylabel(r"状态 $x$", fontsize=12)
ax.set_ylim(150, 420)
ax.grid(alpha=0.3)
ax.legend(fontsize=10, loc="lower right", framealpha=0.9)
ax.set_title("直接配点法：节点、中点与分段上的动力学匹配", fontsize=14)

fig.tight_layout()
fig.savefig(OUT + r"\fig_q34_flow.pdf")
plt.close(fig)
print("fig_q34_flow.pdf 已替换为配点法示意图")
