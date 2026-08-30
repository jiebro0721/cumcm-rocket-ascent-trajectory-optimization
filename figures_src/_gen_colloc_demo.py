"""生成论文图7：直接配点法示意图（简化图名 + 近似与实际曲线差异可见）。

设计：
  - 真实运动轨迹（细黑线）：解析函数（含拐点）
  - 配点近似（彩色分段曲线）：每段二次插值（仅用端点+中点三个拟合点），
    与真实曲线存在可见误差——直观展示"在节点与中点匹配"的配点本质
  - 标注：knot point（空心圈）、mid-point（实心点）、segment（下方括号）
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

OUT = r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2\论文"


def true(t):
    """真实轨迹：含拐点的平滑函数（示意状态）。"""
    return 250 + 150 * (0.5 + 0.5 * np.sin(2.0 * np.pi * t - 0.6)) ** 1.4


def quad_piece(t0, t1, x0, x1, xm):
    """二次插值，过 (t0,x0) (中, xm) (t1,x1)。"""
    ts = np.linspace(t0, t1, 50)
    tm = 0.5 * (t0 + t1)
    # 拉格朗日三点插值
    def L(t, pts, vals):
        out = np.zeros_like(t)
        for i, (pt, val) in enumerate(zip(pts, vals)):
            num = np.ones_like(t)
            den = 1.0
            for j, (ptj, _) in enumerate(zip(pts, vals)):
                if j != i:
                    num = num * (t - ptj)
                    den = den * (pt - ptj)
            out = out + val * num / den
        return out

    pts = [t0, tm, t1]
    vals = [x0, xm, x1]
    return ts, L(ts, pts, vals)


fig, ax = plt.subplots(figsize=(9.6, 5.2))

# 真实轨迹
N = 1001
tau = np.linspace(0.0, 1.0, N)
ax.plot(tau, true(tau), color="#333", lw=1.5, label="真实运动轨迹")

# 配点近似：3 段二次插值
nseg = 3
tk = np.linspace(0.0, 1.0, nseg + 1)
colors = ["#1F4E79", "#b23c3c", "#2E6B4E"]
for i in range(nseg):
    t0, t1 = tk[i], tk[i + 1]
    tm = 0.5 * (t0 + t1)
    x0, xm, x1 = true(t0), true(tm), true(t1)
    ts, xs = quad_piece(t0, t1, x0, x1, xm)
    ax.plot(ts, xs, lw=3.0, color=colors[i],
            label="配点近似（三个拟合点）" if i == 0 else None)

# 节点（空心圈）与中点（实心点）
for t in tk:
    ax.plot(t, true(t), "o", ms=13, mfc="white", mec="#111", mew=2.0, zorder=5)
for t in (tk[:-1] + tk[1:]) / 2:
    ax.plot(t, true(t), "o", ms=6, color="#111", zorder=5)

# segment 括号
ybox = 190
for i in range(nseg):
    ax.plot([tk[i] + 0.03, tk[i] + 0.03, tk[i + 1] - 0.03, tk[i + 1] - 0.03],
            [ybox, ybox - 10, ybox - 10, ybox],
            color="#555", lw=1.2)
    ax.text((tk[i] + tk[i + 1]) / 2, ybox - 22, r"segment $[t_k, t_{k+1}]$",
            ha="center", fontsize=10, color="#555")

# 标注
ax.annotate("knot point 配点节点", xy=(tk[0], true(tk[0])),
            xytext=(0.05, 350), fontsize=11, color="#111",
            arrowprops=dict(arrowstyle="->", color="#555", lw=1.2))
ax.annotate("mid-point 中点", xy=((tk[0] + tk[1]) / 2 + 0.01, true((tk[0] + tk[1]) / 2)),
            xytext=(0.22, 168), fontsize=11, color="#111",
            arrowprops=dict(arrowstyle="->", color="#555", lw=1.2))

ax.set_xlabel(r"燃烧段归一化时间 $\tau=(t-t_2)/t_b$", fontsize=12)
ax.set_ylabel(r"状态 $x$", fontsize=12)
ax.set_ylim(140, 420)
ax.grid(alpha=0.3)
ax.legend(fontsize=10, loc="upper center", framealpha=0.9)

fig.tight_layout()
fig.savefig(OUT + r"\fig_q34_flow.pdf")
plt.close(fig)
print("fig_q34_flow.pdf 已重新生成（近似误差可见）")
