"""生成问题三/四直接配点求解流程图（论文图：fig_q34_flow.pdf）。"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "论文" / "fig_q34_flow.pdf"

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

fig, ax = plt.subplots(figsize=(7.2, 7.6))
ax.set_xlim(0, 10)
ax.set_ylim(0, 15.6)
ax.axis("off")

# 框定义：(文本, 颜色, y, 高, 是否为宽框)
boxes = [
    ("统一多阶段动力学 + 顺行入轨终端流形\n（问题一的物理底座）", "#f7f9fc", 14.4, 1.0),
    ("决策变量：滑行时间 t_c、燃烧时间 t_b、\n状态节点、俯仰角节点 φ_k（问题四再加节流节点 σ_k）", "#f7f9fc", 12.9, 1.0),
    ("Hermite--Simpson 转录：中点公式 + 缺陷约束\n（动力学化为稀疏非线性规划的显式约束）", "#f7f9fc", 11.4, 1.0),
    ("NLP 求解：CasADi 自动微分 + 内点法\n粗网格 (N_c,N_b)=(10,20) 求解", "#fdf3f3", 9.9, 1.0),
    ("解插值到细网格 (20,40) 重解\n（问题四另加 σ≡0.8 低节流独立初值）", "#fdf3f3", 8.4, 1.0),
    ("网格一致性比较 + 独立 DOP853 重积分复验\n（终端高度、径向/切向速度误差验收）", "#f0f4fa", 6.9, 1.0),
    ("输出：局部最优候选解（耗药、控制律、轨迹）", "#eef6ee", 5.4, 1.0),
]

for i, (text, fc, y, h) in enumerate(boxes):
    box = FancyBboxPatch((1.2, y), 7.6, h, boxstyle="round,pad=0.12,rounding_size=0.22",
                         linewidth=1.3, edgecolor="#1F4E79", facecolor=fc)
    ax.add_patch(box)
    ax.text(5.0, y + h / 2, text, ha="center", va="center", fontsize=10.2, color="#111")
    if i < len(boxes) - 1:
        arr = FancyArrowPatch((5.0, y - 0.08), (5.0, y - h + 0.08 + 0.0),
                              arrowstyle="-|>", mutation_scale=16,
                              linewidth=1.2, color="#333")
        # 下一框顶部在 y_next + h_next
        y_next, h_next = boxes[i + 1][2], boxes[i + 1][3]
        arr = FancyArrowPatch((5.0, y - 0.10), (5.0, y_next + h_next + 0.10),
                              arrowstyle="-|>", mutation_scale=16,
                              linewidth=1.2, color="#333")
        ax.add_patch(arr)

ax.text(5.0, 15.35, "问题三、四直接配点求解流程", ha="center", va="center",
        fontsize=13, fontweight="bold", color="#111")
fig.savefig(OUT, bbox_inches="tight")
print("saved:", OUT)
