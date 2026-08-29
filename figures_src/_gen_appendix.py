# 重新生成论文附录源代码（与当前 src/ 一致）
import os

WS = r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2"

FILES = [
    ("src/common.py", "全程仿真公共模块 common.py"),
    ("src/q1_baseline.py", "问题一基准控制策略全程仿真 q1_baseline.py"),
    ("src/q2_inscription.py", "问题二入轨条件有界三维打靶 q2_inscription.py"),
    ("src/q34_direct_collocation.py", "问题三、四直接配点求解器 q34_direct_collocation.py"),
    ("src/run_all.py", "一键复现脚本 run_all.py"),
    ("cross_validation/code/rocket_trajectory_solution.py",
     "交叉验证独立实现 rocket_trajectory_solution.py"),
]

out = []
out.append("% 附录：完整源代码（自动生成，与仓库 src/ 与 cross_validation/ 完全一致）")
out.append("% 注意：lstlisting 内使用 ttfamily 字体需包含中文，由 xeCJK 处理")
out.append("\\clearpage")

for path, title in FILES:
    full = os.path.join(WS, path)
    with open(full, encoding="utf-8") as f:
        code = f.read().rstrip("\n")
    # section 标题中的下划线必须转义为 \_
    title_tex = title.replace("_", "\\_")
    out.append("")
    out.append("\\section{%s}" % title_tex)
    out.append("\\begin{lstlisting}[language=Python]")
    out.append(code)
    out.append("\\end{lstlisting}")

with open(os.path.join(WS, "论文", "appendix_codes.tex"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

print("已生成 appendix_codes.tex, 包含 %d 个文件" % len(FILES))
for path, title in FILES:
    print("  -", path)
