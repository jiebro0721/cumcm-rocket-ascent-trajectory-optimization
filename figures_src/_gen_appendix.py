# 重新生成论文附录源代码（仅保留关键求解实现，与 src/ 主链一致）
import os

WS = r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2"

FILES = [
    ("src/common.py", "公共模块 common.py：物理常数、动力学方程、分段积分器与入轨残差"),
    ("src/q1_baseline.py", "问题一基准弹道仿真 q1_baseline.py"),
    ("src/q2_inscription.py", "问题二入轨打靶 q2_inscription.py"),
    ("src/q34_direct_collocation.py", "问题三、四直接配点求解器 q34_direct_collocation.py"),
]

out = []
out.append("% 附录：关键源代码（自动生成，与仓库 src/ 完全一致）")
out.append("% 注意：lstlisting 内使用 ttfamily 字体需包含中文，由 xeCJK 处理")
out.append("\\clearpage")

for path, title in FILES:
    full = os.path.join(WS, path)
    with open(full, encoding="utf-8") as f:
        code = f.read().rstrip("\n")
    title_tex = title.replace("_", "\\_")
    out.append("")
    out.append("\\section{%s}" % title_tex)
    out.append("\\begin{lstlisting}[language=Python]")
    out.append(code)
    out.append("\\end{lstlisting}")

with open(os.path.join(WS, "论文", "appendix_codes.tex"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")

print("已生成 appendix_codes.tex, 包含 %d 个关键文件" % len(FILES))
for path, title in FILES:
    print("  -", path)
