"""单一公开入口：完整复现论文全部数值结果与图表。

用法（仓库根目录）：
    python src/run_all.py

执行内容：
  Q1  基准策略仿真 -> results/q1_*.csv, figures/q1_*.png
  Q2  有界三维打靶 -> results/q2_*.csv
  Q3  配点 10/20 与 20/40（细网格解插值重解）-> results/q3_collocation_summary.csv
  Q4  配点 20/40（满推力与低节流两类初值）-> results/q4_collocation_summary.csv
  图表 fig_q3_trajectory/controls, fig_q4_trajectory/controls
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print(">>>", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd or ROOT, check=True)


def main() -> None:
    py = sys.executable

    # Q1
    run([py, "src/q1_baseline.py"])
    # Q2
    run([py, "src/q2_inscription.py"])
    # Q3 粗网格 -> 细网格
    run([py, "src/q34_direct_collocation.py",
         "--question", "3", "--n-coast", "10", "--n-burn", "20"])
    run([py, "src/q34_direct_collocation.py",
         "--question", "3", "--n-coast", "20", "--n-burn", "40"])
    # Q4 细网格（两类初值由求解器内部完成）
    run([py, "src/q34_direct_collocation.py",
         "--question", "4", "--n-coast", "20", "--n-burn", "40"])
    # 论文图表
    run([py, str(ROOT / "figures_src" / "_gen_q34_plots.py")])

    print("\n全部复现完成。数值见 results/，图表见 论文/ 与 figures/。")


if __name__ == "__main__":
    main()
