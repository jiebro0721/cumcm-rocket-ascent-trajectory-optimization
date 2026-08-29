# 贡献指南

感谢参与本项目！本仓库为建模竞赛项目的完整实现，包含主实现（`src/`）与独立交叉验证实现（`cross_validation/`）。为保证可复现性与结果一致性，请遵循以下约定。

## 项目结构

| 目录 | 职责 | 提交要求 |
|---|---|---|
| `src/` | 主实现（DOP853 + least_squares + CasADi 配点）| 修改必须跑通全部测试 |
| `cross_validation/` | 独立交叉验证实现（RK4 + scipy 配点）| 修改需说明与主实现的对比 |
| `tests/` | 回归测试 + 交叉验证测试 | 新功能需新增对应测试 |
| `results/` | 论文权威数值结果 | 如变更论文必须同步 |
| `论文/` | LaTeX 论文（main.tex）| 编译通过后提交 main.pdf |
| `legacy/` | 历史版本归档（不参与主流程）| 一般不改动 |

## 开发流程

1. **创建分支**：从 `main` 切出 `feature/xxx` 或 `fix/xxx`；
2. **修改代码**：保持与 `common.py` 的物理模型一致（惯性 ECI 笛卡尔、顺行圆轨道终端条件）；
3. **验证**：
   ```bash
   python -m pytest                # 全部测试须通过
   python src/run_all.py           # 复现全部结果（时间较长可选）
   ```
4. **提交 PR**：描述改动内容、验证结果、对论文/结果的影响；
5. **评审**：Maintainer 独立复现后合并。

## 数值一致性约定

- **权威结果**：`results/` 中的 summary 文件（控制节点 17 位有效数字）；
- **交叉验证**：主实现与 `cross_validation/` 关键数值须一致（Q3 耗药差 < 0.05 kg）；
- **论文数字**：所有表格数据必须能从 `results/` 或 `cross_validation/results/` 追溯，禁止手工录入。

## 论文编译

```bash
cd 论文 && xelatex main.tex && xelatex main.tex   # 需 XeLaTeX
```

如无法编译，请勿提交 `main.pdf`，在 PR 中注明。
