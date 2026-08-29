# 独立交叉验证实现

本目录是项目的**第二套完全独立实现**，与 `src/` 主实现并列，用于交叉验证数值结果（贡献来自队友 wyy-w11）。

## 实现特点（与主实现完全独立）

- **数值方法**：固定步长经典四阶 Runge-Kutta（非 DOP853）；Q2 用终端灵敏度矩阵驱动的 Newton 预测-校正 + 回溯线搜索（非 least_squares）；Q3/Q4 用 Hermite-Simpson 直接配点 + scipy 优化器（非 CasADi/Ipopt）。
- **坐标系**：极坐标/径向-切向分量（主实现为惯性笛卡尔），物理模型等价。

## 验证结论（与主实现一致）

| 项目 | 独立实现 | 主实现（src/q34_direct_collocation.py） | 一致性 |
|---|---|---|---|
| Q1 燃尽 | h=433.527 km, V=8521.18 m/s, γ=−0.686° | 433.53 km / 8521.18 / −0.69° | ✅ |
| Q2 | t_c=103.632, k=−0.0474, t_b=398.049, 耗药 57985.3 | 完全相同 | ✅ |
| Q3 最优 | t_c=4.374, 燃烧 394.353 s, 耗药 57446.933 kg | 57446.934 kg | ✅（差 0.002 kg）|
| Q4 最优 | σ 节点≈1.0, 耗药 57446.933 kg | 57446.935 kg | ✅ |

## 独特贡献（已并入论文 §7.3.2/§7.4.3-7.4.6）

1. **三档网格收敛**（11/21/31 节点，表 tab:q3mesh）：证明结果与网格加密路径无关；
2. **五档节流初值审计**（表 tab:q4audit）：0.6~1.0 五档初值全部收敛到近似满推力解；
3. **恒定节流灵敏度定量结果**（表 tab:q4sens）：σ=0.9/0.8/0.7/0.6 分别增加耗药 47.1/196.2/532.3/1279.4 kg；
4. **速度损失预算**（表 tab:loss）：Q2 重力损失 367.1、转向损失 53.3 m/s；Q3 为 700.0/78.0 m/s——定量解释滑行时间与重力损失的权衡。

## 数据文件

- `code/rocket_trajectory_solution.py`：独立求解器（840 行，可独立运行）
- `results/results.json`：全部结果（Q1-Q4 + 常数/方法说明）
- `results/q3_mesh_convergence.csv`：三档网格收敛（11/21/31 节点）
- `results/q4_constant_throttle_sensitivity.csv`：恒节流灵敏度（σ=0.6~1.0）
- `results/q4_multistart_audit.csv`：多初值审计（五档节流初值）
- `results/second_stage_loss_budget.csv`：速度损失预算（Q2/Q3 重力+转向损失）
- `results/rk4_step_convergence.csv`：RK4 步长收敛性（证明独立传播器精度）
- `results/q2_newton_iterations.csv`：Q2 Newton 预测-校正迭代过程
- `results/q{n}_rk4_trajectory.csv`：各问题 RK4 完整轨迹（与主实现笛卡尔轨迹可比）
- `results/intermediate_phase_states.csv`：各阶段衔接点状态
- `results/figures/*.png`：独立实现的轨迹/控制律图（optimized_controls、q1_baseline_history、trajectory_comparison）
