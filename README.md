# 两级运载火箭发射轨迹建模与入轨优化

全国大学生数学建模竞赛培训项目：两级运载火箭（赤道发射 → 400 km 近地圆轨道）的**动力学建模、基准仿真、入轨打靶与燃料最优控制**完整实现。本仓库包含全部可运行源码、数值结果、插图与论文正文，可端到端复现。

---

## 项目结构

```
.
├── src/                    # 全部源代码（权威实现）
│   ├── common.py           # 物理常数、火箭参数、动力学方程、分段积分器、顺行入轨残差
│   ├── q1_baseline.py      # 问题1：基准控制策略全程仿真
│   ├── q2_inscription.py   # 问题2：入轨条件有界三维打靶（least_squares + 196 组多初值）
│   ├── q34_direct_collocation.py # 问题3/4：Hermite-Simpson直接配点 + CasADi/Ipopt
│   └── run_all.py          # 一键复现：Q1→Q4→图表
├── legacy/shooting_baseline/  # 旧版单重打靶求解器与结果（仅方法对比，非权威）
├── 论文/                   # 参赛论文（LaTeX，cumcmthesis 模板）
│   ├── main.tex            # 论文正文（四问模型、算法、结果、评价）
│   └── main.pdf            # 编译产物（24 页）
├── figures/                # 论文插图（PNG，由 src/ 生成）
├── figures_src/            # 示意图源（SVG 矢量 + 生成脚本）
├── results/                # 权威数值结果（CSV：summary/trajectory/convergence）
├── tests/                  # 回归测试（pytest，8 项）
├── 文献/                   # 参考文献（中外期刊论文 + 开源求解器代码）
└── requirements.txt        # Python 依赖
```

## 运行环境与依赖

- Python ≥ 3.10
- `pip install -r requirements.txt`（numpy、pandas、scipy、matplotlib、casadi）
- 论文编译：XeLaTeX（MiKTeX/TeXLive），`xelatex main.tex` 两遍

## 复现步骤（单一入口）

```bash
# 一键复现论文全部数值结果与图表
python src/run_all.py

# 或分步执行
python src/q1_baseline.py      # 问题1：基准策略仿真
python src/q2_inscription.py   # 问题2：入轨条件有界打靶
python src/q34_direct_collocation.py --question 3 --n-coast 10 --n-burn 20   # Q3 粗网格
python src/q34_direct_collocation.py --question 3 --n-coast 20 --n-burn 40   # Q3 细网格
python src/q34_direct_collocation.py --question 4 --n-coast 20 --n-burn 40   # Q4（两类初值内置）

# 运行回归测试
python -m pytest

# 编译论文（在 论文/ 目录）
cd 论文 && xelatex main.tex && xelatex main.tex
```

**结果链约定**：`results/` 是论文唯一数据源（Q1--Q4 的 summary 与控制节点以 17 位有效数字保存）；`legacy/` 中的旧版单重打靶结果仅供方法对比，不作为论文数字来源。

## 核心建模

**动力学模型**（惯性 ECI 平面分量，含地球自转）：

- 状态 `x = [x, y, vx, vy, m]`，赤道平面内
- 初值 `v0 = ω_E·R_E`（赤道东射约 465 m/s 惯性初速）
- 重力 `−μr/r³`、推力 `(T/m)·û`（俯仰角控制）、阻力用相对速度 `D = −½ρSC_D|v_rel|v_rel`
- 指数大气 `ρ(h) = 1.225·e^(−h/7200)`
- 飞行剖面：垂直起飞 10 s → 程序转弯 0.4°/s → 一子级分离 → 无动力滑行 → 二子级点火入轨
- 顺行入轨条件：`|r| = R_E+400 km`，`v_r = 0`，`v_t = +√(μ/|r|)`（残差以切向速度符号严格排除逆行分支）

**求解算法**：

| 问题 | 方法 |
| :--- | :--- |
| 问题1 | 分段 DOP853 自适应积分 + 事件检测（推进剂耗尽） |
| 问题2 | 有界三维打靶（`least_squares`，变量 [t_c, k, t_b]）+ 完整动力学重积分复验 |
| 问题3 | Hermite-Simpson 直接配点 + 自动微分 + Ipopt；DOP853 独立复验 |
| 问题4 | 与问题3共用配点模型，增加节流变量和累计推进会约束 |

问题3/4使用问题2轨迹暖启动和网格加密（10/20 → 20/40），不再依赖几十组盲目初值。所得候选解均用完整 `Simulator` 端到端复验；非凸问题只声明局部最优，不作全局保证。

## 主要结果

| 问题 | 结果 |
| :--- | :--- |
| 问题1 | 二子级燃尽：h = 433.53 km、V = 8521.18 m/s、γ = −0.69°（进入椭圆轨道而非目标圆轨道） |
| 问题2 | t_c = 103.63 s、k = −0.0474°/s、t_b = 398.05 s（关机 650.74 s），剩余推进剂 4015 kg |
| 问题3 | 20/40网格：t_c = 4.3767 s、燃烧 394.3534 s，消耗 57446.934 kg |
| 问题4 | 20/40网格：σ ≈ 1，消耗 57446.935 kg；当前模型下未发现节流收益 |

## 参考文献

1. Nair V S, Vaidyanathan A. Ascent trajectory design and optimization of a two-stage throttleable liquid rocket. *Advances in Space Research*, 2022.
2. Benedikter B, et al. Convex optimization of launch vehicle ascent trajectory with heat-flux and splash-down constraints. *Journal of Spacecraft and Rockets*, 2022.
3. Betts J T. Survey of numerical methods for trajectory optimization. *Journal of Guidance, Control, and Dynamics*, 1998.
4. 刘超越, 张成. 基于高斯伪谱法的二级助推战术火箭多阶段轨迹优化. 兵工学报, 2019.
5. 胡冬生 等. 含滑行时间约束的真空段弹道设计研究. 宇航总体技术, 2023.
6. 李惠峰 等. 液体火箭上升段制导方法的发展综述. 航天控制, 2023.
7. 余梦伦 等. 弹道学在运载火箭总体设计中的实践与展望. 宇航总体技术, 2023.

详见 `论文/main.tex` 参考文献节与 `文献/` 目录。
