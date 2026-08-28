# 两级运载火箭发射轨迹建模与入轨优化

全国大学生数学建模竞赛培训项目：两级运载火箭（赤道发射 → 400 km 近地圆轨道）的**动力学建模、基准仿真、入轨打靶与燃料最优控制**完整实现。本仓库包含全部可运行源码、数值结果、插图与论文正文，可端到端复现。

---

## 项目结构

```
.
├── src/                    # 全部源代码
│   ├── common.py           # 物理常数、火箭参数、动力学方程、分段积分器、入轨残差
│   ├── q1_baseline.py      # 问题1：基准控制策略全程仿真
│   ├── q2_inscription.py   # 问题2：入轨条件三维打靶（fsolve + 196 组多初值）
│   ├── q3_fuel_opt.py      # 问题3：燃料最省（控制参数化 + SLSQP + 45 初值并行）
│   └── q4_throttle_opt.py  # 问题4：推力节流下燃料最省（88 初值并行）
├── 论文/                   # 参赛论文（LaTeX，cumcmthesis 模板）
│   ├── main.tex            # 论文正文（四问模型、算法、结果、评价）
│   ├── appendix_codes.tex  # 附录：完整源代码（自动生成，与 src/ 一致）
│   └── main.pdf            # 编译产物（正文 + 附录完整代码）
├── figures/                # 论文插图（PNG，由 src/ 生成）
├── figures_src/            # 示意图源（SVG 矢量 + 生成脚本）
├── results/                # 数值结果（CSV）
├── 文献/                   # 参考文献（中外期刊论文 + 开源求解器代码）
├── 项目审阅报告.md          # 第三方审阅报告及整改记录
└── requirements.txt        # Python 依赖
```

## 运行环境与依赖

- Python ≥ 3.10
- `pip install -r requirements.txt`（numpy、pandas、scipy、matplotlib）
- 论文编译：XeLaTeX（MiKTeX/TeXLive），`xelatex main.tex` 两遍

## 复现步骤

```bash
# 1. 依次运行四问求解器
python src/q1_baseline.py      # 问题1：基准策略仿真
python src/q2_inscription.py   # 问题2：入轨条件打靶
python src/q3_fuel_opt.py      # 问题3：燃料最省优化
python src/q4_throttle_opt.py  # 问题4：推力节流优化

# 2. 结果输出
#    results/q{1,2,3,4}_*.csv  数值结果
#    figures/q*.png            插图

# 3. 编译论文（在 论文/ 目录）
cd 论文 && xelatex main.tex && xelatex main.tex
```

## 核心建模

**动力学模型**（惯性 ECI 平面分量，含地球自转）：

- 状态 `x = [x, y, vx, vy, m]`，赤道平面内
- 初值 `v0 = ω_E·R_E`（赤道东射约 465 m/s 惯性初速）
- 重力 `−μr/r³`、推力 `(T/m)·û`（俯仰角控制）、阻力用相对速度 `D = −½ρSC_D|v_rel|v_rel`
- 指数大气 `ρ(h) = 1.225·e^(−h/7200)`
- 飞行剖面：垂直起飞 10 s → 程序转弯 0.4°/s → 一子级分离 → 无动力滑行 → 二子级点火入轨
- 入轨条件：`|r| = R_E+400 km`，`|v| = √(μ/|r|)`，`r·v = 0`

**求解算法**：

| 问题 | 方法 |
| :--- | :--- |
| 问题1 | 分段 DOP853 自适应积分 + 事件检测（推进剂耗尽） |
| 问题2 | 三维打靶：Newton 型 fsolve + 196 组网格多初值 |
| 问题3 | 控制参数化（分段线性俯仰角）+ SLSQP + 45 组初值并行 |
| 问题4 | 控制参数化（俯仰角 + 节流比）+ SLSQP + 88 组初值并行 |

所有最优解均用完整 Simulator 端到端复验，入轨残差 < 10⁻⁵ 才接受。

## 主要结果

| 问题 | 结果 |
| :--- | :--- |
| 问题1 | 二子级燃尽：h = 433.53 km、V = 8521.18 m/s、γ = −0.69°（过冲未入轨） |
| 问题2 | t_c = 103.63 s、k = −0.0474°/s、t_shut = 650.74 s，剩余推进剂 4015 kg |
| 问题3 | t_c = 103.63 s、t_shut = 650.74 s，消耗 57984.8 kg（比 Q2 微省） |
| 问题4 | 最优节流律 σ ≡ 1.0（满推力）→ 无路径约束时节流无增益 |

## 参考文献

1. Nair V S, Vaidyanathan A. Ascent trajectory design and optimization of a two-stage throttleable liquid rocket. *Advances in Space Research*, 2022.
2. Benedikter B, et al. Convex optimization of launch vehicle ascent trajectory with heat-flux and splash-down constraints. *Journal of Spacecraft and Rockets*, 2022.
3. Betts J T. Survey of numerical methods for trajectory optimization. *Journal of Guidance, Control, and Dynamics*, 1998.
4. 刘超越, 张成. 基于高斯伪谱法的二级助推战术火箭多阶段轨迹优化. 兵工学报, 2019.
5. 胡冬生 等. 含滑行时间约束的真空段弹道设计研究. 宇航总体技术, 2023.
6. 李惠峰 等. 液体火箭上升段制导方法的发展综述. 航天控制, 2023.
7. 余梦伦 等. 弹道学在运载火箭总体设计中的实践与展望. 宇航总体技术, 2023.

详见 `论文/main.tex` 参考文献节与 `文献/` 目录。
