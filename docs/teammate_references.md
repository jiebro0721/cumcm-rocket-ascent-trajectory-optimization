# 两级运载火箭发射轨迹建模与仿真 —— 相关文献清单

题目要点：两级液体火箭（一子级 0–10 s 垂直上升 + 程序转弯、无动力滑行、二子级入轨修正），
目标 400 km 近地圆轨道；涉及变质量动力学、旋转地球、指数大气模型（ρ = ρ0·e^(−h/h0)）、气动阻力、
比冲/推力模型；问题 2–4 分别是参数调整入轨、滑行时间与俯仰角速率优化（燃料最省）、
推力节流（60%–100%）下再优化。

---

## 一、动力学建模基础（对应问题 1 的建模与仿真）

1. Cornelisse J. W., Schöyer H. F. R., Wakker K. F. *Rocket Propulsion and Spaceflight Dynamics*. Pitman, 1979.
   —— 运载火箭飞行动力学经典教材：发射弹道、程序转弯/重力转弯、多级火箭质量方程、旋转地球与大气模型。

2. Tewari A. *Atmospheric and Space Flight Dynamics: Modeling and Simulation with MATLAB and Simulink*. Birkhäuser, 2007.
   —— 含变质量刚体动力学、指数大气、升空段轨迹仿真，可直接对应本题数值仿真实现。
   （书目信息：[KAIST 图书馆](http://library.kaist.ac.kr/search/ctlgSearch/posesn/view.do?bibctrlno=329249&ty=B)）

3. Sutton G. P., Biblarz O. *Rocket Propulsion Elements*, 9th ed. Wiley, 2017.
   —— 比冲、推力、质量比、推进剂流量等发动机参数的建模依据。

4. Thomson W. T. *Introduction to Space Dynamics*. Dover, 1986.
   —— 轨道力学基础，入轨圆轨道条件（v = √(μ/r) 等）的来源。

5. 贾沛然, 陈克俊, 何力. 《远程火箭弹道学》. 国防科技大学出版社, 1993.
   —— 中文经典：主动段弹道计算、程序转弯、重力转弯、关机点参数确定，与本题结构几乎一一对应。

6. 陈克俊, 刘鲁华, 孟云鹤. 《远程火箭飞行动力学与制导》. 国防工业出版社.
   —— 上升段动力学建模与制导方案设计。

7. 大气模型参考：*US Standard Atmosphere, 1976*（本题采用的指数模型是其简化形式）。

## 二、重力转弯 / 程序转弯（对应问题 1–2 的转弯段）

8. Culler G. J., Fried B. D. "Universal Gravity Turn Trajectories." *Journal of Applied Physics*, 28(6): 672–676, 1957.
   —— 重力转弯轨迹的原始经典文献。
   （[AIP 链接](https://pubs.aip.org/aip/jap/article-abstract/28/6/672/161544/Universal-Gravity-Turn-Trajectories?redirectedFrom=PDF)）

9. "A trade-off methodology for micro-launchers." *Aerospace Systems*, 2021.
   —— 明确给出两级入轨各阶段（lift-off → pitch over → constant pitch → first-stage burn out → coast → second-stage ignition），
   与本题三段式飞行结构直接对应。
   （[Springer 链接](https://rd.springer.com/article/10.1007/s42401-021-00095-w)）

## 三、轨迹优化理论（对应问题 3–4 的燃料最省）

10. Betts J. T. *Practical Methods for Optimal Control and Estimation Using Nonlinear Programming*, 2nd ed. SIAM, 2010.
    —— 直接法（打靶/配点/伪谱）求解最优控制的权威参考。

11. Longuski J. M., Guzmán J. J., Prussing J. E. *Optimal Control with Aerospace Applications*. Springer, 2014.
    —— 庞特里亚金极大值原理、奇异弧（推力节流时最优解常处于奇异弧上）等航天最优控制理论。

12. Vinh N. X. *Optimal Trajectories in Atmospheric Flight*. Elsevier, 1975.
    —— 大气飞行最小燃料/最小时间上升轨迹的经典专著（Goddard 问题等）。

13. Bryson A. E., Ho Y. C. *Applied Optimal Control*. Hemisphere, 1975.
    —— 最优控制与两点边值问题经典教材。

14. Goddard 火箭问题相关文献：燃料最省垂直上升/上升轨迹的经典问题，是本题问题 3 的简化原型
    （见 Beispiele der Optimalen Steuerung: <https://www.math.uni-bremen.de/zetem/alt/optimmedia/webcontrol/rakete2.html>）。

## 四、运载火箭上升段轨迹优化论文（与问题 2–4 直接相关）

15. Lu P., Griffin B., Dukeman G., Chavez F. "Rapid Optimal Multiburn Ascent Planning and Guidance." AIAA GNC Conference, 2008 (AIAA-2008-6219).
    —— 多级火箭“点火—滑行—再点火”的最优上升段规划与制导，正好覆盖本题“主动段—滑行—入轨修正”结构。
    （[Semantic Scholar](https://www.semanticscholar.org/paper/Rapid-Optimal-Multiburn-Ascent-Planning-and-Lu-Griffin/82d43b6588072cba038dca85912c797fdfbc5b26)）

16. Lu P., Pan B. "Highly Constrained Optimal Launch Ascent Guidance." *Journal of Guidance, Control, and Dynamics*, 2010.
    （[Semantic Scholar](https://www.semanticscholar.org/paper/Highly-Constrained-Optimal-Launch-Ascent-Guidance-Lu-Pan/bc79cd29e67eaddcc741e77db3064633a731c6c8)）

17. Jezewski D. J. "Optimal Analytic Multiburn Trajectories." *Journal of Guidance, Control, and Dynamics*, 15(3), 1992.
    —— 多脉冲（多段点火）最优轨迹的解析方法。
    （[Semantic Scholar](https://www.semanticscholar.org/paper/Optimal-analytic-multiburn-trajectories-Jezewski/596f8836f07686f7632208d39e2faf3acc15846b)）

18. Dukeman G. A. "Atmospheric Ascent Guidance for Rocket-Powered Launch Vehicles." AIAA GNC Conference, 2002 (AIAA-2002-4558).
    —— 大气上升段制导，含过载/动压约束处理。

19. "Launch vehicle ascent trajectory optimization based on modified Gauss Pseudospectral method." IEEE.
    —— 高斯伪谱法用于上升段轨迹优化。
    （[infona](https://www.infona.pl/resource/bwmeta1.element.ieee-art-000006852246)）

20. Garg D., Patterson M., Hager W. W., Rao A. V., Benson D. A., Huntington G. T. "Direct trajectory optimization and costate estimation via an orthogonal collocation method." *Journal of Guidance, Control, and Dynamics*, 33(5), 2010.
    —— 伪谱配点法（GPOPS 系列软件的理论基础，可用于问题 3、4 的直接求解）。

21. "Ascent trajectory optimization with singular arc using linear Gauss pseudospectral model predictive control." *Aerospace Science and Technology*（2026 在线）.
    —— 推力可调（奇异弧）上升段轨迹优化，直接对应问题 4 的推力节流情形。
    （[ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1270963826006978)）

22. "Optimal Conceptual Design of Two-Stage Reusable Rocket Vehicles Including Trajectory Optimization." *Journal of Spacecraft and Rockets*, 41(5), 2004.
    —— 两级火箭总体设计与轨迹优化一体化。
    （[AIAA](https://arc.aiaa.org/doi/10.2514/1.1082)）

23. "The two-point boundary-value problem for rocket trajectories."（间接法求解火箭轨迹两点边值问题）
    （[ULisboa](https://researchportal.ulisboa.pt/en/publications/the-two-point-boundary-value-problem-for-rocket-trajectories/)）

## 五、综述文献

24. "A Review on Trajectory Optimisation Techniques for Launch Vehicles during Ascent Phase." IEEE, 2023.
    —— 上升段轨迹优化方法综述（直接法/间接法/智能算法对比）。
    （[IEEE](https://ieeexplore.ieee.org/document/10039742)）

25. "A Review of Intelligent Trajectory Planning and Optimization for Aerospace Vehicles."
    —— 智能轨迹规划综述。
    （[西北工业大学 Pure](https://pure.nwpu.edu.cn/en/publications/a-review-of-intelligent-trajectory-planning-and-optimization-for-/)）

## 六、中文文献

26. 运载火箭一阶梯度最优弹道在线规划方法. 宇航学报（spacejournal）.
    （[链接](https://www.spacejournal.cn/cn/article/id/ea195688-d010-4b4e-8bc7-c5059b3e6379)）

27. 基于 Gauss 伪谱法的火箭飞行轨迹优化求解方法研究.
    （[CSDN 文库](https://wenku.csdn.net/doc/3d7ckxtdfe)）

28. RBCC 可重复使用运载器上升段轨迹优化设计. 被引 18.
    （[维普](http://dianda.cqvip.com/Qikan/Article/Detail?id=42657007)）

29. 垂直起降火箭入轨回收一体轨迹优化.
    （[国家科技期刊平台](https://search.napstic.cn/literature/periodical/010fxlx202405004)）

---

## 附：可直接使用的求解工具

- **GPOPS-II**（Patterson & Rao, ACM TOMS 2014）：MATLAB 下的 hp-自适应高斯伪谱法最优控制软件，
  适合直接求解问题 3、4（配合 IPOPT/SNOPT 求解 NLP）。
- MATLAB/Simulink 数值积分（ode45 等）完成问题 1 的基准仿真。

> 说明：标有链接的条目已通过检索核实出处；个别条目（如 18、20、22）为学界公认文献，
> 引用前建议在数据库（AIAA/Web of Science/知网）中二次核对卷期页码。
