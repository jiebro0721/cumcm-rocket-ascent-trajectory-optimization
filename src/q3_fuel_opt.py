"""问题 3：设计滑行时间与二子级俯仰角优化控制策略，使二子级燃料消耗最省。

优化问题：
    min  J = 二子级消耗推进剂 = m(t2+) - m(t_shut)   （等价 max m(t_shut)）
    s.t. 入轨条件 3 等式（|r|=a, |v|=sqrt(mu/a), r·v=0）
         俯仰角 phi(t) 以分段线性控制律参数化（N 段节点值）
         关机时刻 t_shut 自由（由入轨条件决定，可提前关机）

求解：控制参数化（分段线性）+ 直接打靶 + SLSQP（带入轨等式约束），
     多初值（含问题 2 的解作为初值）避免局部解。
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

from common import (
    DEG,
    G0,
    MU,
    RES_DIR,
    FIG_DIR,
    RocketParams,
    Simulator,
    orbit_residual,
    thrust_unit,
)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False

N_PHI = 5          # 俯仰角控制节点数（分段线性）
N_GRID_TIMES = 6   # 初值扫描点火前网格数量


class FuelOptimizer:
    """问题 3：燃料最省优化器（控制参数化 + SLSQP）。"""

    def __init__(self, rk: RocketParams | None = None, n_phi: int = N_PHI, rtol: float = 1e-9):
        self.rk = rk or RocketParams()
        self.n_phi = n_phi
        self.rtol = rtol
        self._ign_cache: dict[float, tuple[np.ndarray, float]] = {}

    def ignition_state(self, t_coast: float) -> tuple[np.ndarray, float]:
        key = round(t_coast, 6)
        if key not in self._ign_cache:
            sim = Simulator(self.rk)
            res = sim.simulate(
                t_coast=key, controller=None,
                t_shut=self.rk.t_burn1 + key + 1e-6, rtol=self.rtol,
            )
            seg = res.stages[-1]
            self._ign_cache[key] = (seg.x[0].copy(), float(seg.t[0]))
        return self._ign_cache[key]

    def phi0_of(self, x0: np.ndarray) -> float:
        """点火时刻惯性速度路径角 = 初始俯仰角基准。

        二子级段无大气阻力，不存在相对大气参照，初始俯仰角按惯性速度方向确定。
        """
        return float(np.arctan2(
            x0[2] * x0[0] + x0[3] * x0[1],
            x0[2] * (-x0[1]) + x0[3] * x0[0],
        ))

    # -- 二级段积分（phi 用控制节点线性插值，节点时间线性分布于 [t2, t_shut]）--
    def stage2_end(self, x0, t2, phi0, phi_nodes, t_shut):
        rk = self.rk
        n = self.n_phi
        # 节点时刻（相对 t2），等距分布于 [0, t_shut-t2]
        tau_nodes = np.linspace(0.0, 1.0, n)
        # 节点值为绝对俯仰角 [deg]，首点固定为 phi0，其余节点为优化变量
        phi_abs = np.concatenate(([phi0 / DEG], phi_nodes))   # 长度 n

        def phi_at(tau):
            return np.interp(tau, tau_nodes, phi_abs) * DEG

        def rhs(t, x):
            tau = (t - t2) / (t_shut - t2)
            p = phi_at(tau)
            tx_h, ty_h = thrust_unit(p, x[0], x[1])
            T = rk.Fmax2
            r = np.hypot(x[0], x[1])
            ax = -MU / r**3 * x[0] + T / x[4] * tx_h
            ay = -MU / r**3 * x[1] + T / x[4] * ty_h
            return np.array([x[2], x[3], ax, ay, -T / (rk.Isp2 * G0)])

        sol = solve_ivp(rhs, (t2, max(t_shut, t2 + 1e-3)), x0, method="DOP853",
                        rtol=self.rtol, atol=[1e-2, 1e-2, 1e-4, 1e-4, 1e-2])
        return sol.y[:, -1]

    # -- 目标 + 约束 --------------------------------------------------------
    def evaluate(self, z: np.ndarray):
        """决策向量 z = [t_coast, phi_1..phi_{n-1}, t_shut]。返回 (目标, 残差)。"""
        t_coast = z[0]
        phi_nodes = z[1:-1]
        t_shut = z[-1]
        x0, t2 = self.ignition_state(t_coast)
        phi0 = self.phi0_of(x0)
        m_ignit = x0[4]
        if t_shut < t2 + 1.0 or t_shut > t2 + self.rk.t_burn2:
            # 罚：理想目标给大值
            return 1e6 + abs(t_shut - t2), np.array([1.0, 1.0, 1.0])
        xf = self.stage2_end(x0, t2, phi0, phi_nodes, t_shut)
        prop_used = m_ignit - xf[4]
        res = orbit_residual(xf, self.rk)
        return prop_used, res

    def objective(self, z: np.ndarray) -> float:
        prop_used, res = self.evaluate(z)
        # 罚函数：入轨残差违反时惩罚（SLSQP 也用约束，双保险）
        penalty = 1e3 * float(np.sum(res**2))
        return prop_used + penalty

    def constraint_eq(self, z: np.ndarray) -> np.ndarray:
        _, res = self.evaluate(z)
        return res

    # -- 求解 ---------------------------------------------------------------
    def _solve_one(self, z0: np.ndarray, maxiter: int) -> tuple | None:
        """单初值 SLSQP（供并行调用）。"""
        cons = {"type": "eq", "fun": lambda z: self.constraint_eq(z) * 1e2}
        opt = minimize(
            self.objective, z0, method="SLSQP",
            bounds=[(0.0, 400.0)] + [(-90.0, 90.0)] * (self.n_phi - 1)
                   + [(0.0, 1e3)],
            constraints=[cons],
            options={"ftol": 1e-8, "maxiter": maxiter, "disp": False},
        )
        # 候选解：SLSQP 终值优先，其次初值（若初值本身可行）
        cands = [opt.x]
        prop0, res0 = self.evaluate(z0)
        if np.linalg.norm(res0) < 1e-5:
            cands.append(z0)
        best = None
        for cand in cands:
            prop_used, res = self.evaluate(cand)
            if np.linalg.norm(res) < 1e-5:
                if best is None or prop_used < best[1]:
                    best = (cand, prop_used, float(np.linalg.norm(res)))
        return best

    def solve(self, starts: list[np.ndarray], maxiter: int = 60, workers: int = 8):
        from concurrent.futures import ProcessPoolExecutor
        results = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(self._solve_one, z0, maxiter) for z0 in starts]
            for fut in futures:
                r = fut.result()
                if r is not None:
                    results.append(r)
        results.sort(key=lambda r: r[1])
        return results


def make_starts(rk, n_phi) -> list[np.ndarray]:
    """初值集合：覆盖滑行时间 0~300s、俯仰角律多种形状（线性/凸/凹）。

    特别地，将 Q2 打靶解作为初值之一，保证优化从已知可行解出发。
    """
    s = Simulator(rk)
    starts = []
    bt = rk.t_burn2
    opt0 = FuelOptimizer(rk)

    # 0) Q2 打靶解（读取结果目录，若无则跳过）
    try:
        q2 = pd.read_csv(RES_DIR / "q2_summary.csv").iloc[0]
        tc_q2 = float(q2["t_coast_s"])
        ts_q2 = float(q2["t_shut_s"])
        k_q2 = float(q2["k_deg_s"])
        res0 = s.simulate(t_coast=tc_q2, controller=None,
                          t_shut=rk.t_burn1 + tc_q2 + 1e-6)
        phi0_q2 = opt0.phi0_of(res0.stages[-1].x[0])
        t2_q2 = float(res0.stages[-1].t[0])
        for k_off in [0.0, 0.5, -0.5]:   # k 与其抖动
            z = np.zeros(n_phi + 1)
            z[0] = tc_q2
            for j in range(n_phi - 1):
                tau = (j + 1) / n_phi
                z[1 + j] = phi0_q2 / DEG + (k_q2 + k_off * 0.01) * tau * (ts_q2 - t2_q2)
            z[-1] = ts_q2
            starts.append(z)
    except Exception:
        pass

    # 1) 网格初值
    for tc in [0.0, 60.0, 90.0, 113.0, 130.0, 160.0, 200.0, 250.0, 300.0]:
        res = s.simulate(t_coast=tc, controller=None, t_shut=rk.t_burn1 + tc + 1e-6)
        x0 = res.stages[-1].x[0]
        phi0 = opt0.phi0_of(x0)
        t2 = float(res.stages[-1].t[0])
        for shape in [-1.0, -0.3, 0.0, 0.3, 1.0]:
            z = np.zeros(n_phi + 1)
            z[0] = tc
            # phi(tau) = phi0 + shape * (phi_land - phi0) * tau^p，p=1 线性、p>1 凸、p<1 凹
            # 取终端俯仰角目标 phi_land 为 5deg 附近（gamma 需归零），用形状参数生成
            phi_land = 4.0
            for j in range(n_phi - 1):
                tau = (j + 1) / n_phi
                p = 1.0 if shape == 0 else (2.0 if shape > 0 else 0.5)
                z[1 + j] = phi0 / DEG + (phi_land - phi0 / DEG) * tau**p
            z[-1] = t2 + 0.85 * bt
            starts.append(z)
    return starts


def main() -> None:
    rk = RocketParams()
    opt = FuelOptimizer(rk)
    print("=" * 72)
    print("问题 3：二子级燃料最省（滑行时间 + 俯仰角分段线性控制）")
    print(f"控制节点数 N_phi = {opt.n_phi}（首点固定为初始俯仰角）")
    print("=" * 72)

    starts = make_starts(rk, opt.n_phi)
    print(f"初值数量：{len(starts)}")
    results = opt.solve(starts)
    print(f"收敛解数量：{len(results)}")
    if not results:
        print("未找到满足入轨条件的解！")
        return

    z_best, prop_used, nres = results[0]
    # 双精度重算并输出
    z_best = np.array(z_best, float)
    prop_used, res = opt.evaluate(z_best)
    t_coast, phi_nodes, t_shut = z_best[0], z_best[1:-1], z_best[-1]
    x0, t2 = opt.ignition_state(t_coast)
    phi0 = opt.phi0_of(x0)

    print("-" * 72)
    print(f"最优解：")
    print(f"  滑行时间    t_c    = {t_coast:.4f} s")
    print(f"  关机时刻    t_shut = {t_shut:.4f} s（点火 t2 = {t2:.4f} s，"
          f"燃尽上限 {t2 + rk.t_burn2:.4f} s）")
    print(f"  俯仰角节点（deg）：", end="")
    tau_nodes = np.linspace(0.0, 1.0, opt.n_phi)
    for k, v in zip(tau_nodes, np.concatenate(([phi0 / DEG], phi_nodes))):
        print(f" [{k:.2f}]:{v:.3f}", end="")
    print()
    print(f"  二子级消耗推进剂 = {prop_used:.1f} kg")
    print(f"  关机质量        = {x0[4] - prop_used:.1f} kg（剩余推进剂 "
          f"{x0[4] - prop_used - rk.ms2 - rk.m_payload:.1f} kg）")
    print(f"  入轨残差 |res|   = {nres:.2e}")

    # 与问题 2 解对比
    q2 = pd.read_csv(RES_DIR / "q2_summary.csv").iloc[0]
    print("-" * 72)
    print(f"对比问题 2（恒定俯仰角速率 k = {q2['k_deg_s']:.4f} deg/s）：")
    print(f"  Q2 消耗推进剂 = {rk.mp2 - q2['propellant_remaining_kg']:.1f} kg")
    print(f"  Q3 消耗推进剂 = {prop_used:.1f} kg（节省 {rk.mp2 - q2['propellant_remaining_kg'] - prop_used:.1f} kg）")

    # ---- 完整仿真重新验证 ----
    def controller(t):
        tau = (t - t2) / (t_shut - t2)
        phi_abs = np.concatenate(([phi0 / DEG], phi_nodes))
        return float(np.interp(tau, tau_nodes, phi_abs) * DEG)

    sim = Simulator(rk)
    res = sim.simulate(t_coast=t_coast, controller=controller, t_shut=t_shut)
    fs_ = res.final_state()
    print("-" * 72)
    print("完整仿真验证：")
    print(f"  h = {fs_['h_km']:.4f} km, V = {fs_['v_in_mps']:.3f} m/s, "
          f"gamma = {fs_['gamma_deg']:.5f} deg")
    print(f"  入轨残差 = {orbit_residual(res.x[-1], rk)}")

    # ---- 保存 ----
    pd.DataFrame([{
        "t_coast_s": t_coast, "t_shut_s": t_shut,
        "phi_nodes_deg": " | ".join(f"{v:.4f}" for v in phi_nodes),
        "phi0_deg": phi0 / DEG,
        "propellant_used_kg": prop_used,
        "m_final_kg": fs_["m_kg"],
        "propellant_remaining_kg": fs_["m_kg"] - rk.ms2 - rk.m_payload,
        "residual_norm": nres,
    }]).to_csv(RES_DIR / "q3_summary.csv", index=False)
    df = pd.DataFrame({
        "t_s": res.t, "h_km": res.h / 1e3, "V_mps": res.v_in,
        "gamma_deg": res.gamma / DEG, "m_kg": res.m,
        "phi_deg": res.phi / DEG, "phase": res.phase,
    })
    df.to_csv(RES_DIR / "q3_trajectory.csv", index=False)

    # ---- 绘图：控制律与轨迹 ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax = axes[0, 0]
    ax.plot(res.t / 60, res.h / 1e3, lw=1.6)
    ax.axhline(400, ls="--", color="gray", label="目标轨道")
    ax.set_xlabel("t (min)"); ax.set_ylabel("h (km)"); ax.legend(fontsize=8)
    ax.set_title("高度"); ax.grid(alpha=0.4)
    ax = axes[0, 1]
    ax.plot(res.t / 60, res.v_in / 1e3, lw=1.6, color="tab:red")
    ax.axhline(rk.v_circular / 1e3, ls="--", color="gray", label="圆轨道速度")
    ax.set_xlabel("t (min)"); ax.set_ylabel("V (km/s)"); ax.legend(fontsize=8)
    ax.set_title("惯性速度"); ax.grid(alpha=0.4)
    ax = axes[1, 0]
    ax.plot(res.t / 60, res.gamma / DEG, lw=1.6, color="tab:green")
    ax.axhline(0, ls="--", color="gray")
    ax.set_xlabel("t (min)"); ax.set_ylabel("gamma (deg)")
    ax.set_title("飞行路径角"); ax.grid(alpha=0.4)
    ax = axes[1, 1]
    ax.plot(res.t[res.t >= t2] / 60, res.phi[res.t >= t2] / DEG, lw=1.6, color="tab:purple")
    ax.set_xlabel("t (min)"); ax.set_ylabel("俯仰角 phi (deg)")
    ax.set_title("二子级俯仰角控制律（分段线性）"); ax.grid(alpha=0.4)
    fig.suptitle("问题 3：燃料最省最优轨迹与最优控制", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG_DIR / "q3_fuel_optimal.png", dpi=150)
    plt.close(fig)
    print(f"已保存：{RES_DIR / 'q3_summary.csv'}, {RES_DIR / 'q3_trajectory.csv'}, "
          f"{FIG_DIR / 'q3_fuel_optimal.png'}")


if __name__ == "__main__":
    main()
