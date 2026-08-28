"""问题 4：二子级发动机具备推力节流能力（60%~100% 额定推力）下，
重新设计滑行时间与二子级推力大小、俯仰角的优化控制策略，使燃料消耗最省。

优化问题：
    min  J = 二子级消耗推进剂 = integral mdotdt = integral (sigma(t)·Fmax2)/(Isp2·g0) dt
    s.t. 入轨条件 3 等式（|r|=a, |v|=sqrt(mu/a), r·v=0）
         节流比 sigma(t)  in  [0.6, 1]（分段线性参数化）
         俯仰角 phi(t) 分段线性参数化（首点取点火时刻速度方向）
         关机时刻 t_shut 自由（由入轨条件决定）

求解：控制参数化（phi 与 sigma 均 N 节点分段线性）+ 直接打靶 + SLSQP 多初值。
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

N_SIGMA = 4   # 节流比控制节点数
SIGMA_MIN = 0.6


class ThrottleOptimizer:
    """问题 4：节流 + 俯仰角联合优化（控制参数化 + SLSQP）。"""

    def __init__(self, rk: RocketParams | None = None,
                 n_phi: int = 5, n_sigma: int = N_SIGMA, rtol: float = 1e-9):
        self.rk = rk or RocketParams()
        self.n_phi = n_phi
        self.n_sigma = n_sigma
        self.rtol = rtol
        self._ign_cache: dict[float, tuple[np.ndarray, float]] = {}

    def ignition_state(self, t_coast: float) -> tuple[np.ndarray, float]:
        # 保留优化器施加的微小扰动，使滑行时间方向的数值导数可见。
        key = float(t_coast)
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
        """点火时刻惯性速度路径角 = 初始俯仰角基准（二子级段无大气参照）。"""
        return float(np.arctan2(
            x0[2] * x0[0] + x0[3] * x0[1],
            x0[2] * (-x0[1]) + x0[3] * x0[0],
        ))

    # -- 二级段积分 ----------------------------------------------------------
    def stage2_end(self, x0, t2, phi0, phi_nodes, sigma_nodes, t_shut):
        rk = self.rk
        n_phi, n_sig = self.n_phi, self.n_sigma
        tau_phi = np.linspace(0.0, 1.0, n_phi)
        tau_sig = np.linspace(0.0, 1.0, n_sig)
        phi_abs = np.concatenate(([phi0 / DEG], phi_nodes))       # 长度 n_phi
        sig_arr = np.asarray(sigma_nodes, float)                   # 长度 n_sig

        def interp(tau, tau_nodes, vals):
            return float(np.interp(tau, tau_nodes, vals))

        def rhs(t, x):
            tau = (t - t2) / (t_shut - t2)
            p = interp(tau, tau_phi, phi_abs) * DEG
            sig = float(np.clip(interp(tau, tau_sig, sig_arr), SIGMA_MIN, 1.0))
            tx_h, ty_h = thrust_unit(p, x[0], x[1])
            T = sig * rk.Fmax2
            r = np.hypot(x[0], x[1])
            ax = -MU / r**3 * x[0] + T / x[4] * tx_h
            ay = -MU / r**3 * x[1] + T / x[4] * ty_h
            mdot = -T / (rk.Isp2 * G0)
            return np.array([x[2], x[3], ax, ay, mdot])

        sol = solve_ivp(rhs, (t2, max(t_shut, t2 + 1e-3)), x0, method="DOP853",
                        rtol=self.rtol, atol=[1e-2, 1e-2, 1e-4, 1e-4, 1e-2])
        return sol.y[:, -1]

    # -- 目标 + 约束 --------------------------------------------------------
    def evaluate(self, z: np.ndarray):
        """z = [t_coast, phi_1..phi_{n_phi-1}, sigma_1..sigma_{n_sigma}, t_shut]"""
        t_coast = z[0]
        phi_nodes = z[1:1 + (self.n_phi - 1)]
        sigma_nodes = z[1 + (self.n_phi - 1): 1 + (self.n_phi - 1) + self.n_sigma]
        t_shut = z[-1]
        x0, t2 = self.ignition_state(t_coast)
        phi0 = self.phi0_of(x0)
        # 节流时允许比满推力燃烧更久；真正的容量约束是累计耗药量。
        max_duration = self.rk.t_burn2 / SIGMA_MIN
        if t_shut < t2 + 1.0 or t_shut > t2 + max_duration:
            return 1e6 + abs(t_shut - t2), np.array([1.0, 1.0, 1.0])
        xf = self.stage2_end(x0, t2, phi0, phi_nodes, sigma_nodes, t_shut)
        prop_used = x0[4] - xf[4]
        if prop_used > self.rk.mp2 + 1e-6:
            return 1e6 + (prop_used - self.rk.mp2), np.array([1.0, 1.0, 1.0])
        res = orbit_residual(xf, self.rk)
        return prop_used, res

    def objective(self, z: np.ndarray) -> float:
        prop_used, res = self.evaluate(z)
        penalty = 1e3 * float(np.sum(res**2))
        return prop_used + penalty

    def constraint_eq(self, z: np.ndarray) -> np.ndarray:
        _, res = self.evaluate(z)
        return res

    def _solve_one(self, z0: np.ndarray, maxiter: int, bounds) -> tuple | None:
        cons = {"type": "eq", "fun": lambda z: self.constraint_eq(z) * 1e2}
        opt = minimize(
            self.objective, z0, method="SLSQP", bounds=bounds,
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

    def solve(self, starts: list[np.ndarray], maxiter: int = 80, workers: int = 8):
        from concurrent.futures import ProcessPoolExecutor
        n_phi, n_sig = self.n_phi, self.n_sigma
        bounds = [(0.0, 400.0)] + [(-90.0, 90.0)] * (n_phi - 1) \
                 + [(SIGMA_MIN, 1.0)] * n_sig + [(0.0, 1e3)]
        results = []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(self._solve_one, z0, maxiter, bounds) for z0 in starts]
            for fut in futures:
                r = fut.result()
                if r is not None:
                    results.append(r)
        results.sort(key=lambda r: r[1])
        return results


def make_starts_q4(rk, n_phi, n_sigma) -> list[np.ndarray]:
    """初值：Q3 最优解（全推力）+ 节流 0.6/0.8 + 多滑行时间组合。"""
    s = Simulator(rk)
    starts = []
    opt0 = ThrottleOptimizer(rk, n_phi, n_sigma)
    # 1) 从 Q3 解构造（若存在）：t_c=90, 线性俯仰律, sigma=1
    try:
        q3 = pd.read_csv(RES_DIR / "q3_summary.csv").iloc[0]
        tc_q3 = float(q3["t_coast_s"])
        ts_q3 = float(q3["t_shut_s"])
        phi0_q3 = float(q3["phi0_deg"])
        phi_nodes_q3 = [float(v) for v in str(q3["phi_nodes_deg"]).split(" | ")]
        for sig_level in [0.6, 0.8, 0.9, 1.0]:
            z = np.zeros(1 + (n_phi - 1) + n_sigma + 1)
            z[0] = tc_q3
            z[1:1 + (n_phi - 1)] = phi_nodes_q3[:n_phi - 1]
            z[1 + (n_phi - 1): 1 + (n_phi - 1) + n_sigma] = sig_level
            z[-1] = ts_q3 * (0.95 if sig_level < 1.0 else 1.0)
            starts.append(z)
    except Exception:
        pass
    # 2) 广泛网格
    for tc in [30.0, 60.0, 90.0, 113.0, 150.0, 200.0, 260.0]:
        res = s.simulate(t_coast=tc, controller=None, t_shut=rk.t_burn1 + tc + 1e-6)
        x0 = res.stages[-1].x[0]
        phi0 = opt0.phi0_of(x0)
        t2 = float(res.stages[-1].t[0])
        for sig_level in [0.6, 0.8, 1.0]:
            for k0 in [-0.2, -0.1, -0.05, 0.0]:
                z = np.zeros(1 + (n_phi - 1) + n_sigma + 1)
                z[0] = tc
                phi_land = 4.0
                for j in range(n_phi - 1):
                    tau = (j + 1) / n_phi
                    z[1 + j] = phi0 / DEG + (phi_land - phi0 / DEG) * tau
                for j in range(n_sigma):
                    z[1 + (n_phi - 1) + j] = sig_level
                z[-1] = t2 + 0.85 * rk.t_burn2
                starts.append(z)
    return starts


def main() -> None:
    rk = RocketParams()
    opt = ThrottleOptimizer(rk)
    print("=" * 72)
    print("问题 4：推力节流（60%~100%）下二子级燃料最省")
    print(f"控制：俯仰角 {opt.n_phi} 节点分段线性 + 节流比 {opt.n_sigma} 节点分段线性")
    print("=" * 72)

    starts = make_starts_q4(rk, opt.n_phi, opt.n_sigma)
    print(f"初值数量：{len(starts)}")
    results = opt.solve(starts)
    print(f"收敛解数量：{len(results)}")
    if not results:
        print("未找到满足入轨条件的解！")
        return

    z_best, prop_used, nres = results[0]
    prop_used, res = opt.evaluate(np.array(z_best, float))
    z_best = np.array(z_best)
    t_coast = z_best[0]
    phi_nodes = z_best[1:1 + (opt.n_phi - 1)]
    sigma_nodes = z_best[1 + (opt.n_phi - 1): 1 + (opt.n_phi - 1) + opt.n_sigma]
    t_shut = z_best[-1]
    x0, t2 = opt.ignition_state(t_coast)
    phi0 = opt.phi0_of(x0)

    print("-" * 72)
    print("最优解：")
    print(f"  滑行时间    t_c     = {t_coast:.4f} s")
    print(f"  关机时刻    t_shut  = {t_shut:.4f} s（点火 t2 = {t2:.4f} s，"
          f"燃尽上限 {t2 + rk.t_burn2:.4f} s）")
    print(f"  俯仰角节点（deg）: ", end="")
    tau_phi = np.linspace(0.0, 1.0, opt.n_phi)
    for k, v in zip(tau_phi, np.concatenate(([phi0 / DEG], phi_nodes))):
        print(f" [{k:.2f}]:{v:.3f}", end="")
    print()
    print(f"  节流比节点: ", end="")
    tau_sig = np.linspace(0.0, 1.0, opt.n_sigma)
    for k, v in zip(tau_sig, sigma_nodes):
        print(f" [{k:.2f}]:{v:.3f}", end="")
    print()
    print(f"  消耗推进剂  = {prop_used:.1f} kg")
    print(f"  关机质量    = {x0[4] - prop_used:.1f} kg（剩余推进剂 "
          f"{x0[4] - prop_used - rk.ms2 - rk.m_payload:.1f} kg）")
    print(f"  入轨残差    = {nres:.2e}")

    # 与 Q3 对比
    try:
        q3 = pd.read_csv(RES_DIR / "q3_summary.csv").iloc[0]
        print("-" * 72)
        print(f"对比问题 3（无节流，全推力）：消耗推进剂 {q3['propellant_used_kg']:.1f} kg")
        print(f"            本问节省 {q3['propellant_used_kg'] - prop_used:.1f} kg")
    except Exception:
        print("（Q3 结果未找到，跳过对比）")

    # ---- 完整仿真验证 ----
    def controller(t):
        tau = (t - t2) / (t_shut - t2)
        phi_abs = np.concatenate(([phi0 / DEG], phi_nodes))
        return float(np.interp(tau, tau_phi, phi_abs) * DEG)

    def sigma_fn(t):
        tau = (t - t2) / (t_shut - t2)
        return float(np.clip(np.interp(tau, tau_sig, sigma_nodes), SIGMA_MIN, 1.0))

    sim = Simulator(rk)
    res = sim.simulate(t_coast=t_coast, controller=controller,
                       sigma=sigma_fn, t_shut=t_shut)
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
        "sigma_nodes": " | ".join(f"{v:.4f}" for v in sigma_nodes),
        "phi0_deg": phi0 / DEG,
        "propellant_used_kg": prop_used,
        "m_final_kg": fs_["m_kg"],
        "propellant_remaining_kg": fs_["m_kg"] - rk.ms2 - rk.m_payload,
        "residual_norm": nres,
    }]).to_csv(RES_DIR / "q4_summary.csv", index=False)
    df = pd.DataFrame({
        "t_s": res.t, "h_km": res.h / 1e3, "V_mps": res.v_in,
        "gamma_deg": res.gamma / DEG, "m_kg": res.m,
        "phi_deg": res.phi / DEG, "sigma": res.sigma, "phase": res.phase,
    })
    df.to_csv(RES_DIR / "q4_trajectory.csv", index=False)

    # ---- 绘图 ----
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    ax = axes[0, 0]
    ax.plot(res.t / 60, res.h / 1e3, lw=1.6)
    ax.axhline(400, ls="--", color="gray")
    ax.set_xlabel("t (min)"); ax.set_ylabel("h (km)")
    ax.set_title("高度"); ax.grid(alpha=0.4)
    ax = axes[0, 1]
    ax.plot(res.t / 60, res.v_in / 1e3, lw=1.6, color="tab:red")
    ax.axhline(rk.v_circular / 1e3, ls="--", color="gray")
    ax.set_xlabel("t (min)"); ax.set_ylabel("V (km/s)")
    ax.set_title("惯性速度"); ax.grid(alpha=0.4)
    ax = axes[0, 2]
    ax.plot(res.t / 60, res.gamma / DEG, lw=1.6, color="tab:green")
    ax.axhline(0, ls="--", color="gray")
    ax.set_xlabel("t (min)"); ax.set_ylabel("gamma (deg)")
    ax.set_title("飞行路径角"); ax.grid(alpha=0.4)
    ax = axes[1, 0]
    mask = res.t >= t2
    ax.plot(res.t[mask] / 60, res.phi[mask] / DEG, lw=1.6, color="tab:purple")
    ax.set_xlabel("t (min)"); ax.set_ylabel("俯仰角 phi (deg)")
    ax.set_title("二子级俯仰角控制律"); ax.grid(alpha=0.4)
    ax = axes[1, 1]
    ax.plot(res.t[mask] / 60, res.sigma[mask], lw=1.6, color="tab:orange")
    ax.axhline(0.6, ls="--", color="gray", label="0.6 下限")
    ax.axhline(1.0, ls="--", color="gray", label="1.0 上限")
    ax.set_xlabel("t (min)"); ax.set_ylabel("节流比 sigma")
    ax.set_title("节流比控制律"); ax.legend(fontsize=8); ax.grid(alpha=0.4)
    ax = axes[1, 2]
    ax.plot(res.t[mask] / 60, res.m[mask] / 1e3, lw=1.6, color="tab:brown")
    ax.set_xlabel("t (min)"); ax.set_ylabel("m (t)")
    ax.set_title("二子级质量"); ax.grid(alpha=0.4)
    fig.suptitle("问题 4：推力节流下燃料最省最优解", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG_DIR / "q4_throttle_optimal.png", dpi=150)
    plt.close(fig)
    print(f"已保存：{RES_DIR / 'q4_summary.csv'}, {RES_DIR / 'q4_trajectory.csv'}, "
          f"{FIG_DIR / 'q4_throttle_optimal.png'}")


if __name__ == "__main__":
    main()
