"""问题 2：入轨条件反演 —— 调整滑行时间与二子级俯仰角变化率（可提前关机）。

控制策略（题目给定）：
- 二子级推力俯仰角以恒定速率 k 变化，初始俯仰角与（相对大气）速度方向一致；
- 二子级维持额定全推力；
- 关机时刻 t_shut 由入轨条件决定（可提前关机）；
- 目标：进入 400 km 近地圆轨道。

入轨条件（惯性系，3 个等式）：
    |r| = R_E + H,  |v| = sqrt(mu/|r|),  r·v = 0

未知量 (t_coast, k, t_shut) 共 3 个、方程 3 个 -> 三维打靶（fsolve 多初值）。
点火前轨迹只与 t_coast 有关，缓存之；二级段单独快速积分。
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

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


class InscriptionSolver:
    """问题 2 三维打靶求解器。"""

    def __init__(self, rk: RocketParams | None = None, rtol: float = 1e-9):
        self.rk = rk or RocketParams()
        self.rtol = rtol
        self._ign_cache: dict[float, tuple[np.ndarray, float]] = {}

    # -- 点火前状态缓存（只与 t_coast 有关）-------------------------------
    def ignition_state(self, t_coast: float) -> tuple[np.ndarray, float]:
        # 连续优化变量不能取整作为缓存键，否则有限差分扰动会被抹掉。
        key = float(t_coast)
        if key not in self._ign_cache:
            sim = Simulator(self.rk)
            res = sim.simulate(
                t_coast=key, controller=None,
                t_shut=self.rk.t_burn1 + key + 1e-6,
                rtol=self.rtol,
            )
            seg = res.stages[-1]
            self._ign_cache[key] = (seg.x[0].copy(), float(seg.t[0]))
        return self._ign_cache[key]

    # -- 二级段单次积分 ----------------------------------------------------
    def stage2_end(
        self, x0: np.ndarray, t2: float, phi0: float,
        k_deg_s: float, t_shut: float,
    ) -> np.ndarray:
        k = k_deg_s * DEG
        rk = self.rk

        def rhs(t, x):
            p = phi0 + k * (t - t2)
            tx_h, ty_h = thrust_unit(p, x[0], x[1])
            T = rk.Fmax2
            r = np.hypot(x[0], x[1])
            ax = -MU / r**3 * x[0] + T / x[4] * tx_h
            ay = -MU / r**3 * x[1] + T / x[4] * ty_h
            return np.array([x[2], x[3], ax, ay, -T / (rk.Isp2 * G0)])

        sol = solve_ivp(
            rhs, (t2, t_shut), x0, method="DOP853",
            rtol=self.rtol, atol=[1e-2, 1e-2, 1e-4, 1e-4, 1e-2],
        )
        return sol.y[:, -1]

    # -- 残差 --------------------------------------------------------------
    def residual(self, z: np.ndarray) -> np.ndarray:
        t_coast, k_deg_s, t_shut = z
        x0, t2 = self.ignition_state(t_coast)
        phi0 = np.arctan2(x0[2] * x0[0] + x0[3] * x0[1],
                          x0[2] * (-x0[1]) + x0[3] * x0[0])
        if t_shut > t2 + self.rk.t_burn2:   # 超出燃尽上限 -> 罚
            return np.array([1.0, 1.0, 0.0])
        xf = self.stage2_end(x0, t2, phi0, k_deg_s, t_shut)
        return orbit_residual(xf, self.rk)

    # -- 多初值打靶 ---------------------------------------------------------
    def solve(self, starts=None) -> list[tuple[np.ndarray, float]]:
        if starts is None:
            starts = []
            for tc in [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0]:
                for k in [-0.4, -0.2, -0.1, -0.05, 0.0, 0.1, 0.2]:
                    for ts in [350.0, 450.0, 550.0, 650.0]:
                        starts.append((tc, k, ts))
        sols = {}
        for z0 in starts:
            # 可行域检查：燃尽上限保护
            x0, t2 = self.ignition_state(z0[0])
            if z0[2] > t2 + self.rk.t_burn2:
                continue
            try:
                sol, info, ier, msg = fsolve(
                    self.residual, np.array(z0, float), full_output=True, xtol=1e-11,
                )
                if ier == 1:
                    n = float(np.linalg.norm(self.residual(sol)))
                    key = (round(sol[0], 1), round(sol[1], 4))
                    if key not in sols or n < sols[key][1]:
                        sols[key] = (sol, n)
            except Exception:
                pass
        return sorted(sols.values(), key=lambda kv: kv[1])


def main() -> None:
    rk = RocketParams()
    solver = InscriptionSolver(rk)

    print("=" * 72)
    print("问题 2：入轨打靶（网格多初值 + fsolve）")
    print("=" * 72)
    sols = solver.solve()
    print(f"找到 {len(sols)} 组收敛根（按残差排序）：")
    for i, (sol, n) in enumerate(sols[:6]):
        print(f"  #{i+1}: t_c={sol[0]:8.3f} s  k={sol[1]:8.4f} deg/s  "
              f"t_shut={sol[2]:8.3f} s  |res|={n:.2e}")

    if not sols:
        print("未找到可行解！")
        return

    sol, n = sols[0]
    t_coast, k_deg_s, t_shut = sol

    # ---- 完整重仿验证（用 Simulator 的完整流程，含事件检测与燃尽保护）----
    x0, t2 = solver.ignition_state(t_coast)
    phi0 = np.arctan2(x0[2] * x0[0] + x0[3] * x0[1],
                      x0[2] * (-x0[1]) + x0[3] * x0[0])
    k = k_deg_s * DEG

    def controller(t):
        return phi0 + k * (t - t2)

    sim = Simulator(rk)
    res = sim.simulate(t_coast=t_coast, controller=controller, t_shut=t_shut)
    fs_ = res.final_state()

    print("-" * 72)
    print("最终验证（入轨时刻状态）：")
    print(f"  滑行时间    t_c    = {t_coast:.3f} s")
    print(f"  俯仰角速率  k      = {k_deg_s:.4f} deg/s")
    print(f"  初始俯仰角  phi0   = {phi0 / DEG:.4f} deg（与点火时刻速度方向一致）")
    print(f"  一子级关机  t1     = {res.t1:.3f} s；二子级点火 t2 = {t2:.3f} s")
    print(f"  关机时刻    t_shut = {t_shut:.3f} s（燃尽上限 {t2 + rk.t_burn2:.3f} s）")
    print(f"  轨道高度    h      = {fs_['h_km']:.4f} km（目标 400）")
    print(f"  惯性速度    V      = {fs_['v_in_mps']:.3f} m/s（目标 {rk.v_circular:.3f}）")
    print(f"  飞行路径角  gamma      = {fs_['gamma_deg']:.5f} deg（目标 0）")
    print(f"  关机质量    m      = {fs_['m_kg']:.1f} kg（剩余推进剂 "
          f"{fs_['m_kg'] - rk.ms2 - rk.m_payload:.1f} kg）")
    print(f"  入轨残差           = {orbit_residual(res.x[-1], rk)}")

    # ---- 保存 & 绘图 ----
    df = pd.DataFrame({
        "t_s": res.t, "h_km": res.h / 1e3, "V_mps": res.v_in,
        "gamma_deg": res.gamma / DEG, "m_kg": res.m, "phase": res.phase,
    })
    df.to_csv(RES_DIR / "q2_trajectory.csv", index=False)
    pd.DataFrame([{
        "t_coast_s": t_coast, "k_deg_s": k_deg_s, "t_shut_s": t_shut,
        "phi0_deg": phi0 / DEG,
        "h_final_km": fs_["h_km"], "v_final_mps": fs_["v_in_mps"],
        "gamma_final_deg": fs_["gamma_deg"],
        "m_final_kg": fs_["m_kg"],
        "propellant_remaining_kg": fs_["m_kg"] - rk.ms2 - rk.m_payload,
        "residual_norm": n,
    }]).to_csv(RES_DIR / "q2_summary.csv", index=False)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9))
    ax = axes[0]
    ax.plot(res.t / 60, res.h / 1e3, lw=1.6, color="tab:blue")
    ax.axhline(400, ls="--", color="gray", label="目标轨道 400 km")
    ax.axvline(res.t2 / 60, ls=":", color="k", alpha=0.4, label="二子级点火")
    ax.axvline(res.t_shut / 60, ls=":", color="r", alpha=0.5, label="关机")
    ax.set_xlabel("时间 t (min)"); ax.set_ylabel("高度 h (km)")
    ax.set_title("问题 2 入轨轨迹：h(t)"); ax.legend(fontsize=8); ax.grid(alpha=0.4)

    ax = axes[1]
    ax.plot(res.t / 60, res.v_in / 1e3, lw=1.6, color="tab:red")
    ax.axhline(rk.v_circular / 1e3, ls="--", color="gray",
               label=f"圆轨道速度 {rk.v_circular/1e3:.2f} km/s")
    ax.set_xlabel("时间 t (min)"); ax.set_ylabel("惯性速度 V (km/s)")
    ax.set_title("V(t)"); ax.legend(fontsize=8); ax.grid(alpha=0.4)

    ax = axes[2]
    ax.plot(res.t / 60, res.gamma / DEG, lw=1.6, color="tab:green")
    ax.axhline(0, ls="--", color="gray", label="gamma = 0")
    ax.set_xlabel("时间 t (min)"); ax.set_ylabel("飞行路径角 gamma (deg)")
    ax.set_title("gamma(t)"); ax.legend(fontsize=8); ax.grid(alpha=0.4)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "q2_orbit_inscription.png", dpi=150)
    plt.close(fig)
    print(f"已保存：{RES_DIR / 'q2_trajectory.csv'}, {RES_DIR / 'q2_summary.csv'}, "
          f"{FIG_DIR / 'q2_orbit_inscription.png'}")


if __name__ == "__main__":
    main()
