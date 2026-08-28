"""问题 1：基准控制策略下的全程动力学仿真。

基准策略（题目给定）：
- 一子级：额定最大推力；垂直起飞 10 s 后俯仰角以 0.4deg/s 线性减小（程序转弯），
  直至推进剂耗尽关机分离；
- 二子级：全额推力，推力方向始终与（相对大气）速度方向一致（恒攻角为零）；
- 滑行时间 60 s。
输出：起飞至二子级燃料耗尽全程的高度、速度、飞行路径角、总质量曲线，
以及二子级关机时刻的轨道高度、惯性速度与飞行路径角。
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    DEG,
    R_E,
    RES_DIR,
    FIG_DIR,
    RocketParams,
    Simulator,
    orbit_residual,
)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial"]
plt.rcParams["axes.unicode_minus"] = False


def phase_boundaries(res) -> dict:
    """各阶段边界时刻（供绘图竖线）。"""
    return {
        "垂直段结束": res.stages[0].t[-1],
        "一子级关机": res.t1,
        "二子级点火": res.t2,
        "二子级关机": res.t_shut,
    }


def plot_profiles(res: SimResult, rk: RocketParams, fname: str) -> None:
    """四联图：高度、惯性速度、飞行路径角、总质量（含阶段边界竖线）。"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(res.t / 60.0, res.h / 1e3, lw=1.6, color="tab:blue")
    ax.set_xlabel("时间 t (min)")
    ax.set_ylabel("高度 h (km)")
    ax.set_title("(a) 飞行高度")
    ax.grid(alpha=0.4)

    ax = axes[0, 1]
    ax.plot(res.t / 60.0, res.v_in / 1e3, lw=1.6, color="tab:red")
    ax.axhline(rk.v_circular / 1e3, ls="--", lw=1.0, color="gray",
               label=f"目标圆轨道速度 {rk.v_circular/1e3:.1f} km/s")
    ax.set_xlabel("时间 t (min)")
    ax.set_ylabel("惯性速度 V (km/s)")
    ax.set_title("(b) 惯性速度")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.4)

    ax = axes[1, 0]
    ax.plot(res.t / 60.0, res.gamma / DEG, lw=1.6, color="tab:green", label="惯性 gamma")
    ax.plot(res.t / 60.0, res.gamma_rel / DEG, lw=1.2, ls="--", color="tab:orange",
            label="相对大气 gamma")
    ax.set_xlabel("时间 t (min)")
    ax.set_ylabel("飞行路径角 gamma (deg)")
    ax.set_title("(c) 飞行路径角")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.4)

    ax = axes[1, 1]
    ax.plot(res.t / 60.0, res.m / 1e3, lw=1.6, color="tab:purple")
    ax.set_xlabel("时间 t (min)")
    ax.set_ylabel("总质量 m (t)")
    ax.set_title("(d) 总质量")
    ax.grid(alpha=0.4)

    for ax in axes.ravel():
        for name, tb in phase_boundaries(res).items():
            if tb is not None:
                ax.axvline(tb / 60.0, ls=":", lw=0.8, color="k", alpha=0.5)
    fig.suptitle("问题 1：基准控制策略全程飞行曲线", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG_DIR / fname, dpi=150)
    plt.close(fig)


def plot_trajectory(res: SimResult, fname: str) -> None:
    """地面轨迹（经度 vs 高度）与轨道图。"""
    rx, ry = res.x[:, 0], res.x[:, 1]
    lon = np.arctan2(ry, rx) / DEG
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    ax.plot(lon, res.h / 1e3, lw=1.5, color="tab:blue")
    # 阶段着色（用散点区分）
    for name, color in [("vertical", "tab:blue"), ("pitch-over", "tab:green"),
                        ("coast", "tab:orange"), ("stage2", "tab:red")]:
        mask = np.array([p == name for p in res.phase])
        if mask.any():
            ax.plot(lon[mask], res.h[mask] / 1e3, lw=1.8, color=color, label=name)
    ax.set_xlabel("经度 (deg)")
    ax.set_ylabel("高度 h (km)")
    ax.set_title("发射轨迹（经度-高度）")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.4)

    ax = axes[1]
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(R_E * np.cos(theta) / 1e3, R_E * np.sin(theta) / 1e3,
            lw=1.0, color="gray", label="地球表面")
    ax.plot([0, 0], [-R_E / 1e3, R_E / 1e3], color="k", lw=0.5)
    ax.plot(rx / 1e3, ry / 1e3, lw=1.8, color="tab:red", label="火箭轨迹")
    ax.plot(rx[0] / 1e3, ry[0] / 1e3, "ko", ms=5, label="发射点")
    ax.plot(rx[-1] / 1e3, ry[-1] / 1e3, "k*", ms=12, label="关机点")
    ax.set_aspect("equal")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_title("赤道平面内轨迹")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG_DIR / fname, dpi=150)
    plt.close(fig)


def main() -> None:
    rk = RocketParams()
    sim = Simulator(rk)
    res = sim.simulate(t_coast=rk.t_coast_q1, controller=None, t_shut=None)

    print("=" * 72)
    print("问题 1：基准控制策略仿真结果")
    print("=" * 72)
    print(f"一子级质量流量        mdot1 = {rk.mdot1:10.2f} kg/s")
    print(f"一子级推进剂耗尽时刻  t1    = {rk.t_burn1:10.2f} s")
    print(f"二子级质量流量        mdot2 = {rk.mdot2:10.2f} kg/s")
    print(f"二子级推进剂耗尽时长  tb2   = {rk.t_burn2:10.2f} s")
    print(f"分离时一子级关机质量        = {rk.M0 - rk.mp1:10.0f} kg")
    print(f"二子级点火质量（分离后）    = {rk.m_after_sep:10.0f} kg")
    print("-" * 72)
    print("各阶段边界时刻：")
    for name, tb in phase_boundaries(res).items():
        print(f"  {name:<8s} t = {tb:8.2f} s" if tb is not None else f"  {name:<8s} --")
    print("-" * 72)
    fs = res.final_state()
    print("二子级关机时刻状态（最终）：")
    print(f"  轨道高度 h     = {fs['h_km']:10.2f} km    (目标 400 km)")
    print(f"  惯性速度 V     = {fs['v_in_mps']:10.2f} m/s  (圆轨道 {rk.v_circular:.2f} m/s)")
    print(f"  飞行路径角 gamma   = {fs['gamma_deg']:10.4f} deg (惯性) / "
          f"{fs['gamma_rel_deg']:.4f} deg (相对大气)")
    print(f"  总质量 m       = {fs['m_kg']:10.0f} kg   (结构+载荷 = {rk.ms2 + rk.m_payload:.0f} kg)")
    res_orb = orbit_residual(res.x[-1], rk)
    print(f"  入轨条件残差   = {res_orb}  (0 为精确入轨)")
    print("-" * 72)
    print("说明：基准策略下二子级燃料耗尽即关机，入轨条件不一定满足；")
    print("      超出/不足的差异即问题 2 需要修正的内容。")

    # ---- 保存结果 ----
    df = pd.DataFrame({
        "t_s": res.t,
        "h_km": res.h / 1e3,
        "V_inertial_mps": res.v_in,
        "gamma_inertial_deg": res.gamma / DEG,
        "gamma_rel_deg": res.gamma_rel / DEG,
        "m_kg": res.m,
        "T_N": res.T,
        "phase": res.phase,
    })
    df.to_csv(RES_DIR / "q1_full_trajectory.csv", index=False)

    summary = {
        "t1_burn1_s": rk.t_burn1,
        "t2_ignite_s": res.t2,
        "t_shut_s": res.t_shut,
        "h_final_km": fs["h_km"],
        "v_final_mps": fs["v_in_mps"],
        "gamma_final_deg": fs["gamma_deg"],
        "gamma_rel_final_deg": fs["gamma_rel_deg"],
        "m_final_kg": fs["m_kg"],
        "h_target_km": rk.H_target / 1e3,
        "v_circular_mps": rk.v_circular,
    }
    pd.DataFrame([summary]).to_csv(RES_DIR / "q1_summary.csv", index=False)

    plot_profiles(res, rk, "q1_flight_profiles.png")
    plot_trajectory(res, "q1_trajectory.png")
    print(f"已保存：{RES_DIR / 'q1_full_trajectory.csv'}, {RES_DIR / 'q1_summary.csv'}")
    print(f"已保存：{FIG_DIR / 'q1_flight_profiles.png'}, {FIG_DIR / 'q1_trajectory.png'}")


if __name__ == "__main__":
    main()
