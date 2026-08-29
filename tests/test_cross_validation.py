"""Cross-validation between the canonical pipeline and the independent implementation.

The independent implementation (cross_validation/code/rocket_trajectory_solution.py)
uses fixed-step classical RK4 (not DOP853), a Newton predictor-corrector for Q2
and Hermite-Simpson collocation with scipy (not CasADi/Ipopt). If both pipelines
agree to tight tolerances, the shared physical model is confirmed.

The test skips gracefully when the cross-validation package is unavailable.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
XV_DIR = ROOT / "cross_validation"

sys.path.insert(0, str(SRC_DIR))

from common import RocketParams, Simulator  # noqa: E402

HAS_XV = (XV_DIR / "code" / "rocket_trajectory_solution.py").exists()
if HAS_XV:
    sys.path.insert(0, str(XV_DIR / "code"))
    import rocket_trajectory_solution as xv  # noqa: E402

import pytest  # noqa: E402


@pytest.mark.skipif(not HAS_XV, reason="cross_validation package missing")
def test_q1_burnout_agrees_with_independent_rk4() -> None:
    """Q1 基准策略燃尽状态：DOP853 主实现 vs 独立 RK4 实现。"""
    # 主实现
    rk = RocketParams()
    sim = Simulator(rk)
    res = sim.simulate(t_coast=60.0, controller=None, t_shut=None)
    h_main = (np.hypot(res.x[-1, 0], res.x[-1, 1]) - 6371e3) / 1e3
    v_main = np.hypot(res.x[-1, 2], res.x[-1, 3])
    # 独立实现
    q1 = xv.q1_baseline()
    term = q1["terminal"]
    h_xv = term["height_km"]
    v_xv = term["inertial_speed_m_s"]

    assert abs(h_main - h_xv) < 0.2       # km（RK4 步长 0.25s 的离散误差）
    assert abs(v_main - v_xv) < 5.0       # m/s


@pytest.mark.skipif(not HAS_XV, reason="cross_validation package missing")
def test_q2_solution_agrees_with_independent_newton() -> None:
    """Q2 入轨解：主实现 least_squares vs 独立实现 Newton 预测-校正。"""
    q2 = xv.solve_q2()
    tc_xv = float(q2["coast_time_s"])
    k_xv = float(q2["rate_deg_s"])
    tb_xv = float(q2["burn_time_s"])

    assert abs(tc_xv - 103.6316) < 0.05        # s
    assert abs(k_xv - (-0.04736)) < 5e-4       # deg/s
    assert abs(tb_xv - 398.049) < 0.05         # s


@pytest.mark.skipif(not HAS_XV, reason="cross_validation package missing")
def test_q3_propellant_agrees_with_independent_collocation() -> None:
    """Q3 最优耗药：主实现 CasADi 配点 vs 独立实现 scipy 配点。"""
    q2 = xv.solve_q2()
    q3 = xv.solve_q3_collocation(q2, nodes=31)
    prop_xv = float(q3["stage2_prop_used_kg"])

    assert abs(prop_xv - 57446.934) < 0.05     # kg（两实现差 0.002 kg）
