"""Regression tests for model constraints (collocation-based pipeline).

These tests guard the two root-cause fixes that motivated the direct
collocation rewrite:

1. Coast-time cache keys must preserve finite-difference perturbations
   (SLSQP t_c derivatives previously vanished due to `round(tc, 6)`).
2. Throttled burn duration may exceed the full-thrust duration; capacity is
   governed by propellant used, not wall-clock time.

Legacy SLSQP implementations have been archived under `legacy/shooting_baseline/`
for comparison; the canonical pipeline is `src/q34_direct_collocation.py`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC_DIR))

from common import RocketParams, Simulator, orbit_residual  # noqa: E402


def test_coast_time_perturbation_changes_ignition_state() -> None:
    """A finite-difference step must not be erased by cache-key rounding."""
    rk = RocketParams()
    t_ign = rk.t_burn1 + 103.6316115
    base = Simulator(rk).simulate(t_coast=103.6316115, controller=None,
                                  t_shut=t_ign + 1e-6)
    x_base = base.stages[-1].x[0, :4]
    t_ign2 = rk.t_burn1 + 103.6316115 + 1e-8
    pert = Simulator(rk).simulate(t_coast=103.6316115 + 1e-8, controller=None,
                                  t_shut=t_ign2 + 1e-6)
    x_pert = pert.stages[-1].x[0, :4]

    assert np.linalg.norm(x_pert - x_base) > 0.0


def test_low_throttle_can_burn_longer_than_full_thrust_duration() -> None:
    """Throttle feasibility is governed by propellant used, not wall-clock burn time."""
    rk = RocketParams()
    sim = Simulator(rk)
    t_coast = 103.6316115
    burn_duration = 1.10 * rk.t_burn2
    t_ign = rk.t_burn1 + t_coast

    # sigma = 0.6 且燃烧 1.10 倍满推力时长：不应被截断至满推力燃尽
    result = sim.simulate(
        t_coast=t_coast,
        sigma=0.6,
        t_shut=t_ign + burn_duration,
    )

    assert burn_duration > rk.t_burn2
    assert abs(result.t[-1] - (t_ign + burn_duration)) < 1.0e-6
    # 60% 推力燃烧 1.10*425.61s 耗药 = 0.6*1.10*62000 = 40920 kg < 62000
    prop_used = rk.m_after_sep - result.x[-1, 4]
    assert 0.6 * 1.10 * rk.mp2 - 1.0 < prop_used < 0.6 * 1.10 * rk.mp2 + 1.0


def test_full_simulator_does_not_clip_throttled_burn_to_full_thrust_time() -> None:
    rk = RocketParams()
    sim = Simulator(rk)
    t_coast = 103.6316115
    requested_shutdown = rk.t_burn1 + t_coast + 1.10 * rk.t_burn2

    result = sim.simulate(
        t_coast=t_coast,
        sigma=0.8,
        t_shut=requested_shutdown,
    )

    assert abs(result.t[-1] - requested_shutdown) < 1.0e-6
    assert result.x[-1, 4] > rk.ms2 + rk.m_payload


def test_orbit_residual_excludes_retrograde_orbit() -> None:
    """顺行圆轨道残差为零；逆行圆轨道（|v|=vc、vr=0、vt=-vc）必须被排除。"""
    rk = RocketParams()
    a = rk.a_target
    vc = rk.v_circular

    prograde = np.array([a, 0.0, 0.0, vc, 18000.0])
    res_pro = orbit_residual(prograde, rk)
    assert np.linalg.norm(res_pro) < 1e-12

    retrograde = np.array([a, 0.0, 0.0, -vc, 18000.0])
    res_retro = orbit_residual(retrograde, rk)
    # 半径与径向速度为零，仅切向速度残差为 -2（顺行方向以东为正）
    assert abs(res_retro[0]) < 1e-12
    assert abs(res_retro[1]) < 1e-12
    assert abs(res_retro[2] + 2.0) < 1e-12
