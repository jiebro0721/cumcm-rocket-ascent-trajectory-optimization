"""Regression tests for optimization variables and throttled fuel feasibility."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC_DIR))

from common import RocketParams, Simulator  # noqa: E402
from q3_fuel_opt import FuelOptimizer  # noqa: E402
from q4_throttle_opt import SIGMA_MIN, ThrottleOptimizer  # noqa: E402


def test_coast_time_perturbation_changes_ignition_state() -> None:
    """A finite-difference step must not be erased by cache-key rounding."""
    opt = FuelOptimizer(rtol=1e-10)
    x_base, _ = opt.ignition_state(103.6316115)
    x_perturbed, _ = opt.ignition_state(103.6316115 + 1.0e-8)

    assert np.linalg.norm(x_perturbed[:4] - x_base[:4]) > 0.0


def test_low_throttle_can_burn_longer_than_full_thrust_duration() -> None:
    """Throttle feasibility is governed by propellant used, not wall-clock burn time."""
    rk = RocketParams()
    opt = ThrottleOptimizer(rk=rk, n_phi=5, n_sigma=4)
    t_coast = 103.6316115
    x0, t2 = opt.ignition_state(t_coast)
    phi0_deg = np.degrees(opt.phi0_of(x0))
    burn_duration = 1.10 * rk.t_burn2

    z = np.array(
        [t_coast]
        + [phi0_deg] * (opt.n_phi - 1)
        + [SIGMA_MIN] * opt.n_sigma
        + [t2 + burn_duration],
        dtype=float,
    )
    propellant_used, _ = opt.evaluate(z)

    assert burn_duration > rk.t_burn2
    assert propellant_used < rk.mp2


def test_full_simulator_does_not_clip_throttled_burn_to_full_thrust_time() -> None:
    rk = RocketParams()
    sim = Simulator(rk)
    t_coast = 103.6316115
    requested_shutdown = rk.t_burn1 + t_coast + 1.10 * rk.t_burn2

    result = sim.simulate(
        t_coast=t_coast,
        sigma=SIGMA_MIN,
        t_shut=requested_shutdown,
    )

    assert abs(result.t[-1] - requested_shutdown) < 1.0e-6
    assert result.x[-1, 4] > rk.ms2 + rk.m_payload
