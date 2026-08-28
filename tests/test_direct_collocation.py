"""Small mathematical checks for the direct-collocation transcription."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC_DIR))

from q34_direct_collocation import (  # noqa: E402
    DirectCollocationSolver,
    hermite_simpson_defect_numeric,
    inertial_flight_path_angle_numeric,
    linear_control_integral,
)


def test_hermite_simpson_defect_is_zero_for_constant_derivative() -> None:
    x0 = np.array([1.0, -2.0])
    derivative = np.array([3.0, 4.0])
    step = 0.25
    x1 = x0 + step * derivative

    defect = hermite_simpson_defect_numeric(
        x0, x1, derivative, derivative, derivative, step
    )

    np.testing.assert_allclose(defect, np.zeros(2), atol=1e-14)


def test_linear_control_integral_matches_trapezoidal_rule() -> None:
    sigma = np.array([0.6, 0.8, 1.0])

    integral = linear_control_integral(sigma, duration=100.0)

    assert abs(integral - 80.0) < 1e-12


def test_solver_uses_repository_parameter_interface() -> None:
    solver = DirectCollocationSolver(question=3, n_coast=2, n_burn=2)

    assert solver.r_scale == solver.rk.a_target
    assert solver.m_scale == solver.rk.m_after_sep


def test_flight_path_angle_uses_local_radial_and_prograde_components() -> None:
    state = np.array([2.0, 0.0, 1.0, 1.0, 1.0])

    angle = inertial_flight_path_angle_numeric(state)

    assert abs(angle - np.pi / 4.0) < 1e-14
