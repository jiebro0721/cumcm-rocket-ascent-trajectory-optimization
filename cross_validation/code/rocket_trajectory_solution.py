"""Auditable model and solution for the two-stage launch problem.

Numerical layers are kept separate: fixed-step RK4 for forward simulation,
a Newton predictor-corrector for Question 2, and Hermite-Simpson direct
collocation for the continuous optimal-control problems in Questions 3 and 4.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR.parent / "结果"
FIG_DIR = OUT_DIR / "figures"


@dataclass(frozen=True)
class Constants:
    re: float = 6_371_000.0
    mu: float = 3.986004418e14
    omega: float = 7.2921159e-5
    g0: float = 9.80665
    rho0: float = 1.225
    h0: float = 7_200.0
    cd: float = 0.3
    area: float = 12.5
    target_h: float = 400_000.0


C = Constants()
R_TARGET = C.re + C.target_h
V_CIRC = math.sqrt(C.mu / R_TARGET)

M0, MS1, MP1, ISP1, T1 = 500_000.0, 40_000.0, 380_000.0, 300.0, 7.5e6
MS2, MP2, M_PAYLOAD, ISP2, T2_MAX = 8_000.0, 62_000.0, 10_000.0, 420.0, 6.0e5
TB1 = MP1 * ISP1 * C.g0 / T1
TB2_MAX = MP2 * ISP2 * C.g0 / T2_MAX
M2_INITIAL = MS2 + MP2 + M_PAYLOAD
M2_DRY_PAYLOAD = MS2 + M_PAYLOAD
MDOT2_MAX = T2_MAX / (ISP2 * C.g0)
R_SCALE, V_SCALE, M_SCALE = R_TARGET, V_CIRC, M2_INITIAL


def flight_path_angle(state):
    return math.atan2(state[2], state[3])


def terminal_metrics(state):
    r, _, vr, vt, mass = state
    return {
        "height_km": (r - C.re) / 1_000.0,
        "inertial_speed_m_s": math.hypot(vr, vt),
        "ground_relative_speed_m_s": math.hypot(vr, vt - C.omega * r),
        "radial_speed_m_s": vr,
        "tangential_speed_m_s": vt,
        "flight_path_deg": math.degrees(math.atan2(vr, vt)),
        "mass_kg": mass,
    }


def rhs(_time, state, thrust, isp, theta, include_drag):
    """Variable-mass point dynamics in an Earth-centered inertial frame."""
    r, _phi, vr, vt, mass = state
    drag_r = drag_t = 0.0
    if include_drag:
        height = max(0.0, r - C.re)
        relative_t = vt - C.omega * r
        relative_speed = math.hypot(vr, relative_t)
        if relative_speed > 1e-12:
            density = C.rho0 * math.exp(-height / C.h0)
            drag = 0.5 * density * relative_speed**2 * C.cd * C.area
            drag_r = -drag * vr / relative_speed
            drag_t = -drag * relative_t / relative_speed
    return np.array(
        [
            vr,
            vt / r,
            vt**2 / r - C.mu / r**2 + (thrust * math.sin(theta) + drag_r) / mass,
            -vr * vt / r + (thrust * math.cos(theta) + drag_t) / mass,
            -thrust / (isp * C.g0) if thrust > 0 else 0.0,
        ]
    )


def rk4_segment(state0, duration, control, include_drag, n_steps):
    """Classical fourth-order Runge-Kutta formula on an equal time mesh."""
    if duration <= 0:
        return np.array([0.0]), np.array([state0], dtype=float)
    n_steps = max(1, int(n_steps))
    dt = duration / n_steps
    times = np.linspace(0.0, duration, n_steps + 1)
    states = np.empty((n_steps + 1, len(state0)))
    states[0] = state0
    for index in range(n_steps):
        time, state = times[index], states[index]

        def derivative(local_time, local_state):
            thrust, isp, theta = control(local_time, local_state)
            return rhs(local_time, local_state, thrust, isp, theta, include_drag)

        k1 = derivative(time, state)
        k2 = derivative(time + dt / 2, state + dt * k1 / 2)
        k3 = derivative(time + dt / 2, state + dt * k2 / 2)
        k4 = derivative(time + dt, state + dt * k3)
        states[index + 1] = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return times, states


def stage1_pitch(time):
    return math.pi / 2 if time <= 10 else math.pi / 2 - math.radians(0.4) * (time - 10)


def simulate_stage1(step=0.25):
    initial = np.array([C.re, 0.0, 0.0, C.omega * C.re, M0])

    def control(time, _state):
        return T1, ISP1, stage1_pitch(time)

    times, states = rk4_segment(initial, TB1, control, True, math.ceil(TB1 / step))
    states[-1, 4] = M0 - MP1
    separated = states[-1].copy()
    separated[4] -= MS1
    return times, states, separated


STAGE1_TIME, STAGE1_STATE, STAGE1_SEPARATED = simulate_stage1()


def simulate_coast(state0, duration, n_steps=320):
    return rk4_segment(state0, duration, lambda _t, _x: (0.0, ISP2, 0.0), False, n_steps)


def interpolate_control(values, tau):
    return float(np.interp(np.clip(tau, 0.0, 1.0), np.linspace(0.0, 1.0, len(values)), values))


def simulate_stage2(state0, duration, theta_function, throttle_function=None, n_steps=1200):
    throttle_function = throttle_function or (lambda _time: 1.0)

    def control(time, _state):
        throttle = float(np.clip(throttle_function(time), 0.0, 1.0))
        return throttle * T2_MAX, ISP2, float(theta_function(time))

    return rk4_segment(state0, duration, control, False, n_steps)


def combine_segments(segments):
    time_parts, state_parts, offset = [], [], 0.0
    for index, (times, states) in enumerate(segments):
        if index:
            times, states = times[1:], states[1:]
        time_parts.append(times + offset)
        state_parts.append(states)
        offset += float(times[-1]) if len(times) else 0.0
    return np.concatenate(time_parts), np.vstack(state_parts)


def q1_baseline(step=0.25):
    t1, y1, separated = simulate_stage1(step)
    tc, yc = simulate_coast(separated, 60.0, math.ceil(60.0 / step))

    def control(_time, state):
        return T2_MAX, ISP2, flight_path_angle(state)

    tb, yb = rk4_segment(yc[-1], TB2_MAX, control, False, math.ceil(TB2_MAX / step))
    times, states = combine_segments([(t1, y1), (tc, yc), (tb, yb)])
    return {
        "name": "Q1 baseline",
        "coast_time_s": 60.0,
        "burn_time_s": TB2_MAX,
        "stage2_prop_used_kg": MP2,
        "time_s": times,
        "state": states,
        "terminal": terminal_metrics(states[-1]),
    }


def terminal_residual(state):
    return np.array([(state[0] - R_TARGET) / 1_000.0, state[2] / 10.0, (state[3] - V_CIRC) / 10.0])


def q2_terminal(decision, detailed=False):
    coast_time, rate_deg_s, burn_time = decision
    _, coast = simulate_coast(STAGE1_SEPARATED, coast_time, 240)
    theta0 = flight_path_angle(coast[-1])
    theta = lambda time: theta0 + math.radians(rate_deg_s) * time
    burn_times, burn_states = simulate_stage2(coast[-1], burn_time, theta, n_steps=1600 if detailed else 700)
    return burn_times, burn_states, theta0


def solve_q2():
    """Newton predictor-corrector using an explicit terminal sensitivity matrix."""
    decision = np.array([100.0, -0.05, 398.0])
    lower, upper = np.array([0.0, -0.35, 40.0]), np.array([1_200.0, 0.20, TB2_MAX])
    increments = np.array([0.05, 2e-5, 0.05])
    records = []

    def residual_at(value):
        return terminal_residual(q2_terminal(value)[1][-1])

    for iteration in range(20):
        residual = residual_at(decision)
        jacobian = np.empty((3, 3))
        for column in range(3):
            plus, minus = decision.copy(), decision.copy()
            plus[column] += increments[column]
            minus[column] -= increments[column]
            jacobian[:, column] = (residual_at(plus) - residual_at(minus)) / (2 * increments[column])
        correction = np.linalg.solve(jacobian, -residual)
        records.append(
            {
                "iteration": iteration,
                "coast_time_s": decision[0],
                "pitch_rate_deg_s": decision[1],
                "burn_time_s": decision[2],
                "residual_norm": np.linalg.norm(residual),
                "jacobian_condition": np.linalg.cond(jacobian),
                "scaled_correction_norm": np.linalg.norm(correction / np.array([100.0, 0.05, 100.0])),
            }
        )
        if np.linalg.norm(residual) < 1e-9:
            break
        for exponent in range(12):
            candidate = np.clip(decision + (0.5**exponent) * correction, lower, upper)
            if np.linalg.norm(residual_at(candidate)) < np.linalg.norm(residual):
                decision = candidate
                break
        else:
            raise RuntimeError("Q2 Newton corrector could not reduce the terminal residual")

    _, states, theta0 = q2_terminal(decision, detailed=True)
    return {
        "name": "Q2 Newton targeting",
        "coast_time_s": decision[0],
        "rate_deg_s": decision[1],
        "burn_time_s": decision[2],
        "stage2_prop_used_kg": MDOT2_MAX * decision[2],
        "theta0_rad": theta0,
        "iterations": records,
        "terminal": terminal_metrics(states[-1]),
        "residual_norm": float(np.linalg.norm(terminal_residual(states[-1]))),
    }


def coast_initial_scaled(coast_time):
    _, states = simulate_coast(STAGE1_SEPARATED, coast_time, 160)
    return states[-1, [0, 2, 3]] / np.array([R_SCALE, V_SCALE, V_SCALE])


def vacuum_rhs_physical(state, theta, throttle):
    r, vr, vt, mass = state
    thrust = throttle * T2_MAX
    return np.array(
        [
            vr,
            vt**2 / r - C.mu / r**2 + thrust * math.sin(theta) / mass,
            -vr * vt / r + thrust * math.cos(theta) / mass,
            -thrust / (ISP2 * C.g0),
        ]
    )


def q3_scaled_rhs(state, theta, tau, burn_time):
    r, vr, vt = state * np.array([R_SCALE, V_SCALE, V_SCALE])
    mass = M2_INITIAL - MDOT2_MAX * burn_time * tau
    derivative = vacuum_rhs_physical(np.array([r, vr, vt, mass]), theta, 1.0)[:3]
    return burn_time * derivative / np.array([R_SCALE, V_SCALE, V_SCALE])


def q4_scaled_rhs(state, theta, throttle, burn_time):
    physical_state = state * np.array([R_SCALE, V_SCALE, V_SCALE, M_SCALE])
    derivative = vacuum_rhs_physical(physical_state, theta, throttle)
    return burn_time * derivative / np.array([R_SCALE, V_SCALE, V_SCALE, M_SCALE])


def hermite_simpson_q3(states, angles, burn_time):
    count, defects = len(angles), []
    step = 1.0 / (count - 1)
    for index in range(count - 1):
        tau_left, tau_right = index * step, (index + 1) * step
        f_left = q3_scaled_rhs(states[index], angles[index], tau_left, burn_time)
        f_right = q3_scaled_rhs(states[index + 1], angles[index + 1], tau_right, burn_time)
        state_mid = (states[index] + states[index + 1]) / 2 + step * (f_left - f_right) / 8
        f_mid = q3_scaled_rhs(state_mid, (angles[index] + angles[index + 1]) / 2, tau_left + step / 2, burn_time)
        defects.append(states[index + 1] - states[index] - step * (f_left + 4 * f_mid + f_right) / 6)
    return np.concatenate(defects)


def hermite_simpson_q4(states, angles, throttles, burn_time):
    count, defects = len(angles), []
    step = 1.0 / (count - 1)
    for index in range(count - 1):
        f_left = q4_scaled_rhs(states[index], angles[index], throttles[index], burn_time)
        f_right = q4_scaled_rhs(states[index + 1], angles[index + 1], throttles[index + 1], burn_time)
        state_mid = (states[index] + states[index + 1]) / 2 + step * (f_left - f_right) / 8
        f_mid = q4_scaled_rhs(
            state_mid,
            (angles[index] + angles[index + 1]) / 2,
            (throttles[index] + throttles[index + 1]) / 2,
            burn_time,
        )
        defects.append(states[index + 1] - states[index] - step * (f_left + 4 * f_mid + f_right) / 6)
    return np.concatenate(defects)


def q3_seed_from_q2(q2, nodes):
    coast_time, burn_time = q2["coast_time_s"], q2["burn_time_s"]
    _, coast = simulate_coast(STAGE1_SEPARATED, coast_time, 320)
    theta0, rate = flight_path_angle(coast[-1]), math.radians(q2["rate_deg_s"])
    theta = lambda time: theta0 + rate * time
    _, burn = simulate_stage2(coast[-1], burn_time, theta, n_steps=nodes - 1)
    states = burn[:, [0, 2, 3]] / np.array([R_SCALE, V_SCALE, V_SCALE])
    angles = np.array([theta(time) for time in np.linspace(0.0, burn_time, nodes)])
    return coast_time, burn_time, states, angles


def refine_q3_seed(solution, nodes):
    old_tau, new_tau = np.linspace(0, 1, solution["nodes"]), np.linspace(0, 1, nodes)
    states = np.column_stack(
        [np.interp(new_tau, old_tau, solution["states_scaled"][:, column]) for column in range(3)]
    )
    angles = np.interp(new_tau, old_tau, solution["angles_rad"])
    return solution["coast_time_s"], solution["burn_time_s"], states, angles


def solve_q3_collocation(q2, nodes, previous=None):
    coast0, burn0, states0, angles0 = (
        q3_seed_from_q2(q2, nodes) if previous is None else refine_q3_seed(previous, nodes)
    )
    initial = np.r_[coast0 / 300.0, burn0 / TB2_MAX, states0.ravel(), angles0]

    def unpack(vector):
        return (
            vector[0] * 300.0,
            vector[1] * TB2_MAX,
            vector[2 : 2 + 3 * nodes].reshape(nodes, 3),
            vector[2 + 3 * nodes :],
        )

    def equality(vector):
        coast_time, burn_time, states, angles = unpack(vector)
        return np.r_[
            states[0] - coast_initial_scaled(coast_time),
            hermite_simpson_q3(states, angles, burn_time),
            states[-1] - np.array([1.0, 0.0, 1.0]),
        ]

    lower = np.r_[
        0.0,
        300.0 / TB2_MAX,
        np.tile([0.93, -0.4, 0.0], nodes),
        np.full(nodes, math.radians(-45.0)),
    ]
    upper = np.r_[4.0, 1.0, np.tile([1.25, 0.5, 1.25], nodes), np.full(nodes, math.radians(90.0))]
    history = []

    def callback(vector):
        history.append(
            {
                "iteration": len(history),
                "burn_time_s": vector[1] * TB2_MAX,
                "max_equality_residual": np.max(np.abs(equality(vector))),
            }
        )

    result = minimize(
        lambda vector: vector[1],
        initial,
        method="SLSQP",
        bounds=list(zip(lower, upper)),
        constraints={"type": "eq", "fun": equality},
        callback=callback,
        options={"ftol": 1e-10, "maxiter": 800, "disp": False},
    )
    if not result.success:
        raise RuntimeError(f"Q3 collocation failed on {nodes} nodes: {result.message}")
    coast_time, burn_time, states, angles = unpack(result.x)
    return {
        "nodes": nodes,
        "coast_time_s": coast_time,
        "burn_time_s": burn_time,
        "stage2_prop_used_kg": MDOT2_MAX * burn_time,
        "states_scaled": states,
        "angles_rad": angles,
        "max_collocation_defect": float(np.max(np.abs(equality(result.x)))),
        "iterations": int(result.nit),
        "history": history,
    }


def verify_q3(solution, steps=10_000):
    coast_time, burn_time, angles = solution["coast_time_s"], solution["burn_time_s"], solution["angles_rad"]
    _, coast = simulate_coast(STAGE1_SEPARATED, coast_time, 2_000)
    theta = lambda time: interpolate_control(angles, time / burn_time)
    times, states = simulate_stage2(coast[-1], burn_time, theta, n_steps=steps)
    residual = terminal_residual(states[-1])
    return {
        "times": times,
        "states": states,
        "terminal": terminal_metrics(states[-1]),
        "scaled_terminal_residual": residual,
        "terminal_residual_norm": float(np.linalg.norm(residual)),
    }


def q4_state_angle_bounds(nodes):
    lower = np.r_[
        0.0,
        300.0 / TB2_MAX,
        np.tile([0.93, -0.4, 0.0, M2_DRY_PAYLOAD / M_SCALE], nodes),
        np.full(nodes, math.radians(-45.0)),
    ]
    upper = np.r_[
        4.0,
        1.0 / 0.6,
        np.tile([1.25, 0.5, 1.25, 1.0], nodes),
        np.full(nodes, math.radians(90.0)),
    ]
    return lower, upper


def solve_q4_fixed_throttle(q3, throttle):
    """Optimize the trajectory with one prescribed throttle level for sensitivity analysis."""
    nodes = q3["nodes"]
    burn0 = q3["burn_time_s"] / throttle
    tau = np.linspace(0.0, 1.0, nodes)
    mass = (M2_INITIAL - MDOT2_MAX * throttle * burn0 * tau) / M_SCALE
    states0 = np.column_stack([q3["states_scaled"], mass])
    initial = np.r_[q3["coast_time_s"] / 300.0, burn0 / TB2_MAX, states0.ravel(), q3["angles_rad"]]

    def unpack(vector):
        return (
            vector[0] * 300.0,
            vector[1] * TB2_MAX,
            vector[2 : 2 + 4 * nodes].reshape(nodes, 4),
            vector[2 + 4 * nodes :],
        )

    def equality(vector):
        coast_time, burn_time, states, angles = unpack(vector)
        throttles = np.full(nodes, throttle)
        return np.r_[
            states[0] - np.r_[coast_initial_scaled(coast_time), 1.0],
            hermite_simpson_q4(states, angles, throttles, burn_time),
            states[-1, :3] - np.array([1.0, 0.0, 1.0]),
        ]

    lower, upper = q4_state_angle_bounds(nodes)
    result = minimize(
        lambda vector: -unpack(vector)[2][-1, 3],
        initial,
        method="SLSQP",
        bounds=list(zip(lower, upper)),
        constraints={"type": "eq", "fun": equality},
        options={"ftol": 1e-10, "maxiter": 1_200, "disp": False},
    )
    if not result.success:
        raise RuntimeError(f"Q4 fixed-throttle solve failed at eta={throttle}: {result.message}")
    coast_time, burn_time, states, angles = unpack(result.x)
    return {
        "seed_throttle": throttle,
        "coast_time_s": coast_time,
        "burn_time_s": burn_time,
        "stage2_prop_used_kg": M2_INITIAL - states[-1, 3] * M_SCALE,
        "states_scaled": states,
        "angles_rad": angles,
        "max_collocation_defect": float(np.max(np.abs(equality(result.x)))),
        "iterations": int(result.nit),
    }


def release_q4_throttle(fixed_solution):
    """Release all throttle nodes from a feasible constant-throttle solution."""
    nodes = len(fixed_solution["angles_rad"])
    seed_throttle = fixed_solution["seed_throttle"]
    initial = np.r_[
        fixed_solution["coast_time_s"] / 300.0,
        fixed_solution["burn_time_s"] / TB2_MAX,
        fixed_solution["states_scaled"].ravel(),
        fixed_solution["angles_rad"],
        np.full(nodes, seed_throttle),
    ]

    def unpack(vector):
        return (
            vector[0] * 300.0,
            vector[1] * TB2_MAX,
            vector[2 : 2 + 4 * nodes].reshape(nodes, 4),
            vector[2 + 4 * nodes : 2 + 5 * nodes],
            vector[2 + 5 * nodes :],
        )

    def equality(vector):
        coast_time, burn_time, states, angles, throttles = unpack(vector)
        return np.r_[
            states[0] - np.r_[coast_initial_scaled(coast_time), 1.0],
            hermite_simpson_q4(states, angles, throttles, burn_time),
            states[-1, :3] - np.array([1.0, 0.0, 1.0]),
        ]

    lower, upper = q4_state_angle_bounds(nodes)
    lower = np.r_[lower, np.full(nodes, 0.6)]
    upper = np.r_[upper, np.ones(nodes)]
    result = minimize(
        lambda vector: -unpack(vector)[2][-1, 3],
        initial,
        method="SLSQP",
        bounds=list(zip(lower, upper)),
        constraints={"type": "eq", "fun": equality},
        options={"ftol": 1e-10, "maxiter": 1_500, "disp": False},
    )
    if not result.success:
        raise RuntimeError(f"Q4 released-throttle solve failed from eta={seed_throttle}: {result.message}")
    coast_time, burn_time, states, angles, throttles = unpack(result.x)
    return {
        "name": "Q4 throttle optimal control",
        "nodes": nodes,
        "seed_throttle": seed_throttle,
        "coast_time_s": coast_time,
        "burn_time_s": burn_time,
        "stage2_prop_used_kg": M2_INITIAL - states[-1, 3] * M_SCALE,
        "states_scaled": states,
        "angles_rad": angles,
        "throttles": throttles,
        "max_collocation_defect": float(np.max(np.abs(equality(result.x)))),
        "iterations": int(result.nit),
    }


def solve_q4_collocation(q3):
    fixed_solutions = [solve_q4_fixed_throttle(q3, throttle) for throttle in [1.0, 0.9, 0.8, 0.7, 0.6]]
    candidates = [release_q4_throttle(solution) for solution in fixed_solutions]
    best = min(candidates, key=lambda solution: solution["stage2_prop_used_kg"])
    best["constant_throttle_audit"] = [
        {
            "throttle": solution["seed_throttle"],
            "coast_time_s": solution["coast_time_s"],
            "burn_time_s": solution["burn_time_s"],
            "stage2_prop_used_kg": solution["stage2_prop_used_kg"],
            "max_collocation_defect": solution["max_collocation_defect"],
            "iterations": solution["iterations"],
        }
        for solution in fixed_solutions
    ]
    best["multistart_audit"] = [
        {
            "seed_throttle": solution["seed_throttle"],
            "coast_time_s": solution["coast_time_s"],
            "burn_time_s": solution["burn_time_s"],
            "stage2_prop_used_kg": solution["stage2_prop_used_kg"],
            "throttle_min": float(np.min(solution["throttles"])),
            "throttle_max": float(np.max(solution["throttles"])),
            "max_collocation_defect": solution["max_collocation_defect"],
            "iterations": solution["iterations"],
        }
        for solution in candidates
    ]
    return best


def verify_q4(solution, steps=10_000):
    coast_time, burn_time = solution["coast_time_s"], solution["burn_time_s"]
    _, coast = simulate_coast(STAGE1_SEPARATED, coast_time, 2_000)
    theta = lambda time: interpolate_control(solution["angles_rad"], time / burn_time)
    throttle = lambda time: interpolate_control(solution["throttles"], time / burn_time)
    times, states = simulate_stage2(coast[-1], burn_time, theta, throttle, steps)
    residual = terminal_residual(states[-1])
    return {
        "times": times,
        "states": states,
        "terminal": terminal_metrics(states[-1]),
        "scaled_terminal_residual": residual,
        "terminal_residual_norm": float(np.linalg.norm(residual)),
    }


def build_full_trajectory(coast_time, burn_time, theta, throttle=None, stage2_steps=2_000):
    coast_times, coast = simulate_coast(STAGE1_SEPARATED, coast_time, 500)
    burn_times, burn = simulate_stage2(coast[-1], burn_time, theta, throttle, stage2_steps)
    return combine_segments([(STAGE1_TIME, STAGE1_STATE), (coast_times, coast), (burn_times, burn)])


def trajectory_frame(times, states):
    return pd.DataFrame(
        {
            "time_s": times,
            "height_km": (states[:, 0] - C.re) / 1_000.0,
            "inertial_speed_m_s": np.hypot(states[:, 2], states[:, 3]),
            "ground_relative_speed_m_s": np.hypot(states[:, 2], states[:, 3] - C.omega * states[:, 0]),
            "flight_path_deg": np.degrees(np.arctan2(states[:, 2], states[:, 3])),
            "mass_kg": states[:, 4],
            "radial_speed_m_s": states[:, 2],
            "tangential_speed_m_s": states[:, 3],
        }
    )


def rk4_convergence_table():
    rows = []
    for step in [1.0, 0.5, 0.25]:
        terminal = q1_baseline(step)["terminal"]
        rows.append(
            {
                "step_s": step,
                "final_height_km": terminal["height_km"],
                "final_speed_m_s": terminal["inertial_speed_m_s"],
                "final_flight_path_deg": terminal["flight_path_deg"],
                "final_mass_kg": terminal["mass_kg"],
            }
        )
    return pd.DataFrame(rows)


def loss_budget(times, states, theta_function, throttle_function=None):
    throttle_function = throttle_function or (lambda _time: 1.0)
    theta = np.array([theta_function(time) for time in times])
    throttle = np.array([throttle_function(time) for time in times])
    gamma = np.arctan2(states[:, 2], states[:, 3])
    thrust_acceleration = throttle * T2_MAX / states[:, 4]
    ideal = float(np.trapz(thrust_acceleration, times))
    steering = float(np.trapz(thrust_acceleration * (1 - np.cos(theta - gamma)), times))
    gravity = float(np.trapz(C.mu / states[:, 0] ** 2 * np.sin(gamma), times))
    speed_gain = float(
        math.hypot(states[-1, 2], states[-1, 3]) - math.hypot(states[0, 2], states[0, 3])
    )
    return {
        "ideal_delta_v_m_s": ideal,
        "actual_speed_gain_m_s": speed_gain,
        "gravity_loss_m_s": gravity,
        "steering_loss_m_s": steering,
        "balance_error_m_s": ideal - speed_gain - gravity - steering,
    }


def plot_results(q1, trajectories, angle_series, q4_tau, q4_throttle):
    frame = trajectory_frame(q1["time_s"], q1["state"])
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), dpi=180)
    columns = ["height_km", "inertial_speed_m_s", "flight_path_deg", "mass_kg"]
    labels = ["height / km", "inertial speed / m s$^{-1}$", "flight-path angle / deg", "mass / kg"]
    for axis, column, label in zip(axes.ravel(), columns, labels):
        axis.plot(frame["time_s"], frame[column])
        axis.set(xlabel="time / s", ylabel=label)
        axis.grid(True, alpha=0.3)
    fig.suptitle("Q1 baseline trajectory: fixed-step RK4")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "q1_baseline_history.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), dpi=180)
    for label, (times, states) in trajectories.items():
        data = trajectory_frame(times, states)
        axes[0].plot(data["time_s"], data["height_km"], label=label)
        axes[1].plot(data["time_s"], data["inertial_speed_m_s"], label=label)
        axes[2].plot(data["time_s"], data["flight_path_deg"], label=label)
    axes[0].axhline(400.0, color="k", lw=0.8, ls="--")
    axes[1].axhline(V_CIRC, color="k", lw=0.8, ls="--")
    axes[2].axhline(0.0, color="k", lw=0.8, ls="--")
    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=8)
    axes[0].set_ylabel("height / km")
    axes[1].set_ylabel("speed / m s$^{-1}$")
    axes[2].set(ylabel="flight-path angle / deg", xlabel="time / s")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "trajectory_comparison.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7), dpi=180)
    for label, tau, angle in angle_series:
        axes[0].plot(tau, angle, label=label)
    axes[0].set_ylabel("pitch angle / deg")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].plot(q4_tau, q4_throttle)
    axes[1].set(xlabel="normalized second-stage burn time", ylabel="throttle ratio", ylim=(0.55, 1.05))
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "optimized_controls.png", bbox_inches="tight")
    plt.close(fig)


def json_safe(value):
    if isinstance(value, dict):
        return {
            key: json_safe(item)
            for key, item in value.items()
            if key not in {"states_scaled", "angles_rad", "times", "states", "time_s", "state"}
        }
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    q1, q2 = q1_baseline(0.25), solve_q2()
    mesh_solutions, previous = [], None
    for nodes in [11, 21, 31]:
        solution = solve_q3_collocation(q2, nodes, previous)
        solution["verification"] = verify_q3(solution, 4_000 if nodes < 31 else 10_000)
        mesh_solutions.append(solution)
        previous = solution
    q3 = mesh_solutions[-1]
    q4 = solve_q4_collocation(q3)
    q4["verification"] = verify_q4(q4, 10_000)

    q2_theta = lambda time: q2["theta0_rad"] + math.radians(q2["rate_deg_s"]) * time
    q3_theta = lambda time: interpolate_control(q3["angles_rad"], time / q3["burn_time_s"])
    q4_theta = lambda time: interpolate_control(q4["angles_rad"], time / q4["burn_time_s"])
    q4_throttle = lambda time: interpolate_control(q4["throttles"], time / q4["burn_time_s"])
    trajectories = {
        "Q1": (q1["time_s"], q1["state"]),
        "Q2": build_full_trajectory(q2["coast_time_s"], q2["burn_time_s"], q2_theta),
        "Q3": build_full_trajectory(q3["coast_time_s"], q3["burn_time_s"], q3_theta),
        "Q4": build_full_trajectory(q4["coast_time_s"], q4["burn_time_s"], q4_theta, q4_throttle),
    }

    summary_rows = [
        {"case": "Q1", "coast_time_s": 60.0, "second_stage_burn_s": TB2_MAX, "stage2_prop_used_kg": MP2, **q1["terminal"]},
        {"case": "Q2", "coast_time_s": q2["coast_time_s"], "second_stage_burn_s": q2["burn_time_s"], "stage2_prop_used_kg": q2["stage2_prop_used_kg"], **q2["terminal"]},
        {"case": "Q3", "coast_time_s": q3["coast_time_s"], "second_stage_burn_s": q3["burn_time_s"], "stage2_prop_used_kg": q3["stage2_prop_used_kg"], **q3["verification"]["terminal"]},
        {"case": "Q4", "coast_time_s": q4["coast_time_s"], "second_stage_burn_s": q4["burn_time_s"], "stage2_prop_used_kg": q4["stage2_prop_used_kg"], **q4["verification"]["terminal"]},
    ]
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "summary.csv", index=False)
    for label, (times, states) in trajectories.items():
        trajectory_frame(times, states).to_csv(OUT_DIR / f"{label.lower()}_trajectory.csv", index=False)

    phase_rows = [
        {"case": "COMMON", "phase": "lift_off", "time_s": 0.0, **terminal_metrics(STAGE1_STATE[0])},
        {"case": "COMMON", "phase": "stage1_burnout_before_separation", "time_s": TB1, **terminal_metrics(STAGE1_STATE[-1])},
        {"case": "COMMON", "phase": "stage1_separation_after_jettison", "time_s": TB1, **terminal_metrics(STAGE1_SEPARATED)},
    ]
    for label, coast_time, burn_time, theta, throttle in [
        ("Q2", q2["coast_time_s"], q2["burn_time_s"], q2_theta, None),
        ("Q3", q3["coast_time_s"], q3["burn_time_s"], q3_theta, None),
        ("Q4", q4["coast_time_s"], q4["burn_time_s"], q4_theta, q4_throttle),
    ]:
        _, coast = simulate_coast(STAGE1_SEPARATED, coast_time, 500)
        _, burn = simulate_stage2(coast[-1], burn_time, theta, throttle, 3_000)
        phase_rows.append({"case": label, "phase": "second_stage_ignition", "time_s": TB1 + coast_time, **terminal_metrics(coast[-1])})
        phase_rows.append({"case": label, "phase": "terminal_shutdown", "time_s": TB1 + coast_time + burn_time, **terminal_metrics(burn[-1])})
    pd.DataFrame(phase_rows).to_csv(OUT_DIR / "intermediate_phase_states.csv", index=False)

    pd.DataFrame(q2["iterations"]).to_csv(OUT_DIR / "q2_newton_iterations.csv", index=False)
    pd.DataFrame(
        [
            {
                "nodes": solution["nodes"],
                "coast_time_s": solution["coast_time_s"],
                "burn_time_s": solution["burn_time_s"],
                "propellant_kg": solution["stage2_prop_used_kg"],
                "max_collocation_defect": solution["max_collocation_defect"],
                "independent_terminal_residual_norm": solution["verification"]["terminal_residual_norm"],
                "optimizer_iterations": solution["iterations"],
            }
            for solution in mesh_solutions
        ]
    ).to_csv(OUT_DIR / "q3_mesh_convergence.csv", index=False)

    q3_tau, q4_tau = np.linspace(0, 1, q3["nodes"]), np.linspace(0, 1, q4["nodes"])
    pd.DataFrame(
        {
            "node": np.arange(q3["nodes"]),
            "tau": q3_tau,
            "radius_m": q3["states_scaled"][:, 0] * R_SCALE,
            "radial_speed_m_s": q3["states_scaled"][:, 1] * V_SCALE,
            "tangential_speed_m_s": q3["states_scaled"][:, 2] * V_SCALE,
            "pitch_deg": np.degrees(q3["angles_rad"]),
        }
    ).to_csv(OUT_DIR / "q3_collocation_nodes.csv", index=False)
    pd.DataFrame(
        {
            "node": np.arange(q4["nodes"]),
            "tau": q4_tau,
            "radius_m": q4["states_scaled"][:, 0] * R_SCALE,
            "radial_speed_m_s": q4["states_scaled"][:, 1] * V_SCALE,
            "tangential_speed_m_s": q4["states_scaled"][:, 2] * V_SCALE,
            "mass_kg": q4["states_scaled"][:, 3] * M_SCALE,
            "pitch_deg": np.degrees(q4["angles_rad"]),
            "throttle": q4["throttles"],
        }
    ).to_csv(OUT_DIR / "q4_collocation_nodes.csv", index=False)
    pd.DataFrame(q4["constant_throttle_audit"]).to_csv(
        OUT_DIR / "q4_constant_throttle_sensitivity.csv", index=False
    )
    pd.DataFrame(q4["multistart_audit"]).to_csv(OUT_DIR / "q4_multistart_audit.csv", index=False)
    rk4_convergence_table().to_csv(OUT_DIR / "rk4_step_convergence.csv", index=False)

    loss_rows = []
    for label, coast_time, burn_time, theta, throttle in [
        ("Q2", q2["coast_time_s"], q2["burn_time_s"], q2_theta, None),
        ("Q3", q3["coast_time_s"], q3["burn_time_s"], q3_theta, None),
        ("Q4", q4["coast_time_s"], q4["burn_time_s"], q4_theta, q4_throttle),
    ]:
        _, coast = simulate_coast(STAGE1_SEPARATED, coast_time, 1_000)
        burn_times, burn = simulate_stage2(coast[-1], burn_time, theta, throttle, 6_000)
        loss_rows.append({"case": label, **loss_budget(burn_times, burn, theta, throttle)})
    pd.DataFrame(loss_rows).to_csv(OUT_DIR / "second_stage_loss_budget.csv", index=False)

    angle_series = [
        ("Q2 constant rate", np.linspace(0, 1, 200), np.degrees([q2_theta(t) for t in np.linspace(0, q2["burn_time_s"], 200)])),
        ("Q3 collocation", q3_tau, np.degrees(q3["angles_rad"])),
        ("Q4 collocation", q4_tau, np.degrees(q4["angles_rad"])),
    ]
    plot_results(q1, trajectories, angle_series, q4_tau, q4["throttles"])

    payload = {
        "constants": {**asdict(C), "target_radius_m": R_TARGET, "circular_speed_m_s": V_CIRC, "tb1_s": TB1, "tb2_max_s": TB2_MAX},
        "methods": {
            "forward_integration": "fixed-step classical RK4",
            "q2": "Newton predictor-corrector with a finite-difference sensitivity matrix and line search",
            "q3_q4": "Hermite-Simpson direct collocation followed by independent RK4 verification",
        },
        "Q1": json_safe(q1),
        "Q2": json_safe(q2),
        "Q3": json_safe(q3),
        "Q4": json_safe(q4),
    }
    (OUT_DIR / "results.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print("\nQ3 mesh convergence")
    print(pd.read_csv(OUT_DIR / "q3_mesh_convergence.csv").to_string(index=False))


if __name__ == "__main__":
    main()
