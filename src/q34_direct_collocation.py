"""Questions 3 and 4 solved by Hermite--Simpson direct collocation.

The nonlinear program exposes states, controls, coast duration, and burn duration
as decision variables.  DOP853 is used only after optimization to independently
re-integrate the resulting control law.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import casadi as ca
import numpy as np
import pandas as pd

from common import DEG, G0, MU, RES_DIR, RocketParams, Simulator, orbit_residual, thrust_unit


SIGMA_MIN = 0.6


def inertial_flight_path_angle_numeric(state: np.ndarray) -> float:
    """Inertial flight-path angle measured from the local prograde tangent."""
    x, y, vx, vy = np.asarray(state, dtype=float)[:4]
    return float(np.arctan2(vx * x + vy * y, vx * (-y) + vy * x))


def hermite_simpson_defect_numeric(x0, x1, f0, f1, fmid, step):
    """Return the Hermite--Simpson interval defect for numerical checks."""
    return np.asarray(x1) - np.asarray(x0) - step * (
        np.asarray(f0) + 4.0 * np.asarray(fmid) + np.asarray(f1)
    ) / 6.0


def linear_control_integral(nodes: np.ndarray, duration: float) -> float:
    """Integrate equally spaced, linearly interpolated control nodes."""
    values = np.asarray(nodes, dtype=float)
    if values.size < 2:
        return float(values[0] * duration)
    return float(np.trapezoid(values, dx=duration / (values.size - 1)))


@dataclass
class CollocationSolution:
    question: int
    n_coast: int
    n_burn: int
    coast_duration: float
    burn_duration: float
    coast_states: np.ndarray
    burn_states: np.ndarray
    phi_deg: np.ndarray
    sigma: np.ndarray
    propellant_used: float
    objective: float
    max_defect: float


class DirectCollocationSolver:
    """Nondimensional Hermite--Simpson transcription for Q3 or Q4."""

    def __init__(self, question: int, n_coast: int = 10, n_burn: int = 20):
        if question not in (3, 4):
            raise ValueError("question must be 3 or 4")
        self.question = question
        self.n_coast = n_coast
        self.n_burn = n_burn
        self.rk = RocketParams()
        self.r_scale = self.rk.a_target
        self.v_scale = self.rk.v_circular
        self.m_scale = self.rk.m_after_sep

    def scale_state(self, x: np.ndarray) -> np.ndarray:
        return np.array([
            x[0] / self.r_scale,
            x[1] / self.r_scale,
            x[2] / self.v_scale,
            x[3] / self.v_scale,
            x[4] / self.m_scale,
        ])

    def unscale_states(self, x: np.ndarray) -> np.ndarray:
        scales = np.array([self.r_scale, self.r_scale, self.v_scale,
                           self.v_scale, self.m_scale])[:, None]
        return np.asarray(x) * scales

    def _dynamics(self, x, phi, sigma):
        rx, ry, vx, vy, mass = ca.vertsplit(x)
        radius = ca.sqrt(rx * rx + ry * ry)
        tx = ca.cos(phi) * (-ry / radius) + ca.sin(phi) * (rx / radius)
        ty = ca.cos(phi) * (rx / radius) + ca.sin(phi) * (ry / radius)

        gravity_scale = MU / (self.r_scale**2 * self.v_scale)
        thrust_scale = self.rk.Fmax2 / (self.m_scale * self.v_scale)
        return ca.vertcat(
            self.v_scale / self.r_scale * vx,
            self.v_scale / self.r_scale * vy,
            -gravity_scale * rx / radius**3 + thrust_scale * sigma * tx / mass,
            -gravity_scale * ry / radius**3 + thrust_scale * sigma * ty / mass,
            -sigma * self.rk.Fmax2 / (self.rk.Isp2 * G0 * self.m_scale),
        )

    def _coast_dynamics(self, x):
        rx, ry, vx, vy, _mass = ca.vertsplit(x)
        radius = ca.sqrt(rx * rx + ry * ry)
        gravity_scale = MU / (self.r_scale**2 * self.v_scale)
        return ca.vertcat(
            self.v_scale / self.r_scale * vx,
            self.v_scale / self.r_scale * vy,
            -gravity_scale * rx / radius**3,
            -gravity_scale * ry / radius**3,
            0.0,
        )

    @staticmethod
    def _hs_constraint(opti, x0, x1, f0, f1, fmid_fn, step):
        xmid = 0.5 * (x0 + x1) + step * (f0 - f1) / 8.0
        fmid = fmid_fn(xmid)
        opti.subject_to(x1 - x0 - step * (f0 + 4.0 * fmid + f1) / 6.0 == 0)

    def _reference_trajectory(self):
        q2 = pd.read_csv(RES_DIR / "q2_summary.csv").iloc[0]
        tc = float(q2["t_coast_s"])
        shut = float(q2["t_shut_s"])
        k = float(q2["k_deg_s"]) * DEG

        probe = Simulator(self.rk).simulate(
            t_coast=tc, controller=None,
            t_shut=self.rk.t_burn1 + tc + 1e-6,
        )
        ignition = probe.stages[-1].x[0]
        t2 = self.rk.t_burn1 + tc
        phi0 = np.arctan2(
            ignition[2] * ignition[0] + ignition[3] * ignition[1],
            ignition[2] * (-ignition[1]) + ignition[3] * ignition[0],
        )

        def controller(t):
            return phi0 + k * (t - t2)

        result = Simulator(self.rk).simulate(
            t_coast=tc, controller=controller, t_shut=shut, rtol=1e-10,
        )
        return tc, shut - t2, phi0, k, result

    @staticmethod
    def _sample_state(result, query_times):
        return np.vstack([
            np.interp(query_times, result.t, result.x[:, j])
            for j in range(result.x.shape[1])
        ])

    def solve(self, warm: CollocationSolution | None = None, print_level: int = 0):
        opti = ca.Opti()
        nc, nb = self.n_coast, self.n_burn
        tc = opti.variable()
        tb = opti.variable()
        xc = opti.variable(5, nc + 1)
        xb = opti.variable(5, nb + 1)
        phi = opti.variable(1, nb + 1)
        sigma = opti.variable(1, nb + 1) if self.question == 4 else ca.DM.ones(1, nb + 1)

        # Fixed state immediately after first-stage separation.
        baseline = Simulator(self.rk).simulate(
            t_coast=0.0, controller=None,
            t_shut=self.rk.t_burn1 + 1e-6, rtol=1e-10,
        )
        x_sep = self.scale_state(baseline.stages[-1].x[0])
        opti.subject_to(xc[:, 0] == x_sep)
        opti.subject_to(xb[:, 0] == xc[:, -1])

        opti.subject_to(opti.bounded(0.0, tc, 400.0))
        burn_upper = self.rk.t_burn2 if self.question == 3 else self.rk.t_burn2 / SIGMA_MIN
        opti.subject_to(opti.bounded(1.0, tb, burn_upper))
        opti.subject_to(opti.bounded(-0.5 * np.pi, phi, 0.5 * np.pi))
        opti.subject_to(xc[4, :] == 1.0)
        opti.subject_to(opti.bounded(
            (self.rk.ms2 + self.rk.m_payload) / self.m_scale,
            xb[4, :], 1.0,
        ))
        if self.question == 4:
            opti.subject_to(opti.bounded(SIGMA_MIN, sigma, 1.0))

        hc = tc / nc
        for k_idx in range(nc):
            x0, x1 = xc[:, k_idx], xc[:, k_idx + 1]
            f0, f1 = self._coast_dynamics(x0), self._coast_dynamics(x1)
            self._hs_constraint(opti, x0, x1, f0, f1,
                                self._coast_dynamics, hc)

        hb = tb / nb
        for k_idx in range(nb):
            x0, x1 = xb[:, k_idx], xb[:, k_idx + 1]
            p0, p1 = phi[0, k_idx], phi[0, k_idx + 1]
            s0, s1 = sigma[0, k_idx], sigma[0, k_idx + 1]
            f0 = self._dynamics(x0, p0, s0)
            f1 = self._dynamics(x1, p1, s1)
            pmid, smid = 0.5 * (p0 + p1), 0.5 * (s0 + s1)
            self._hs_constraint(
                opti, x0, x1, f0, f1,
                lambda xmid, p=pmid, s=smid: self._dynamics(xmid, p, s), hb,
            )

        rf = xb[0:2, -1]
        vf = xb[2:4, -1]
        opti.subject_to(ca.sumsqr(rf) == 1.0)
        opti.subject_to(ca.dot(rf, vf) == 0.0)
        opti.subject_to(rf[0] * vf[1] - rf[1] * vf[0] == 1.0)

        propellant_used = (1.0 - xb[4, -1]) * self.m_scale
        opti.minimize(propellant_used)

        tc0, tb0, phi0, k0, reference = self._reference_trajectory()
        if warm is None:
            coast_t = np.linspace(self.rk.t_burn1, self.rk.t_burn1 + tc0, nc + 1)
            burn_t = np.linspace(self.rk.t_burn1 + tc0,
                                 self.rk.t_burn1 + tc0 + tb0, nb + 1)
            opti.set_initial(tc, tc0)
            opti.set_initial(tb, tb0)
            opti.set_initial(xc, self.scale_state(self._sample_state(reference, coast_t)))
            opti.set_initial(xb, self.scale_state(self._sample_state(reference, burn_t)))
            opti.set_initial(phi, phi0 + k0 * (burn_t - burn_t[0]))
            if self.question == 4:
                opti.set_initial(sigma, 1.0)
        else:
            opti.set_initial(tc, warm.coast_duration)
            opti.set_initial(tb, warm.burn_duration)
            old_c = np.linspace(0.0, 1.0, warm.coast_states.shape[1])
            old_b = np.linspace(0.0, 1.0, warm.burn_states.shape[1])
            new_c = np.linspace(0.0, 1.0, nc + 1)
            new_b = np.linspace(0.0, 1.0, nb + 1)
            xc_guess = np.vstack([np.interp(new_c, old_c, row) for row in warm.coast_states])
            xb_guess = np.vstack([np.interp(new_b, old_b, row) for row in warm.burn_states])
            opti.set_initial(xc, xc_guess)
            opti.set_initial(xb, xb_guess)
            opti.set_initial(phi, np.interp(new_b, old_b, warm.phi_deg) * DEG)
            if self.question == 4:
                opti.set_initial(sigma, np.interp(new_b, old_b, warm.sigma))

        options = {
            "expand": True,
            "ipopt.print_level": print_level,
            "ipopt.max_iter": 1200,
            "ipopt.tol": 1e-8,
            "ipopt.acceptable_tol": 1e-6,
            "print_time": bool(print_level),
        }
        opti.solver("ipopt", options)
        solved = opti.solve()

        xc_value = np.asarray(solved.value(xc))
        xb_value = np.asarray(solved.value(xb))
        phi_value = np.asarray(solved.value(phi)).reshape(-1) / DEG
        sigma_value = (np.asarray(solved.value(sigma)).reshape(-1)
                       if self.question == 4 else np.ones(nb + 1))
        stats = solved.stats()
        return CollocationSolution(
            question=self.question,
            n_coast=nc,
            n_burn=nb,
            coast_duration=float(solved.value(tc)),
            burn_duration=float(solved.value(tb)),
            coast_states=xc_value,
            burn_states=xb_value,
            phi_deg=phi_value,
            sigma=sigma_value,
            propellant_used=float(solved.value(propellant_used)),
            objective=float(solved.value(propellant_used)),
            max_defect=float(stats.get("iterations", {}).get("inf_pr", [np.nan])[-1]),
        )

    def validate(self, solution: CollocationSolution, rtol: float = 1e-11):
        t2 = self.rk.t_burn1 + solution.coast_duration
        t3 = t2 + solution.burn_duration
        tau_nodes = np.linspace(0.0, 1.0, solution.n_burn + 1)

        def controller(t):
            tau = np.clip((t - t2) / solution.burn_duration, 0.0, 1.0)
            return np.interp(tau, tau_nodes, solution.phi_deg) * DEG

        def sigma_fn(t):
            tau = np.clip((t - t2) / solution.burn_duration, 0.0, 1.0)
            return float(np.interp(tau, tau_nodes, solution.sigma))

        result = Simulator(self.rk).simulate(
            t_coast=solution.coast_duration,
            controller=controller,
            sigma=sigma_fn,
            t_shut=t3,
            rtol=rtol,
        )
        final = result.x[-1]
        radius = np.hypot(final[0], final[1])
        vr = (final[0] * final[2] + final[1] * final[3]) / radius
        vt = (-final[1] * final[2] + final[0] * final[3]) / radius
        return result, {
            "height_error_m": radius - self.rk.a_target,
            "radial_velocity_error_mps": vr,
            "tangential_velocity_error_mps": vt - self.rk.v_circular,
            "residual_norm": float(np.linalg.norm(orbit_residual(final, self.rk))),
            "propellant_used_kg": self.rk.m_after_sep - final[4],
        }


def save_summary(solution: CollocationSolution, validation: dict):
    row = {
        "question": solution.question,
        "n_coast": solution.n_coast,
        "n_burn": solution.n_burn,
        "t_coast_s": solution.coast_duration,
        "burn_duration_s": solution.burn_duration,
        "t_shut_s": RocketParams().t_burn1 + solution.coast_duration + solution.burn_duration,
        "propellant_used_collocation_kg": solution.propellant_used,
        "propellant_used_reintegration_kg": validation["propellant_used_kg"],
        "height_error_m": validation["height_error_m"],
        "radial_velocity_error_mps": validation["radial_velocity_error_mps"],
        "tangential_velocity_error_mps": validation["tangential_velocity_error_mps"],
        "residual_norm": validation["residual_norm"],
        "phi_nodes_deg": " | ".join(f"{v:.6f}" for v in solution.phi_deg),
        "sigma_nodes": " | ".join(f"{v:.6f}" for v in solution.sigma),
    }
    path = RES_DIR / f"q{solution.question}_collocation_summary.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--question", type=int, choices=(3, 4), required=True)
    parser.add_argument("--n-coast", type=int, default=10)
    parser.add_argument("--n-burn", type=int, default=20)
    parser.add_argument("--print-level", type=int, default=0)
    args = parser.parse_args()

    solver = DirectCollocationSolver(args.question, args.n_coast, args.n_burn)
    solution = solver.solve(print_level=args.print_level)
    _trajectory, validation = solver.validate(solution)
    path = save_summary(solution, validation)
    print(f"Q{args.question} Hermite-Simpson: coast={solution.coast_duration:.6f} s, "
          f"burn={solution.burn_duration:.6f} s, fuel={solution.propellant_used:.3f} kg")
    print("DOP853 reintegration:", validation)
    print("saved:", path)


if __name__ == "__main__":
    main()
