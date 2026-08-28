"""两级运载火箭发射轨迹建模：公共工具模块（问题 1~4 共用）。

提供：物理常数、火箭参数、指数大气模型、平面点质量动力学（惯性 ECI 分量，
含地球自转与大气阻力）、分段数值积分器、入轨条件残差、绘图与结果输出。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# 物理常数（SI）
# ---------------------------------------------------------------------------
MU = 3.986004418e14        # 地球引力常数 [m^3/s^2]
R_E = 6371.0e3             # 地球半径（平均）[m]
OMEGA_E = 7.2921159e-5     # 地球自转角速度 [rad/s]
G0 = 9.80665               # 海平面重力加速度 [m/s^2]
RHO0 = 1.225               # 海平面大气密度 [kg/m^3]
H_SCALE = 7200.0           # 指数大气标高 [m]

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "figures"
RES_DIR = ROOT / "results"
FIG_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

DEG = np.pi / 180.0


# ---------------------------------------------------------------------------
# 火箭参数（题目给定）
# ---------------------------------------------------------------------------
@dataclass
class RocketParams:
    """两级液体运载火箭参数（题面给定值）。"""

    # 一子级
    M0: float = 500_000.0        # 起飞质量（含有效载荷与二子级）[kg]
    ms1: float = 40_000.0        # 一子级结构质量 [kg]
    mp1: float = 380_000.0       # 一子级推进剂质量 [kg]
    Isp1: float = 300.0          # 真空比冲 [s]
    Fmax1: float = 7.5e6         # 额定最大推力 [N]
    S: float = 12.5              # 气动阻力参考面积 [m^2]
    CD: float = 0.3              # 阻力系数（忽略马赫数变化）
    # 二子级
    ms2: float = 8_000.0         # 二子级结构质量 [kg]
    mp2: float = 62_000.0        # 二子级推进剂质量 [kg]
    m_payload: float = 10_000.0  # 有效载荷质量 [kg]
    Isp2: float = 420.0          # 真空比冲 [s]
    Fmax2: float = 6.0e5         # 额定最大推力 [N]
    # 任务
    H_target: float = 400.0e3    # 目标圆轨道高度 [m]
    t_vertical: float = 10.0     # 垂直起飞时长 [s]
    pitch_rate_q1: float = 0.4 * DEG  # 问题 1 程序转弯俯仰角角速率 [rad/s]
    t_coast_q1: float = 60.0     # 问题 1 滑行时间 [s]

    @property
    def mdot1(self) -> float:
        """一子级质量流量 = F/(Isp·g0)（真空、全推力）[kg/s]。"""
        return self.Fmax1 / (self.Isp1 * G0)

    @property
    def mdot2(self) -> float:
        """二子级额定质量流量 [kg/s]。"""
        return self.Fmax2 / (self.Isp2 * G0)

    @property
    def t_burn1(self) -> float:
        """一子级推进剂耗尽时间 [s]（一级全推力恒额定，故为定值）。"""
        return self.mp1 / self.mdot1

    @property
    def t_burn2(self) -> float:
        """二子级推进剂耗尽时间 [s]。"""
        return self.mp2 / self.mdot2

    @property
    def m_after_sep(self) -> float:
        """一子级分离瞬间（关机后）剩余质量 [kg]：M0 - mp1 - ms1。"""
        return self.M0 - self.mp1 - self.ms1

    @property
    def a_target(self) -> float:
        """目标圆轨道地心距 [m]。"""
        return R_E + self.H_target

    @property
    def v_circular(self) -> float:
        """目标圆轨道惯性速度 sqrt(mu/a) [m/s]。"""
        return np.sqrt(MU / self.a_target)


# ---------------------------------------------------------------------------
# 大气模型
# ---------------------------------------------------------------------------
def rho_air(h: np.ndarray | float) -> np.ndarray | float:
    """指数大气密度模型 rho = rho0 * exp(-h / h0)。"""
    return RHO0 * np.exp(-h / H_SCALE)


# ---------------------------------------------------------------------------
# 平面点质量动力学（惯性 ECI 笛卡尔分量，赤道平面内）
# ---------------------------------------------------------------------------
def rel_velocity(rx: np.ndarray, ry: np.ndarray, vx: np.ndarray, vy: np.ndarray):
    """相对大气速度 v_rel = v_in - omega_E x r（赤道平面：omega_E x r = om*[-y, x]）。"""
    return vx + OMEGA_E * ry, vy - OMEGA_E * rx


def flight_path_angle(vrx: np.ndarray, vry: np.ndarray, rx: np.ndarray, ry: np.ndarray):
    """相对速度飞行路径角 gamma = atan2(v·r_hat, v·t_hat) [rad]。

    r_hat = r/|r|，t_hat = z_hat x r_hat（东向）。
    垂直上升段 v~=0 时输出无定义，入轨段（高空）有效。
    """
    r = np.hypot(rx, ry)
    r_dot = (vrx * rx + vry * ry) / r          # 径向分量
    t_dot = (vrx * (-ry) + vry * rx) / r       # 切向（东向）分量
    return np.arctan2(r_dot, t_dot)


def gamma_inertial(vx: np.ndarray, vy: np.ndarray, rx: np.ndarray, ry: np.ndarray):
    """惯性速度飞行路径角（与入轨条件 r·v=0 直接对应）[rad]。"""
    r = np.hypot(rx, ry)
    r_dot = (vx * rx + vy * ry) / r
    t_dot = (vx * (-ry) + vy * rx) / r
    return np.arctan2(r_dot, t_dot)


def drag_accel(
    vrx: np.ndarray, vry: np.ndarray,
    rx: np.ndarray, ry: np.ndarray,
    S: float, CD: float,
) -> tuple[np.ndarray, np.ndarray]:
    """气动阻力加速度 [m/s^2]：a_D = -0.5*rho*S*CD*|v_rel|*v_rel / m。

    返回 (ax, ay)（惯性分量）。此处不含质量，由调用者再除以 m。
    """
    r = np.hypot(rx, ry)
    h = r - R_E
    rho = rho_air(h)
    vrel = np.hypot(vrx, vry)
    q = 0.5 * rho * S * CD * vrel
    return -q * vrx, -q * vry


def thrust_unit(phi: np.ndarray, rx: np.ndarray, ry: np.ndarray):
    """由俯仰角 phi（推力与当地水平面的夹角）构造推力单位向量 û（惯性分量）。

    û = cos(phi) * t_hat + sin(phi) * r_hat。
    """
    r = np.hypot(rx, ry)
    rx_h, ry_h = rx / r, ry / r              # 径向单位向量
    tx_h, ty_h = -ry / r, rx / r             # 东向单位向量
    return tx_h * np.cos(phi) + rx_h * np.sin(phi), ty_h * np.cos(phi) + ry_h * np.sin(phi)


# ---------------------------------------------------------------------------
# 分段数值积分器
# ---------------------------------------------------------------------------
@dataclass
class Segment:
    """一段飞行轨迹的积分结果。"""

    t: np.ndarray
    x: np.ndarray   # [x, y, vx, vy, m]
    vrx: np.ndarray
    vry: np.ndarray
    gamma: np.ndarray     # 惯性速度飞行路径角
    gamma_rel: np.ndarray  # 相对速度飞行路径角
    phi: np.ndarray  # 俯仰角
    T: np.ndarray    # 推力
    sigma: np.ndarray
    name: str = ""

    @property
    def h(self) -> np.ndarray:
        return np.hypot(self.x[:, 0], self.x[:, 1]) - R_E

    @property
    def v_in(self) -> np.ndarray:
        return np.hypot(self.x[:, 2], self.x[:, 3])

    @property
    def m(self) -> np.ndarray:
        return self.x[:, 4]

    @property
    def time(self) -> np.ndarray:
        return self.t


@dataclass
class SimResult:
    """完整仿真结果（各阶段拼接）。"""

    t: np.ndarray
    x: np.ndarray             # 分段拼接后的状态
    vrx: np.ndarray
    vry: np.ndarray
    gamma: np.ndarray         # 惯性速度飞行路径角
    gamma_rel: np.ndarray     # 相对速度飞行路径角
    phi: np.ndarray
    T: np.ndarray
    sigma: np.ndarray
    phase: list[str] = field(default_factory=list)  # 每点的阶段名
    t1: float | None = None     # 一子级燃尽时刻
    t2: float | None = None     # 二子级点火时刻
    t_shut: float | None = None  # 二子级关机（实际）时刻
    stages: list[Segment] = field(default_factory=list)

    @property
    def h(self) -> np.ndarray:
        return np.hypot(self.x[:, 0], self.x[:, 1]) - R_E

    @property
    def v_in(self) -> np.ndarray:
        return np.hypot(self.x[:, 2], self.x[:, 3])

    @property
    def m(self) -> np.ndarray:
        return self.x[:, 4]

    @property
    def gamma_deg(self) -> np.ndarray:
        return self.gamma / DEG

    def final_state(self) -> dict:
        return {
            "h_km": float(np.hypot(self.x[-1, 0], self.x[-1, 1]) / 1e3 - R_E / 1e3),
            "v_in_mps": float(np.hypot(self.x[-1, 2], self.x[-1, 3])),
            "gamma_deg": float(self.gamma[-1] / DEG),
            "gamma_rel_deg": float(self.gamma_rel[-1] / DEG),
            "m_kg": float(self.x[-1, 4]),
        }


class Simulator:
    """飞行段积分器：垂直段 -> 程序转弯段 -> 分离 -> 滑行段 -> 二级动力段。"""

    def __init__(self, rk: RocketParams | None = None):
        self.rk = rk or RocketParams()

    # -- 各段右端函数 ------------------------------------------------------
    def _rhs(
        self,
        t: float, x: np.ndarray,
        T_const: float, sigma: float,
        phi: float, along_vel: bool,
        drag_S: float, drag_CD: float,
        Isp: float,
    ) -> np.ndarray:
        rx, ry, vx, vy, m = x
        T = T_const * sigma
        if along_vel:
            vrx, vry = rel_velocity(rx, ry, vx, vy)
            vrel = np.hypot(vrx, vry)
            tx_h, ty_h = (vrx / vrel, vry / vrel) if vrel > 1e-9 else (0.0, 1.0)
        else:
            tx_h, ty_h = thrust_unit(phi, rx, ry)
        r = np.hypot(rx, ry)
        ax = -MU / r**3 * rx + T / m * tx_h
        ay = -MU / r**3 * ry + T / m * ty_h
        if drag_S > 0.0:
            vrx, vry = rel_velocity(rx, ry, vx, vy)
            dax, day = drag_accel(vrx, vry, rx, ry, drag_S, drag_CD)
            ax += dax / m
            ay += day / m
        return np.array([vx, vy, ax, ay, -T / (Isp * G0)])

    # -- 单段积分 ----------------------------------------------------------
    def _integrate(
        self,
        t0: float, tf: float, y0: np.ndarray,
        T_const: float, sigma: float,
        phi_fn, along_vel: bool,
        drag_S: float, drag_CD: float,
        rtol: float, atol: np.ndarray,
        name: str,
    ) -> Segment:
        def rhs(t, x):
            phi = phi_fn(t) if phi_fn is not None else 0.0
            return self._rhs(
                t, x, T_const, sigma, phi, along_vel, drag_S, drag_CD, self.rk.Isp1
            )

        sol = solve_ivp(
            rhs, (t0, tf), y0, method="DOP853",
            rtol=rtol, atol=atol, dense_output=True,
        )
        t = np.linspace(t0, tf, max(2, int((tf - t0) * 20) + 1))
        x = sol.sol(t).T
        vrx, vry = rel_velocity(x[:, 0], x[:, 1], x[:, 2], x[:, 3])
        gamma = gamma_inertial(x[:, 2], x[:, 3], x[:, 0], x[:, 1])
        gamma_rel = flight_path_angle(vrx, vry, x[:, 0], x[:, 1])
        if along_vel:
            phi = np.where(
                np.hypot(vrx, vry) > 1e-9,
                np.arctan2(vrx * x[:, 0] + vry * x[:, 1],
                           vrx * (-x[:, 1]) + vry * x[:, 0]),
                np.nan,
            )
        else:
            phi = np.array([phi_fn(ti) if phi_fn is not None else np.nan for ti in t])
        T = np.full_like(t, T_const * sigma)
        return Segment(t, x, vrx, vry, gamma, gamma_rel, phi, T, np.full_like(t, sigma), name)

    # -- 全流程 ------------------------------------------------------------
    def simulate(
        self,
        t_coast: float = 60.0,
        controller=None,
        t_shut: float | None = None,
        sigma: float | callable = 1.0,
        rtol: float = 1e-9,
        drag_after_sep: bool = False,
    ) -> SimResult:
        """按六段剖面仿真。

        controller(t) -> phi [rad]：二级动力段俯仰角程序；None 表示攻角为零
        （真空二级段推力方向与惯性速度方向一致）。
        sigma: 二级节流比（常数或函数 t->sigma，默认 1.0）。
        t_shut: 二级关机时刻（绝对时间）。None 表示按满推力燃尽时间关机；
        若给定更晚时刻，则由干质量事件限制累计推进剂消耗。
        返回 SimResult（全阶段拼接）。
        """
        rk = self.rk
        st = []
        atol = np.array([1e-2, 1e-2, 1e-4, 1e-4, 1e-2])

        # --- 1) 垂直起飞段 [0, t_vert]：推力沿当地径向，全推力，有阻力 ---
        y0 = np.array([R_E, 0.0, 0.0, OMEGA_E * R_E, rk.M0])
        phi = 0.5 * np.pi
        seg = self._integrate(
            0.0, rk.t_vertical, y0,
            rk.Fmax1, 1.0, lambda t: phi, False,
            rk.S, rk.CD, rtol, atol, "vertical",
        )
        st.append(seg)

        # --- 2) 程序转弯段 [t_vert, t_burn1]：phi 以 0.4deg/s 线性减小 ---
        phi_fn = lambda t: 0.5 * np.pi - rk.pitch_rate_q1 * (t - rk.t_vertical)
        seg = self._integrate(
            rk.t_vertical, rk.t_burn1, seg.x[-1],
            rk.Fmax1, 1.0, phi_fn, False,
            rk.S, rk.CD, rtol, atol, "pitch-over",
        )
        st.append(seg)

        # --- 3) 一子级分离（质量跳变）---
        x1 = seg.x[-1].copy()
        x1[4] -= rk.ms1

        # --- 4) 无动力滑行段 [t1, t1+t_coast]：无阻力 ---
        t1 = rk.t_burn1
        t2 = t1 + t_coast
        seg = self._integrate(
            t1, t2, x1,
            0.0, 1.0, lambda t: 0.0, False,
            0.0, 0.0, rtol, atol, "coast",
        )
        st.append(seg)

        # --- 5) 二级动力段 [t2, t_shut] ---
        mdot2 = rk.mdot2
        sigma_fn = sigma if callable(sigma) else (lambda t: float(sigma))
        if t_shut is None:
            t_shut = t2 + rk.t_burn2   # 燃尽

        def mass_flow(t):
            return -sigma_fn(t) * rk.Fmax2 / (rk.Isp2 * G0)

        def rhs2(t, x):
            sig = sigma_fn(t)
            T = sig * rk.Fmax2
            if controller is None:
                # 推力方向与惯性速度方向一致（二子级段无大气阻力，无相对大气参照）
                vx, vy = x[2], x[3]
                vmag = np.hypot(vx, vy)
                tx_h, ty_h = (vx / vmag, vy / vmag) if vmag > 1e-9 else (0.0, 1.0)
            else:
                p = controller(t)
                tx_h, ty_h = thrust_unit(p, x[0], x[1])
            r = np.hypot(x[0], x[1])
            ax = -MU / r**3 * x[0] + T / x[4] * tx_h
            ay = -MU / r**3 * x[1] + T / x[4] * ty_h
            return np.array([x[2], x[3], ax, ay, mass_flow(t)])

        sol = solve_ivp(
            rhs2, (t2, t_shut), seg.x[-1], method="DOP853",
            rtol=rtol, atol=atol,
            events=[lambda t, x: x[4] - (rk.ms2 + rk.m_payload)],
            dense_output=True,
        )
        if sol.t_events[0].size > 0:
            t_shut = float(sol.t_events[0][0])
        tt = np.linspace(t2, t_shut, max(2, int((t_shut - t2) * 20) + 1))
        xs = sol.sol(tt).T
        sig_arr = np.array([sigma_fn(ti) for ti in tt])
        vrx, vry = rel_velocity(xs[:, 0], xs[:, 1], xs[:, 2], xs[:, 3])
        gamma = gamma_inertial(xs[:, 2], xs[:, 3], xs[:, 0], xs[:, 1])
        gamma_rel = flight_path_angle(vrx, vry, xs[:, 0], xs[:, 1])
        if controller is None:
            phi_arr = np.where(
                np.hypot(xs[:, 2], xs[:, 3]) > 1e-9,
                np.arctan2(xs[:, 2] * xs[:, 0] + xs[:, 3] * xs[:, 1],
                           xs[:, 2] * (-xs[:, 1]) + xs[:, 3] * xs[:, 0]),
                np.nan,
            )
        else:
            phi_arr = np.array([controller(ti) for ti in tt])
        T_arr = sig_arr * rk.Fmax2
        seg2 = Segment(tt, xs, vrx, vry, gamma, gamma_rel, phi_arr, T_arr, sig_arr, "stage2")
        st.append(seg2)

        # --- 拼接 ---
        t = np.concatenate([s.t for s in st])
        x = np.vstack([s.x for s in st])
        vrx = np.concatenate([s.vrx for s in st])
        vry = np.concatenate([s.vry for s in st])
        gamma = np.concatenate([s.gamma for s in st])
        gamma_rel = np.concatenate([s.gamma_rel for s in st])
        phi = np.concatenate([s.phi for s in st])
        T = np.concatenate([s.T for s in st])
        sig = np.concatenate([s.sigma for s in st])
        names = []
        for s in st:
            names += [s.name] * len(s.t)
        return SimResult(t, x, vrx, vry, gamma, gamma_rel, phi, T, sig, names,
                         t1, t2, t_shut, st)


# ---------------------------------------------------------------------------
# 入轨条件残差（问题 2~4 共用）
# ---------------------------------------------------------------------------
def orbit_residual(xf: np.ndarray, rk: RocketParams) -> np.ndarray:
    """由关机时刻惯性状态计算入轨条件残差（无量纲，顺行圆轨道）。

    三个分量依次为：半径偏差、径向速度偏差、切向速度偏差。
    切向速度取顺行方向（局部东向），当 v_t<vc 时残差为负，
    从而严格排除逆行圆轨道分支。
    xf = [x, y, vx, vy, m]（惯性）。
    """
    r = np.hypot(xf[0], xf[1])
    rx, ry = xf[0] / r, xf[1] / r          # 径向单位向量
    tx, ty = -ry, rx                       # 东向（顺行）单位向量
    vr = xf[2] * rx + xf[3] * ry
    vt = xf[2] * tx + xf[3] * ty
    a = rk.a_target
    vc = rk.v_circular
    return np.array([
        (r - a) / R_E,
        vr / vc,
        (vt - vc) / vc,
    ])


def fmt_residual(res: np.ndarray) -> str:
    return "[" + ", ".join(f"{v:+.2e}" for v in res) + "]"
