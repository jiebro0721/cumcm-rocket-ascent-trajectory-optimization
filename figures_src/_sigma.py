"""节流下限灵敏度（全部 5 档并行）。"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2\src")
from common import RocketParams

rk = RocketParams()


def run_sigma(smin):
    import q4_throttle_opt as q4
    old = q4.SIGMA_MIN
    q4.SIGMA_MIN = smin
    opt = q4.ThrottleOptimizer(rk)
    starts = q4.make_starts_q4(rk, opt.n_phi, opt.n_sigma)
    starts = starts[:40]
    results = opt.solve(starts, maxiter=60, workers=8)
    q4.SIGMA_MIN = old
    if not results:
        return None
    z, prop, nres = results[0]
    sigma_nodes = z[1 + (opt.n_phi - 1): 1 + (opt.n_phi - 1) + opt.n_sigma]
    return smin, float(prop), float(nres), list(np.round(sigma_nodes, 3))


if __name__ == "__main__":
    for s in [0.6, 0.7, 0.8, 0.9, 1.0]:
        r = run_sigma(s)
        if r:
            print("sigma_min=%.1f  最优燃料=%.1f kg  残差=%.2e  sigma*=%s" % (r[0], r[1], r[2], r[3]))
        else:
            print("sigma_min=%.1f  无收敛解" % s)
