"""节点收敛性：N_phi = 5/7/9（充分初值并行）。"""
import sys
import numpy as np

sys.path.insert(0, r"C:\Users\Hanamiya\Desktop\2026国赛\培训\机理2\src")
from common import RocketParams

rk = RocketParams()


def run_nphi(n_phi, n_start):
    from q3_fuel_opt import FuelOptimizer, make_starts
    opt = FuelOptimizer(rk, n_phi=n_phi)
    starts = make_starts(rk, n_phi)[:n_start]
    results = opt.solve(starts, maxiter=60, workers=8)
    if not results:
        return None
    z, prop, nres = results[0]
    return n_phi, float(prop), float(nres)


if __name__ == "__main__":
    for n in [5, 7, 9]:
        r = run_nphi(n, 40)
        if r:
            print("N_phi=%d  最优燃料=%.1f kg  残差=%.2e" % (r[0], r[1], r[2]))
        else:
            print("N_phi=%d  无收敛解" % n)
