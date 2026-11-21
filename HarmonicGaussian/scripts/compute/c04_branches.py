"""Whole interwell and intrawell branches of the Gaussian barrier and of the quartic model matched
to (a, k_w), at f = 0.083, 0.15, 0.3 and 0.6 (Fig. 5).

Rows of each branch: (Omega, A, Pe, max|mu|, max x, min x).
Output: data/branches.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params, A_WELL
from harmonic_gaussian.orbits import march, newton_shoot, continuation
from harmonic_gaussian.branches import solve_case, solve_duffing, hb_state, pack, DUFF
from harmonic_gaussian.storage import main, pick, log

FS = [0.083, 0.15, 0.3, 0.6]


def interwell(p, z, Wmax, Astop):
    up = continuation(
        z, p, direction=+1, Wmax=Wmax, Wmin=0.3, ds=0.05, ds_max=0.25, maxpts=8000 if Wmax > 3 else 6000, Astop=Astop
    )
    dn = continuation(z, p, direction=-1, Wmax=Wmax, Wmin=0.3, ds=0.05, ds_max=0.25, maxpts=6000)
    return np.vstack([pack(dn)[::-1], pack(up)])


def intrawell(p, Wstart):
    """Intrawell branch in the right well, continued down in frequency from Wstart."""
    pp = p.copy()
    pp[5] = Wstart
    z = march(np.array([A_WELL, 0.0, 0.0]), pp, 512, 300)
    rec = continuation(z, pp, direction=-1, Wmax=Wstart + 0.01, Wmin=0.3, ds=0.02, ds_max=0.15, maxpts=8000)
    return pack(rec)[::-1]


def compute():
    out = {}
    k1, k3 = DUFF["(a,k_w)"]
    for f in pick(FS, [0.15]):
        o = solve_case(f=f, down=False)  # Gaussian: restart from harmonic balance at the peak
        p = params("gauss", f=f)
        p[5] = o["W_peak"]
        z, ok, _ = newton_shoot(hb_state(p, o["W_peak"])[0], p, 512)
        out[("gauss", f)] = dict(high=interwell(p, z, 2.2, 1.0), intra=intrawell(params("gauss", f=f), 1.8))
        log("gauss", f, "done")
        od = solve_duffing(k1, k3, f)  # quartic: harmonic balance below its jump
        pd = params("duffing", f=f, s=(k1, k3, 0, 0))
        z = None
        for fr in [0.6, 0.7, 0.5, 0.8, 0.4]:
            W0 = max(fr * od["W_fold"], 0.9)
            pd[5] = W0
            r = hb_state(pd, W0, Amin=0.3, Amax=300)
            if r is None:
                continue
            zz, ok, _ = newton_shoot(r[0], pd, 512)
            if ok and np.max(np.abs(zz)) > 1.0:
                z = zz
                break
        out[("duff", f)] = dict(
            high=interwell(pd, z, 12, 0.95), intra=intrawell(params("duffing", f=f, s=(k1, k3, 0, 0)), 9.0)
        )
        log("duffing", f, "done")
    return out


def report(d):
    for key, v in sorted(d.items()):
        hb, ib = v["high"], v["intra"]
        st = hb[:, 3] < 1 + 1e-6
        print(
            f"{key[0]:6s} f={key[1]:<6} interwell: {len(hb)} points, Omega {hb[:, 0].min():.3f}-{hb[:, 0].max():.3f}, max Pe {hb[st, 2].max():.4e};"
            f"  intrawell: Omega {ib[:, 0].min():.2f}-{ib[:, 0].max():.2f}, {int(np.sum(ib[:, 5] > 0))} of {len(ib)} points with min x > 0"
        )


if __name__ == "__main__":
    main("branches", compute, report)
