"""Capture from rest onto the resonant orbit near its own peak frequency (three damping levels and
four hardening levels), and slow frequency sweeps from Omega = 0.4 up to the peak frequency
(Fig. 8f), for the symmetric barrier and the barrier offset by delta = 0.1.

Inputs: data/scaling.pkl, data/hardening.pkl. Output: data/capture_at_peak.pkl.
"""

import numpy as np
from harmonic_gaussian.model import params, ctot, beff, Wlin, wells
from harmonic_gaussian.theory import peak_law, M0_GAUSS
from harmonic_gaussian.ensembles import run_fixed, run_ramp, rest_starts, interwell
from harmonic_gaussian.storage import main, pick, load, log

NPH = pick(64, 8)


def cell(p, W, Apk):
    """Fractions of starts on the resonant orbit, on any periodic interwell orbit, and irregular; mean power."""
    x0, y0, phi = rest_starts(NPH, *wells(params("gauss")))
    o = run_fixed(x0, y0, phi, p, W, pick(1000, 200), 40, 256 if Apk < 15 else 512)
    inter = interwell(o[:, 0], o[:, 1])
    per = o[:, 3] < 1e-3
    big = np.maximum(o[:, 0], -o[:, 1]) >= 0.8 * Apk
    return dict(
        cap_res=float(np.mean(inter & per & big)),
        cap_any=float(np.mean(inter & per)),
        irregular=float(np.mean(inter & ~per)),
        Pe=float(np.mean(o[:, 2])),
    )


def compute():
    S = load("scaling")
    Hd = load("hardening")
    out = dict(cases={})
    for xi in [0.005, 0.01, 0.02]:
        c = ctot(xi, 0.05)
        rows = []
        for o in S["gauss"][xi]:
            if not o.get("ok") or not (2 <= o["f"] / c <= 30):
                continue
            p = params("gauss", xi=xi, lam=0.05, f=o["f"])
            cells = [cell(p, o["W_peak"] + dW, o["A_peak"]) for dW in (-0.02, -0.01, 0.0)]
            rows.append(dict(f=o["f"], fc=o["f"] / c, W_peak=o["W_peak"], A_peak=o["A_peak"], cells=cells))
            log(f'xi={xi} f/c={o["f"]/c:5.1f}', [round(q["cap_res"], 3) for q in cells])
        out["cases"][("xi", xi)] = rows
    for k3 in Hd["k3s"][3:7]:
        rows = []
        for o in Hd["res"]:
            if abs(o["k3h"] - k3) > 1e-12 or o["f"] > 0.53 or not o.get("ok"):
                continue
            p = params("gauss", xi=0.01, lam=0.05, f=o["f"], k3h=k3)
            cells = [cell(p, o["W_peak"] + dW, o["A_peak"]) for dW in (-0.02, -0.01, 0.0)]
            rows.append(
                dict(f=o["f"], fc=o["f"] / ctot(0.01, 0.05), W_peak=o["W_peak"], A_peak=o["A_peak"], cells=cells)
            )
            log(f'k3={k3:.2e} f={o["f"]:.4f}', [round(q["cap_res"], 3) for q in cells])
        out["cases"][("k3", k3)] = rows
    # sweeps from Omega = 0.4 to the peak frequency of the area law, 128 starts, time step 0.02
    Wl = Wlin(0.01, 0.05)
    be = beff(0.05)
    c = ctot(0.01, 0.05)
    ramps = []
    for delta in [0.0, 0.1]:
        x0, y0, phi = rest_starts(NPH, *wells(params("gauss", delta=delta)))
        for f in pick([0.12, 0.15, 0.2], [0.15]):
            Wp = peak_law(f, Wl, c, be, M0_GAUSS)
            pp = params("gauss", f=f, delta=delta)
            for Tr in pick([250, 500, 1000, 2000, 4000, 8000, 16000], [250, 1000]):
                o = run_ramp(x0, phi, pp, 0.4, Wp, float(Tr), 2000.0, 0.02, 251.0)
                ramps.append(
                    dict(delta=delta, f=f, W1=Wp, Tr=Tr, cap=float(np.mean(np.maximum(o[:, 0], -o[:, 1]) > 3)))
                )
                log(f'sweep delta={delta} f={f} Tr={Tr}: {ramps[-1]["cap"]:.3f}')
    out["ramps"] = ramps
    return out


def report(R):
    print("Capture from rest near the peak frequency (Omega_peak - 0.02, - 0.01, + 0)")
    for k, rows in R["cases"].items():
        print(k)
        for r in rows:
            cs = r["cells"]
            print(
                f"  f={r['f']:.4f} f/c={r['fc']:5.1f} W_peak={r['W_peak']:.4f}  capture(max over 3 freqs)={max(q['cap_any'] for q in cs):.2f}  irregular at peak={cs[2]['irregular']:.2f}"
            )
    print("Slow sweeps from Omega = 0.4 to the peak frequency: fraction on the resonant branch")
    for r in R["ramps"]:
        print(f"  delta={r['delta']} f={r['f']} W1={r['W1']:.4f} Tr={r['Tr']}: {r['cap']:.3f}")


if __name__ == "__main__":
    main("capture_at_peak", compute, report)
