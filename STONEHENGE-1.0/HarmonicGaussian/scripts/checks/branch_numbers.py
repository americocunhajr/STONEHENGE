"""Numbers of Sects. 5-7 taken from the stored branch data: peak-gap law ratios, peak power against
its linear value, spread among barrier shapes and with offset, peak and fold with hardening, and the
operating frequency of the calibrated quartic models relative to the Gaussian."""

import numpy as np
from harmonic_gaussian.model import beff, ctot, Wlin
from harmonic_gaussian.theory import M0_GAUSS as M0
from harmonic_gaussian.storage import load

S = load("scaling")
SH = load("shapes")
H = load("hardening")
chi, kap, lam = 0.05, 0.5, 0.05
be = beff(lam)
print("== peak-gap law ratios (continuation gap / law)")
for xi in [0.005, 0.01, 0.02, 0.04]:
    c = ctot(xi, lam)
    Wl = Wlin(xi, lam, bounds=(0.9, 1.1))
    r = np.array(
        [
            (o["f"] / c, (Wl**2 - o["W_peak"] ** 2) / ((2 * M0 / np.pi) * (c * np.sqrt(be) / o["f"]) ** 3))
            for o in S["gauss"][xi]
            if o.get("ok")
        ]
    )
    m = r[:, 0] >= 10
    print(
        f" xi={xi} Wlin={Wl:.6f}  f/c>=10: {r[m, 1].min():.3f}-{r[m, 1].max():.3f} (f/c up to {r[:, 0].max():.0f});  all: {r[:, 1].min():.3f}-{r[:, 1].max():.3f} (f/c from {r[:, 0].min():.1f})"
    )
print("== peak power vs linear value at f=0.6355, xi=0.01")
o = [o for o in S["gauss"][0.01] if o["f"] == 0.6355][0]
c = ctot(0.01, lam)
print(
    " Pe",
    o["Pe_peak"],
    " chi*lam*kap*f^2/(2c^2*beff)",
    chi * lam * kap * 0.6355**2 / (2 * c**2 * be),
    " A_p",
    o["A_peak"],
    " f/(c sqrt beff)",
    0.6355 / (c * np.sqrt(be)),
)
print("== shapes: peak spread and Pe spread; folds")
for i, f in enumerate(SH["fs"]):
    W = [SH[k][i]["W_peak"] for k in ["gauss", "lorentz", "dipole"]]
    P = [SH[k][i]["Pe_peak"] for k in ["gauss", "lorentz", "dipole"]]
    Wf = [SH[k][i].get("W_fold", np.nan) for k in ["gauss", "dipole", "lorentz"]]
    print(
        f"  f={f:<7} peak spread {max(W)-min(W):.1e}  Pe rel spread {(max(P)-min(P))/np.mean(P):.1e}  folds G/D/L {Wf[0]:.4f} {Wf[1]:.4f} {Wf[2]:.4f}"
    )
print("== offset delta=0.1 vs symmetric, f>=0.0868")
dW = []
dF = []
dP = []
for g, b in zip(SH["gauss"], SH["gauss_d0.1"]):
    if g["f"] < 0.08:
        continue
    dW.append(abs(g["W_peak"] - b["W_peak"]))
    dP.append(abs(g["Pe_peak"] - b["Pe_peak"]) / g["Pe_peak"])
    if np.isfinite(g.get("W_fold", np.nan)):
        dF.append(abs(g["W_fold"] - b["W_fold"]))
print(f"  max |dW_peak| {max(dW):.1e}  max |dW_fold| {max(dF):.1e}  max rel dPe {max(dP):.1e}")
print("== hardening k3=1e-4: peak and fold at f>=0.4")
for o in H["res"]:
    if abs(o["k3h"] - 1e-4) < 1e-9 and o["f"] > 0.4:
        print(f"  f={o['f']} peak {o['W_peak']:.4f} fold {o.get('W_fold', np.nan):.4f}")
print("== calibrated Duffing / Gaussian operating-frequency ratio, f>=0.09")
G = {o["f"]: o["W_peak"] for o in S["gauss"][0.01]}
rat = [
    S["calib"][n][i]["W_fold"] / G[f]
    for n in S["calib"]
    for i, f in enumerate(S["fs"])
    if f >= 0.09 and S["calib"][n][i].get("ok")
]
print(f"  min {min(rat):.2f} max {max(rat):.2f}")
spread = []
for i, f in enumerate(S["fs"]):
    v = [S["calib"][n][i]["W_fold"] for n in S["calib"] if S["calib"][n][i].get("ok")]
    if len(v) == 3:
        spread.append((max(v) - min(v)) / min(v))
print("  spread among calibrations (max-min)/min: %.2f-%.2f" % (min(spread), max(spread)))
