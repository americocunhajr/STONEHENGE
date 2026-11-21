"""Numbers of Sect. 8.1 from the map of final states (Fig. 8a-d)."""

import numpy as np
from harmonic_gaussian.analysis import map_fields
from harmonic_gaussian.storage import load

X = load("final_states")
M = map_fields(X)
Ws, fs, P1, PS, PC, PI, GAIN, TS = (M[k] for k in ["Ws", "fs", "P1", "PS", "PC", "PI", "GAIN", "TS"])
print(f'grid {len(Ws)} x {len(fs)} (rows done), {X["nstarts"]} starts per cell, horizon {X["horizon"]} periods')
r = np.where(((P1 + PS + PC) > 0).any(1))[0]
print("first row with any interwell start: f=%.4f (previous %.4f)" % (fs[r[0]], fs[r[0] - 1]))
i = np.where((P1 > 0).any(1))[0][0]
j = np.where(P1[i] > 0)[0]
print(
    "first row with capture: f=%.4f (previous %.4f), Omega %.3f-%.3f, max P1 %.3f"
    % (fs[i], fs[i - 1], Ws[j[0]], Ws[j[-1]], P1[i].max())
)
jW = int(np.argmin(abs(Ws - 1.0)))


def band(i, q):
    """Contiguous run of cells with P1 >= q containing Omega = 0.97, or the extreme cells when there is none."""
    ok = P1[i] >= q
    j0 = int(np.argmin(abs(Ws - 0.97)))
    if not ok[j0]:
        j = np.where(ok)[0]
        return ("no run at 0.97", round(Ws[j[0]], 3), round(Ws[j[-1]], 3)) if len(j) else None
    lo = hi = j0
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(Ws) - 1 and ok[hi + 1]:
        hi += 1
    return (round(Ws[lo], 3), round(Ws[hi], 3))


for f0 in [0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 0.8, 1.0]:
    i = int(np.argmin(abs(fs - f0)))
    print(
        "f=%.3f  P1>=0.5 %s  P1>=0.95 %s  at Omega=1: P1 %.2f sub %.2f chaos %.2f intra %.2f  max gain %.2f at %.3f"
        % (
            fs[i],
            band(i, 0.5),
            band(i, 0.95),
            P1[i, jW],
            PS[i, jW],
            PC[i, jW],
            PI[i, jW],
            GAIN[i].max(),
            Ws[np.argmax(GAIN[i])],
        )
    )
i = np.where(P1[:, jW] >= 0.5)[0][0]
print("P1(Omega=1)>=0.5 first at f=%.4f (previous %.4f)" % (fs[i], fs[i - 1]))
m = fs >= 0.1
print(
    "cells with f>=0.1: max share period 2,3 interwell %.2f, cells where it exceeds 0.25: %d of %d"
    % (PS[m].max(), (PS[m] > 0.25).sum(), PS[m].size)
)
i, j = np.unravel_index(np.argmax(PS), PS.shape)
print("  largest subharmonic share at f=%.3f Omega=%.3f" % (fs[i], Ws[j]))
dom = (PC > 0.5).any(1)
print(
    "chaotic interwell dominant (>0.5) in %d cells; f range %.3f-%.3f"
    % ((PC > 0.5).sum(), fs[dom].min(), fs[dom].max())
)
i, j = np.unravel_index(GAIN.argmax(), GAIN.shape)
print("max gain %.1f at f=%.3f Omega=%.3f" % (GAIN.max(), fs[i], Ws[j]))
for f0 in [0.05, 0.15, 0.3, 0.6]:
    i = int(np.argmin(abs(fs - f0)))
    j = np.where(GAIN[i] > 1)[0]
    runs = [
        (round(Ws[a], 3), round(Ws[b], 3))
        for a, b in zip(j[np.r_[True, np.diff(j) > 1]], j[np.r_[np.diff(j) > 1, True]])
    ]
    print("f=%.3f gain>1 for Omega in" % fs[i], runs)
for f0 in [0.22, 0.23, 0.25, 0.3, 0.5]:
    i = int(np.argmin(abs(fs - f0)))
    print(
        "settling median at Omega=1, f=%.4f: %.0f periods; band-median over captured cells %.0f"
        % (fs[i], TS[i, jW], np.nanmedian(TS[i]))
    )
print("settling: 5th, 50th and 95th percentiles of the cell medians:", np.round(np.nanpercentile(TS, [5, 50, 95]), 0))
