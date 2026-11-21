"""Regression tests: the library reproduces entries of the stored data used in the paper.

Run from the root folder with ``python -m pytest tests`` or ``python -m tests.test_reproduction``
(about one minute). Trajectory-based tests use cases far from fractal basin boundaries, where
outcomes do not depend on round-off.
"""

import numpy as np
from harmonic_gaussian.model import params, dU, A_WELL, ctot, wells
from harmonic_gaussian.orbits import flow
from harmonic_gaussian.branches import solve_case
from harmonic_gaussian.ensembles import run_capture, run_map, rest_starts, interwell
from harmonic_gaussian.theory import omega_inter, M0_GAUSS
from harmonic_gaussian.storage import load


def test_force_and_wells():
    p = params("gauss")
    assert abs(dU(A_WELL, p)) < 1e-14
    assert np.allclose(wells(p), (-A_WELL, A_WELL), atol=1e-12)


def test_area_law():
    U = lambda x: x * x / 2 + np.exp(-x * x)
    A = 12.0
    deficit = 1 - omega_inter(U, A) ** 2
    assert abs(deficit / (2 * M0_GAUSS / (np.pi * A**3)) - 1) < 0.02


def test_peak_and_fold_of_stored_branch():
    o = [o for o in load("scaling")["gauss"][0.01] if abs(o["f"] - 0.3366) < 1e-9][0]
    q = solve_case(f=o["f"], xi=0.01)
    assert abs(q["W_peak"] - o["W_peak"]) < 1e-9 and abs(q["W_fold"] - o["W_fold"]) < 1e-9
    assert abs(q["W_low"] - o["W_low"]) < 1e-9 and q["low_type"] == o["low_type"]


def test_monodromy_of_damped_flow_contracts_volume():
    p = params("gauss", f=0.0, W=1.0)
    _, M, *_ = flow(np.array([0.5, 0.0, 0.0]), p, 512)
    # Liouville: det M = exp(-(2 xi + lam) T) for the linear damping of (x, x', v)
    assert abs(np.linalg.det(M) - np.exp(-(2 * 0.01 + 0.05) * 2 * np.pi)) < 1e-9


def test_capture_times_at_f030():
    T = [r for r in load("crisis") if r["f"] == 0.3][0]["T"]
    x0, y0, phi = rest_starts(64, -A_WELL)
    o = run_capture(x0[:8], y0[:8], phi[:8], params("gauss", f=0.3), 1.0, 2000, 256, 0.5 * 0.3 / ctot(0.01, 0.05))
    assert np.array_equal(o[:, 0], T[:8])


def test_final_state_cell():
    X = load("final_states")
    i = int(np.argmin(abs(X["fs"] - 0.3)))
    j = int(np.argmin(abs(X["Ws"] - 1.0)))
    x0, y0, phi = rest_starts(32, -A_WELL)
    n = x0.size
    o = run_map(x0, y0, phi, np.full(n, X["Ws"][j]), np.full(n, X["fs"][i]), params("gauss"), 1000, 128, 1e-5)
    p1 = np.mean((o[:, 0] == 1) & (o[:, 1] == 1) & interwell(o[:, 3], o[:, 4]))
    assert p1 == X["K"][i, j, 0]
    assert abs(np.median(o[:, 5]) - X["K"][i, j, 6]) < 1e-12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(name, "ok")
