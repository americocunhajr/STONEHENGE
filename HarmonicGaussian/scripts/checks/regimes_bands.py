"""Operating regimes and expected power (Sect. 9, Fig. 10): capture thresholds from the capture
curves, cusp from the fold curve, and the expected power from rest integrated over three frequency
bands relative to the linear harvester."""

from harmonic_gaussian.analysis import capture_thresholds, band_ratio
from harmonic_gaussian.storage import load

P = load("expected_power")
Ws, fs, PE = P["Ws"], P["fs"], P["PE"]
for W0, W1 in [(0.4, 1.4), (0.9, 1.1), (0.4, 0.9)]:
    print((W0, W1), [f"{f}:{band_ratio(Ws, PE[i], f, W0, W1):.3f}" for i, f in enumerate(fs)])
print("cusp f_c", load("fold_curve")["fc"])
print(
    "capture thresholds at 50% (c: f below - f at or above):",
    "  ".join(f"{c:.5f}: {a:.4f}-{b:.4f}" for c, a, b in capture_thresholds(load("capture_curves"))),
)
