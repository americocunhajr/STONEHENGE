# Harmonic-Gaussian vibration energy harvesting

Code and data that reproduce every figure and every computed number of

> A. Cunha Jr, P. Wolszczak, G. Litak, *Harmonic-Gaussian vibration energy harvesting*.

This folder is part of [STONEHENGE](https://stonehengesuite.org), the Suite for Nonlinear Analysis of
Energy Harvesting Systems ([americocunhajr/STONEHENGE](https://github.com/americocunhajr/STONEHENGE)).
It models a piezoelectric cantilever made bistable by magnets. The magnetic force is a localized barrier
(Gaussian or point-dipole) on the linear beam, and it is compared with the quartic (Duffing) model.

## Quick start

```bash
pip install -r requirements.txt          # numpy, scipy, numba, matplotlib, mpmath
python run_all.py figures                # Figures 1-10 into figures/, from the stored data (about 30 s)
python run_all.py results                # every printed number into results/ (about 1 min)
python -m pytest tests                   # regression tests against the stored data (a few seconds)
```

Run these commands from this folder; a single script runs as, e.g., `python -m scripts.figures.fig08_access`.
Python 3.10 or later is needed. The figures use matplotlib's
bundled Computer Modern font. Figure 1 is a TikZ drawing; building it needs `pdflatex` with `circuitikz`,
and it is skipped when `pdflatex` is missing.

## Layout

```
harmonic_gaussian/     the library (imported by every script)
scripts/compute/       c01 ... c23: compute and store a dataset, then print its summary
scripts/figures/       fig01 ... fig10: draw one figure of the paper from the stored data
scripts/checks/        analytic checks and numbers of the paper computed from the stored data
data/                  stored datasets used in the paper (pickled dicts of numpy arrays, 11 MB)
figures/               figures of the paper (PDF)
results/               printed output of every compute and check script
tests/                 regression tests
run_all.py             driver
```

## The library

| Module | Content |
|---|---|
| `model` | Parameter vector and baseline values; restoring force and stiffness of the Gaussian, Lorentzian, dipole-type, point-dipole magnet pair and quartic barriers, with optional hardening and offset (compiled with Numba); well positions; linear beam and circuit: total damping, `beta_eff`, `Omega_lin`, power of the linear harvester. |
| `orbits` | Period map with variational equations (RK4 in the rescaled time), Newton shooting, pseudo-arclength continuation in Omega, Floquet multipliers, refinement of peak and fold, two-parameter continuation of the fold curve and its cusp. |
| `branches` | Harmonic-balance initial guesses and the solver that returns the peak of power, the fold and the low end of the stable interwell band (`solve_case`, `solve_duffing`); barrier shapes of equal area; quartic calibrations. |
| `ensembles` | Parallel ensembles of starts from rest: fixed frequency (`run_fixed`), slow sweeps (`run_ramp`), first-capture times (`run_capture`), Poincaré sections (`strobe`), final states with early stopping (`run_map`). |
| `theory` | Free-oscillation frequencies by quadrature (double precision and mpmath), single-harmonic stiffness, harmonics of the barrier force, the Q invariant, the similarity parameter and the peak law. |
| `analysis` | Post-processing shared by a figure and the script that prints its numbers: crisis lifetimes and fit, capture thresholds, band-averaged power, fields of the map of final states. |
| `storage`, `style` | Paths, loading and saving, quick mode; figure style. |

The parameter vector `p = params(shape, xi, chi, kap, lam, f, W, s, k3h, delta)` carries everything the
compiled kernels need. The baseline is `xi = 0.01`, `chi = 0.05`, `kap = 0.5`, `lam = 0.05` and the
Gaussian barrier `exp(-x^2)` on `x^2/2`.

## Figures

| Figure | Script | Data (compute script) |
|---|---|---|
| 1 Device and model | `fig01_device.tex` | none |
| 2 Potential and phase plane | `fig02_concept.py` | none |
| 3 Q invariant; identification of the area | `fig03_statics.py` | `identification` (c13) |
| 4 Backbones and the area law | `fig04_backbones.py` | none (quadrature) |
| 5 Interwell branches | `fig05_response.py` | `branches` (c04) |
| 6 Similarity scaling and bandwidth | `fig06_scaling.py` | `similarity` (c06, from c05), `magnet` (c14) |
| 7 Quartic models | `fig07_quartic.py` | `scaling` (c01), `quartic_fits` (c12), `harmonics` (c11) |
| 8 Reaching the resonance from rest | `fig08_access.py` | `final_states` (c23), `capture_curves` (c22), `capture_at_peak` (c09), `scaling` (c01), `lowend` (c07), `melnikov` (c10) |
| 9 Boundary crisis | `fig09_mechanism.py` | `basins` (c21), `crisis_fine` (c17) |
| 10 Regimes and expected power | `fig10_regimes.py` | `fold_curve` (c18), `capture_curves` (c22), `expected_power` (c19) |

## Datasets

| Script | Dataset | Content |
|---|---|---|
| c01 | `scaling` | Peak and fold against forcing: Gaussian at four damping levels, quartic models |
| c02 | `hardening` | Peak and fold with beam hardening |
| c03 | `shapes` | Peak and fold for three barrier shapes of equal area, and with offset |
| c04 | `branches` | Whole interwell and intrawell branches at four forcing levels |
| c05 | `branches_eps` | Branches for the similarity analysis (damping, load, barrier size, quartics) |
| c06 | `similarity` | Peak gap, fold distance, bandwidth against the similarity parameter |
| c07 | `lowend` | Low end of the stable band (pitchfork) |
| c08 | `cusp` | Turning points near the cusp |
| c09 | `capture_at_peak` | Capture from rest near the peak frequency; slow frequency sweeps |
| c10 | `melnikov` | Melnikov threshold against frequency |
| c11 | `harmonics` | Harmonics of computed orbits against the closed form |
| c12 | `quartic_fits` | Quartics fitted to the Gaussian force and their predictions |
| c13 | `identification` | Error of the area identified from a noisy static force |
| c14 | `magnet` | Magnet pair: area, backbone, peak and fold |
| c15 | `magnet_capture` | Magnet pair: capture from rest |
| c16 | `crisis` | Capture times from rest at Omega = 1, up to 20 000 periods |
| c17 | `crisis_fine` | Capture times near the crisis, 256 starts, up to 40 000 periods |
| c18 | `fold_curve` | Fold curve in the (Omega, f) plane and its cusp |
| c19 | `expected_power` | Expected power from rest against frequency |
| c20 | `design` | Design example with the magnet pair and a matched load |
| c21 | `basins` | Basins of the resonant orbit and a Poincaré section |
| c22 | `capture_curves` | Capture against forcing at Omega = 1 for three damping levels |
| c23 | `final_states` | Map of final states, power and settling time (181 x 121 cells) |

Each compute script stores `data/<dataset>.pkl` and prints a summary, which `run_all.py results` regenerates
from the stored data with `--report`. The checks in `scripts/checks/` print the remaining numbers of the
paper: the backbone table (Supplementary Table S1), algebraic tails, harmonics, next-order backbone term,
magnet statics, cantilever nonlinearity, fold laws, cusp fit, step-size convergence, crisis fit, band power
ratios and the statistics of the map of final states.

## Recomputing the data

```bash
python run_all.py compute c22 c23        # selected datasets (scripts run in numerical order)
python run_all.py compute                # all datasets: several hours on two cores
python run_all.py compute --quick        # all scripts on small grids, into data_quick/ (about 2 min)
python run_all.py list                   # all scripts
```

A recomputation overwrites `data/`. Quick mode checks that every script runs; its numbers are not those
of the paper. Every calculation is deterministic: ensembles use uniformly spaced phases, and the one
Monte Carlo study uses a fixed seed. On the same platform a rerun reproduces the stored data, apart from
peak frequencies re-maximized by an optimizer (c01, c05), which agree to about 1e-9. On another platform,
last-digit differences in the math library may change the outcome of individual starts that lie on fractal
basin boundaries (Supplementary Sect. S2.2), but not the statistics.

The compiled functions are cached by Numba. Numba does not notice changes made to `model.py` from
functions compiled in other modules, so after editing it run `python run_all.py clean-cache`.

## Citation

If you use this code, please cite the paper above and the STONEHENGE software paper:

> J. P. Norenberg et al., STONEHENGE -- Suite for nonlinear analysis of energy harvesting systems,
> *Software Impacts* 10 (2021) 100161, doi:10.1016/j.simpa.2021.100161.

## License

MIT, as the STONEHENGE repository.
