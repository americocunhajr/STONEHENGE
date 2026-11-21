# STONEHENGE 1.0 — Module guide

**STONEHENGE** (Suite for Nonlinear Analysis of Energy Harvesting Systems) collects numerical tools for studying piezoelectric–magneto–elastic vibration energy harvesters. This directory contains the original MATLAB/C++ modules and the Python `HarmonicGaussian` research package.

This README is an **index of the existing modules**, not a new execution framework. Module names and source files are preserved. Each module is intended to be used from its own directory; its `main_*.m`, `main.cpp`, or `run_all.py` file is generally the entry point.

## Module directory

The following groups organize the modules *conceptually*. No directory has been moved or renamed.

### Response analysis and nonlinear dynamics

| Module | Language | Purpose | Main entry point |
|---|---|---|---|
| [TimeSeriesAnalysis-1.0](TimeSeriesAnalysis-1.0/) | MATLAB | Integrate the nonlinear electromechanical initial-value problem and plot the response, including asymmetry and nonlinear piezoelectric coupling. | `main_pmehna_ivp.m` |
| [FourierAnalysis-1.0](FourierAnalysis-1.0/) | MATLAB | Compute frequency-domain representations of simulated time series, including FFTs and power spectral densities. | `main_piezomagbeam_fourier.m` |
| [PoincareSections-1.0](PoincareSections-1.0/) | MATLAB | Construct phase portraits and stroboscopic Poincaré sections across excitation amplitudes. | `main_PhasePortait_Poincare.m` |
| [ResonanceCurves-1.0](ResonanceCurves-1.0/) | MATLAB | Simulate excitation-frequency sweeps and produce nonlinear resonance-response curves. | `main_piezomagbeam_resocurves.m` |
| [BifurcationAnalysis-1.0](BifurcationAnalysis-1.0/) | MATLAB | Sweep excitation amplitude at several frequencies and plot three-dimensional bifurcation diagrams, including upward and downward sweeps. | `main_3D_bifurcation.m` |
| [BasinsAttraction-1.0](BasinsAttraction-1.0/) | C++ / MATLAB | Integrate grids of initial states, classify asymptotic responses, and visualize basins of attraction. | `main.cpp` |

### Visualization

| Module | Language | Purpose | Main entry point |
|---|---|---|---|
| [DynamicAnimation-1.0](DynamicAnimation-1.0/) | MATLAB | Animate the motion of symmetric and asymmetric piezo-magneto-elastic beams in inertial and moving reference frames; includes GIF output. | `main_piezomagbeam_ivp_v2.m`; `main_piezomagbeam_asymmetric_ivp.m` |

### Control, performance, and design

| Module | Language | Purpose | Main entry point |
|---|---|---|---|
| [ChaosControl-1.0](ChaosControl-1.0/) | MATLAB / Simulink | Explore control of chaotic harvester responses using OGY and extended time-delayed feedback (ETDF). | `main_OGY_R01.m`; `main_ETDF.m` |
| [HarvesterOpt-1.0](HarvesterOpt-1.0/) | MATLAB | Optimize mean harvested power with cross-entropy search or direct/grid search, including 2- and 4-variable and noise-perturbed cases. | `HarvesterOpt_CrossEntropy2var.m`; `HarvesterOpt_DirectSearch4var.m` |
| [DisplacementAmplifier-1.0](DisplacementAmplifier-1.0/) | MATLAB | Simulate a base-excited, two-degree-of-freedom bistable harvester with a mechanical displacement-amplification mechanism. | `main_bistable_amplification_base.m` |

### Uncertainty and stochastic dynamics

| Module | Language | Purpose | Main entry point |
|---|---|---|---|
| [Sobol-1.0](Sobol-1.0/) | MATLAB / UQLab | Estimate global Sobol sensitivity indices by Monte Carlo and polynomial-chaos-expansion approaches. | `main_pmeh_sobol.m` |
| [StochasticSimulations-1.0](StochasticSimulations-1.0/) | MATLAB | Simulate the bistable harvester under colored-noise disturbances and analyze trajectories, attractors, spectra, and harvested power. | `HarvesterColoredNoise_ivp_main.m` |

### Research-specific reproducibility package

| Module | Language | Purpose | Main entry point |
|---|---|---|---|
| [HarmonicGaussian](HarmonicGaussian/) | Python | Reproduce *Harmonic-Gaussian vibration energy harvesting*: localized Gaussian/point-dipole magnetic barriers, continuation, folds, Floquet stability, capture, electrical power, and manuscript figures. | `run_all.py` |

For installation, tests, datasets, individual computations, and figure reproduction, consult the **[HarmonicGaussian README](HarmonicGaussian/README.md)**. Its detailed instructions intentionally remain local to that module.

## What is inside each module?

### [BasinsAttraction-1.0](BasinsAttraction-1.0/)

A C++ implementation of basin mapping for a piezo-magneto-elastic harvester. `main.cpp` sets the initial-displacement/initial-velocity grid, forcing and simulation parameters, and attractor-classification workflow. `header.cpp` and `header.h` implement numerical and auxiliary functions. MATLAB files such as `ploting.m`, `graph_type1.m`, and `graph_contourf_pnt.m` visualize the results. **Note:** this module mixes C++ computation with MATLAB post-processing.

### [BifurcationAnalysis-1.0](BifurcationAnalysis-1.0/)

The driver `main_3D_bifurcation.m` specifies excitation-frequency and forcing-amplitude ranges. `bifurcation_3d.m` computes upward/downward parameter sweeps; `plot_3D_bifurcation.m` visualizes the bifurcation data in three dimensions. `piezomagbeam_asymmetric.m` implements the modeled beam dynamics.

### [ChaosControl-1.0](ChaosControl-1.0/)

Two control approaches are represented: the **Ott–Grebogi–Yorke (OGY)** method in `main_OGY_R01.m`, and **extended time-delayed feedback (ETDF)** in `main_ETDF.m`. The latter invokes the Simulink model `Harvest_DuffingCL.slx`. `piezomagbeam.m`, `plot_time_series.m`, and `z1test.m` support simulation, plotting, and the 0–1 test for chaos. **Dependency:** Simulink is needed for the ETDF example.

### [DisplacementAmplifier-1.0](DisplacementAmplifier-1.0/)

A focused numerical example in `main_bistable_amplification_base.m`. It models a base-excited two-mass system with coupling springs/dampers, nonlinear bistable stiffness, and piezoelectric voltage dynamics. The script integrates the five-state model and plots the resulting response.

### [DynamicAnimation-1.0](DynamicAnimation-1.0/)

Time-domain simulations and animations of the piezo-magneto-elastic beam. The main scripts cover symmetric and asymmetric configurations. The `plot_piezomagbeam_animation_*.m` routines provide inertial-frame and moving-frame visualizations, with an option for GIF output. Model functions include `pmeh_eom.m` and `piezomagbeam_asymmetric.m`.

### [FourierAnalysis-1.0](FourierAnalysis-1.0/)

Fourier-domain characterization of the harvester response. `main_piezomagbeam_fourier.m` generates or processes response histories; `post_piezomagbeam_fourier.m` handles post-processing. `signal_fft_1sided.m`, `signal_psd.m`, and `signal_psd_1sided.m` implement spectral estimators, while `graph_type2_psd.m` creates the associated plots.

### [HarmonicGaussian](HarmonicGaussian/)

An installable Python research package with a reusable `harmonic_gaussian/` library, `scripts/compute/` numerical workflows, `scripts/figures/` figure generators, numerical checks, regression `tests/`, stored `data/`, `figures/`, `results/`, and the `run_all.py` driver. Unlike the older modules, it studies a harmonic beam with a **localized** Gaussian or point-dipole magnetic barrier and compares this model with quartic approximations. See its [module-specific documentation](HarmonicGaussian/README.md).

### [HarvesterOpt-1.0](HarvesterOpt-1.0/)

Optimization tools for maximizing the harvester's mean electrical power. `HarvesterOpt_CrossEntropy*.m` provides cross-entropy-based design search; `HarvesterOpt_DirectSearch*.m` performs deterministic grid/direct search. The filenames distinguish two- and four-design-variable cases and variants with noise. `PiezoMagBeam_ObjFunc*.m` and `PiezoMagBeam_RHS.m` define the objectives and dynamical model. Auxiliary routines plot results and search evolution.

### [PoincareSections-1.0](PoincareSections-1.0/)

Phase-space and stroboscopic analysis for a bistable harvester. `main_PhasePortait_Poincare.m` compares responses over several forcing amplitudes. `phase_portait.m` generates phase-portrait samples, `poincare.m` samples once per excitation period after transients, and `piezomagbeam_asymmetric.m` provides the governing dynamics.

### [ResonanceCurves-1.0](ResonanceCurves-1.0/)

Excitation-frequency response of the classical bistable piezo-magneto-elastic model. `main_piezomagbeam_resocurves.m` is the sweep driver; `piezomagbeam.m` supplies the equations of motion; and `plot_resocurves.m` visualizes the response curves.

### [Sobol-1.0](Sobol-1.0/)

Global sensitivity analysis of the harvesting model using Sobol indices. `main_pmeh_sobol.m` defines uncertain parameters and starts the analysis. `sobol_mc.m` calculates Monte Carlo-based estimates, and `sobol_pce.m` uses a polynomial-chaos metamodel. `piezo_magneto_harvesting.m` defines the model outputs and `plot_sobol.m` plots the indices. **Dependency:** the scripts call UQLab (`uqlab`, `uq_createModel`, `uq_createAnalysis`).

### [StochasticSimulations-1.0](StochasticSimulations-1.0/)

Responses of a bistable harvester subjected to colored-noise excitation. `HarvesterColoredNoise_ivp_main.m` runs the numerical study, `HarvesterColoredNoise.m` defines the stochastic forcing contribution to the ODE, and `HarvesterColoredNoise_ivp_post.m` supports post-processing. `fredholm_expcorr_eig.m` computes eigenpairs for an exponential covariance kernel; the remaining helpers process random variables, power, trajectories, attractors, and plots.

### [TimeSeriesAnalysis-1.0](TimeSeriesAnalysis-1.0/)

A minimal time-domain initial-value analysis. `main_pmehna_ivp.m` sets the parameters and plots displacement, velocity, and voltage histories. `pmehna_eom.m` integrates the electromechanically coupled equations with potential asymmetry, bias, and nonlinear piezoelectric coupling.

## Languages, tools, and reproducibility

- **MATLAB:** most legacy modules consist of scripts and local functions. Start with the indicated main script from the corresponding directory. Some scripts use hard-coded output names or assume a local working directory.
- **C++:** `BasinsAttraction-1.0` contains a standalone source program and supporting C++ files; the results are plotted using MATLAB scripts.
- **Simulink:** needed by the ETDF example in `ChaosControl-1.0`.
- **UQLab:** needed by the sensitivity calculations in `Sobol-1.0`.
- **Python:** `HarmonicGaussian` has its own dependency specification, package metadata, tests, and reproducibility instructions.

The legacy MATLAB/C++ code has been **catalogued, not refactored or execution-tested** as part of this README. No claim is made that all scripts run without external dependencies or environment adjustments. For general STONEHENGE documentation, see the [repository root README](../README.md) and [user manual](../UserManual-1.0.pdf).

## Naming and organization

The `-1.0` suffix is part of the historical names of the original modules. `HarmonicGaussian` is a newer, paper-specific module. Keeping the existing paths avoids breaking references in scripts, publications, and external links.

Any future renaming or physical regrouping should be performed as a separate, reviewed change. The tables above provide logical organization without altering code or directory paths.
