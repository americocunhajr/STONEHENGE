---
layout: default
title: "STONEHENGE Modules — Bistable Energy Harvesting and Nonlinear Dynamics"
description: "Discover 13 open-source MATLAB, C++ and Python modules for piezoelectric vibration energy harvesting, nonlinear dynamics, bifurcations, chaos control, stochastic excitation and sensitivity analysis."
permalink: /modules/
---

# STONEHENGE scientific software modules

**STONEHENGE — Suite for Nonlinear Analysis of Energy Harvesting Systems** is an open-source research software collection for modeling, simulating, optimizing and analyzing nonlinear vibration energy harvesters. It covers bistable piezoelectric–magneto–elastic systems and includes MATLAB, C++ and Python implementations.

Each module has its own source-code directory and README. The 13 modules below are organized by their primary scientific use; they are not a single executable application.

[Visit the STONEHENGE homepage]({{ '/' | relative_url }}) · [Browse the complete source repository](https://github.com/americocunhajr/STONEHENGE) · [Read the module guide](https://github.com/americocunhajr/STONEHENGE/blob/master/STONEHENGE-1.0/README.md)

## Nonlinear dynamics and response analysis

### Time-series analysis

[**TimeSeriesAnalysis**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/TimeSeriesAnalysis) (MATLAB) integrates the coupled mechanical and piezoelectric initial-value problem and studies displacement, velocity and voltage histories. The model permits potential asymmetry and nonlinear electromechanical coupling.

### Fourier analysis

[**FourierAnalysis**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/FourierAnalysis) (MATLAB) analyzes vibration and electrical-response spectra using Fourier transforms and power spectral densities. It includes one-sided FFT and PSD routines.

### Poincaré sections

[**PoincareSections**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/PoincareSections) (MATLAB) constructs phase portraits and stroboscopic Poincaré sections of forced bistable harvesters, helping distinguish periodic responses and more complicated dynamics.

### Resonance and frequency response

[**ResonanceCurves**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/ResonanceCurves) (MATLAB) computes nonlinear resonance curves under harmonic excitation by varying the driving frequency.

### Bifurcation analysis

[**BifurcationAnalysis**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/BifurcationAnalysis) (MATLAB) generates bifurcation diagrams from forcing-amplitude sweeps at different frequencies, including forward/backward sweeps and three-dimensional visualization.

### Basins of attraction

[**BasinsOfAttraction**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/BasinsOfAttraction) (C++ and MATLAB) simulates multiple initial conditions to classify long-time responses and map basins of attraction. The C++ routines perform the numerical calculations and MATLAB scripts provide visualization.

## Visualization, control and design

### Dynamic animation

[**DynamicAnimation**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/DynamicAnimation) (MATLAB) visualizes the motion of piezo-magneto-elastic beams in inertial and moving frames, including symmetric/asymmetric models and GIF output.

### Chaos control

[**ChaosControl**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/ChaosControl) (MATLAB and Simulink) studies stabilization of nonlinear and chaotic harvester responses using the Ott–Grebogi–Yorke (OGY) method and extended time-delayed feedback (ETDF).

### Energy-harvester optimization

[**HarvesterOptimization**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/HarvesterOptimization) (MATLAB) seeks designs that maximize mean electrical power using cross-entropy optimization and direct/grid search, with deterministic and noise-perturbed excitation variants.

### Displacement amplification

[**DisplacementAmplifier**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/DisplacementAmplifier) (MATLAB) studies a bistable energy harvester coupled to a displacement-amplification mechanism and driven by base excitation.

## Uncertainty quantification and stochastic systems

### Global sensitivity analysis

[**SobolSensitivity**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/SobolSensitivity) (MATLAB with UQLab) estimates Sobol variance-based sensitivity indices using Monte Carlo sampling and polynomial-chaos-expansion metamodels.

### Stochastic excitation

[**StochasticSimulations**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/StochasticSimulations) (MATLAB) simulates bistable energy harvesters under colored-noise disturbances and analyzes their trajectories, attractors, and electrical output.

## Harmonic–Gaussian energy harvesting research

[**HarmonicGaussian**](https://github.com/americocunhajr/STONEHENGE/tree/master/STONEHENGE-1.0/HarmonicGaussian) (Python) reproduces numerical results for *Harmonic-Gaussian vibration energy harvesting*. It studies localized Gaussian and point-dipole magnetic barriers, nonlinear resonance, periodic-orbit continuation, Floquet stability, capture, and electrical power, with comparisons to Duffing-type potentials.

The [HarmonicGaussian module README](https://github.com/americocunhajr/STONEHENGE/blob/master/STONEHENGE-1.0/HarmonicGaussian/README.md) documents its Python dependencies, reproducibility scripts, stored data, figures and regression tests.

## Questions about STONEHENGE

### Is STONEHENGE free and open source?

Yes. STONEHENGE is distributed under the [MIT license](https://github.com/americocunhajr/STONEHENGE/blob/master/LICENSE). Some module workflows require separately installed research software such as MATLAB, Simulink or UQLab.

### Which programming languages are supported?

The legacy numerical modules use MATLAB and C++; the HarmonicGaussian module uses Python. Dependencies and execution commands vary by module.

### What physical systems does STONEHENGE model?

Its principal application is vibration energy harvesting using bistable, electromechanically coupled piezoelectric–magneto–elastic oscillators. The modules cover forced nonlinear response, parameter sensitivity, stochastic excitation, chaotic dynamics and performance optimization.

### How should the software be cited?

Cite the associated [STONEHENGE paper in *Software Impacts* (2021)](https://doi.org/10.1016/j.simpa.2021.100161), DOI: **10.1016/j.simpa.2021.100161**. Individual research modules may have additional publication references in their own documentation.
