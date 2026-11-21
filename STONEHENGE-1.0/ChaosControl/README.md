# Chaos control

Contains two approaches to controlling nonlinear and chaotic responses of the harvester: OGY control and extended time-delayed feedback (ETDF).

## Files and entry points

- main_OGY_R01.m — OGY control experiment
- main_ETDF.m — ETDF control experiment
- Harvest_DuffingCL.slx — Simulink model invoked by ETDF
- piezomagbeam.m, plot_time_series.m, z1test.m — dynamics, plotting, and 0–1 chaos diagnostic

## Running and dependencies

The ETDF example requires Simulink in addition to MATLAB; the OGY and ETDF files are separate starting points.

The existing numerical code was not modified or independently execution-tested as part of this documentation update.

[Back to the STONEHENGE-1.0 module guide](../README.md).
