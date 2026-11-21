# Stochastic simulations

Simulates harvester trajectories and electrical output in the presence of colored-noise disturbances.

## Files and entry points

- HarvesterColoredNoise_ivp_main.m — principal experiment
- HarvesterColoredNoise.m — ODE right-hand side with time-varying noise
- HarvesterColoredNoise_ivp_post.m and post_piezomagbeam_attractors_cnoise.m — post-processing
- fredholm_expcorr_eig.m — exponential-covariance integral-operator eigenpairs
- HarvesterPower.m, randvar_*.m, graph_*.m, plot_*.m — analysis utilities

## Running and dependencies

Inspect noise-generation and numerical sampling settings in the MATLAB main script before running.

The existing numerical code was not modified or independently execution-tested as part of this documentation update.

[Back to the STONEHENGE-1.0 module guide](../README.md).
