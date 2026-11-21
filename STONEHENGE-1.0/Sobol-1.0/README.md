# Sobol sensitivity analysis

Performs variance-based global sensitivity analysis of harvester outputs, using Monte Carlo and polynomial-chaos-based Sobol indices.

## Files and entry points

- main_pmeh_sobol.m — uncertain parameters and analysis driver
- sobol_mc.m — Monte Carlo Sobol analysis
- sobol_pce.m — polynomial-chaos Sobol analysis
- piezo_magneto_harvesting.m — response model
- plot_sobol.m — plots

## Running and dependencies

Requires MATLAB with UQLab (calls uqlab, uq_createModel, uq_createAnalysis). The driver initializes a pseudorandom seed.

The existing numerical code was not modified or independently execution-tested as part of this documentation update.

[Back to the STONEHENGE-1.0 module guide](../README.md).
