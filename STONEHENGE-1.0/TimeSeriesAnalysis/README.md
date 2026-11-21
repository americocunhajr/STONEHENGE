# Time-series analysis

Integrates the nonlinear electromechanical initial-value problem and plots displacement, velocity, and voltage, allowing for asymmetry and nonlinear coupling.

## Files and entry points

- main_pmehna_ivp.m — simulation parameters, integration and time-history plots
- pmehna_eom.m — governing differential equations and ODE integration

## Running and dependencies

The MATLAB driver is the entry point; the model contains bias/asymmetry and nonlinear piezoelectric coupling terms.

The existing numerical code was not modified or independently execution-tested as part of this documentation update.

[Back to the STONEHENGE-1.0 module guide](../README.md).
