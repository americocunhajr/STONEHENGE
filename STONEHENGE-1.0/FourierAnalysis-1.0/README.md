# Fourier analysis

Characterizes the frequency content of computed harvester trajectories using one-sided Fourier transforms and power spectral density estimates.

## Files and entry points

- main_piezomagbeam_fourier.m — simulation and spectrum driver
- post_piezomagbeam_fourier.m — post-processing
- signal_fft_1sided.m — one-sided FFT
- signal_psd.m, signal_psd_1sided.m — PSD calculations
- graph_type2_psd.m — spectral plotting
- piezomagbeam.m — dynamical model

## Running and dependencies

The sampling settings in the MATLAB driver affect the frequency resolution and spectral estimates.

The existing numerical code was not modified or independently execution-tested as part of this documentation update.

[Back to the STONEHENGE-1.0 module guide](../README.md).
