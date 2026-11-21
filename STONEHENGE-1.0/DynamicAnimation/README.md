# Dynamic animation

Simulates and animates the motion of symmetric and asymmetric piezo-magneto-elastic cantilevers in inertial and moving reference frames.

## Files and entry points

- main_piezomagbeam_ivp_v2.m — nominal simulation driver
- main_piezomagbeam_asymmetric_ivp.m — asymmetric simulation driver
- plot_piezomagbeam_animation_inertial*.m — inertial-frame views
- plot_piezomagbeam_animation_mobile*.m — moving-frame views, including GIF output
- pmeh_eom.m and piezomagbeam_asymmetric.m — dynamics

## Running and dependencies

Run the appropriate MATLAB main script to prepare the response before plotting or animating.

The existing numerical code was not modified or independently execution-tested as part of this documentation update.

[Back to the STONEHENGE-1.0 module guide](../README.md).
