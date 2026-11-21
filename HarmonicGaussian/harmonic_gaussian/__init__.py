"""Harmonic-Gaussian vibration energy harvesting.

Library behind the paper by A. Cunha Jr, P. Wolszczak and G. Litak. Modules:

    model       parameters, restoring forces of the barrier models, linear beam and circuit
    orbits      period map with variational equations, shooting, continuation, Floquet, fold curve
    branches    harmonic balance and the solver for peak, fold and low end of the interwell branch
    ensembles   parallel ensembles from rest: fixed frequency, sweeps, capture times, sections, maps
    theory      free-oscillation frequencies, harmonics of the barrier force, similarity parameter
    analysis    post-processing shared by a figure and the script printing its numbers
    storage     paths of data, figures and results; quick mode
    style       figure style
"""

__version__ = "1.0.0"
