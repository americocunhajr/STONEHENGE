# Harvester optimization

Optimizes mean electrical output power of a bistable harvester using cross-entropy optimization or direct/grid search, including noisy-excitation cases.

## Files and entry points

- HarvesterOpt_CrossEntropy2var.m / HarvesterOpt_CrossEntropy4var.m — cross-entropy design search
- HarvesterOpt_DirectSearch2var.m / HarvesterOpt_DirectSearch4var.m — direct/grid design search
- *Noise.m — noisy-excitation variants
- PiezoMagBeam_ObjFunc*.m / PiezoMagBeam_RHS.m — objective and dynamics
- CrossEntropyOpt_Example.m and graph_*.m — example and visualization

## Running and dependencies

Each HarvesterOpt_*.m MATLAB script specifies a separate design problem and serves as its entry point.

The existing numerical code was not modified or independently execution-tested as part of this documentation update.

[Back to the STONEHENGE-1.0 module guide](../README.md).
