"""Geometric and inertial nonlinearity of the first mode of the inextensional cantilever (Crespo da
Silva and Glynn 1978), used to size the beam hardening k3 in Sects. 6 and 9.4."""

import numpy as np
from scipy.integrate import cumulative_trapezoid as ct, trapezoid

b = 1.87510407
sg = (np.cosh(b) + np.cos(b)) / (np.sinh(b) + np.sin(b))
s = np.linspace(0, 1, 200001)
phi = np.cosh(b * s) - np.cos(b * s) - sg * (np.sinh(b * s) - np.sin(b * s))
nrm = np.sqrt(trapezoid(phi**2, s))
c = 1 / nrm
P1 = c * b * (np.sinh(b * s) + np.sin(b * s) - sg * (np.cosh(b * s) - np.cos(b * s)))  # phi'
P2 = c * b**2 * (np.cosh(b * s) + np.cos(b * s) - sg * (np.sinh(b * s) + np.sin(b * s)))  # phi''
a1 = 2 * trapezoid(P1**2 * P2**2, s)  # geometric
I1 = ct(ct(P1**2, s, initial=0), s, initial=0)
I1 = I1 - I1[-1]
a2 = -trapezoid(P1**2 * I1, s)  # inertial
cb = 3 * a1 / (8 * b**4) - a2 / 4
print("analytic: alpha1=", a1, " alpha2=", a2, " backbone coeff a^2: %.4f, (w/l)^2: %.4f" % (cb, cb / 4))
