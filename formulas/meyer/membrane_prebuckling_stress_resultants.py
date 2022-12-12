"""
Contains membrane pre-buckling stress resultants for Meyer's shell buckling analysis. These were derived by Flugge, and
assume membrane condition of pre-buckling stress in the cylinder. This neglects any pre-buckling bending due to edge
constraint at the edge of the cylinder due to end closures in actual usage or due to the loading fixture under test
conditions.

Inclusion of these pre-buckling bending effects for the case of uniform compression results in somewhat lower critical
loads than those calculated here.

Source(s): Buckling of Stiffened Cylindrical Shells Subjected to Combined Axial Compression, Normal Pressure, Bending (R.R. Meyer)
"""

import numpy as np


# Let the circular cylinder be loaded by the following external forces:
# F = applied axial load
# M = applied end moment
# V = applied shearing load
# p = uniform normal surface pressure

# Note that:
# N_a : axial stress resultant maximum (phi = 0)
# N_b : bending stress resultant maximum (phi = 0)
# V_o : shear stress resultant maximum (phi = 0)
# N_x_o: membrane pre-buckling stress resultant, in the x-direction (axially)
# N_phi_o: membrane pre-buckling stress resultant, in the phi-direction (circumferentially)


def _N_a(F: float, R: float) -> float:
    """
    Axial pressure
    """
    return F / (np.pi * R ** 2)


def _N_b(M: float, R: float) -> float:
    """
    Bending pressure
    """
    return M / (np.pi * R ** 2)


def _V_o(V: float, R: float) -> float:
    """
    Shear pressure
    """
    return V / (np.pi * R)


