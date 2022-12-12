"""
Contains formulas related to membrane stresses

Assumption: External loads don't cause any additional deformation, such that we can consider these loads alone for the
dome-cylinder connection region.

Source(s): https://engineeringlibrary.org/reference/simple-thin-pressure-vessels-air-force-stress-manual [Section 8.3.1]
--> https://apps.dtic.mil/sti/pdfs/AD0759199.pdf

 In a thin pressure vessel, no stresses other than those tangential to the surface are present at points sufficiently
 removed from discontinuities in the curvature, slope, or thickness of the wall. These tangential or membrane
 stresses are constant throughout the thickness of the shell. At points near discontinuities, such as the junction of
 a cylinder and its head discontinuity, stresses must be superposed upon the membrane stresses in order to obtain the
 total stress.
"""
from typing import Tuple, Any

import numpy as np


# Arrange the membrane stress formulas into a system of equations to be solved. N_mmer is the stress in the direction of
# the meridian line times the shell thickness and N_mt is the stress in the direction of a circle or rotation times the
# shell thickness. dp is the pressure difference between the inside and outside the vessel.
# A * N_mmer + B * N_mt = dp
# C * N_mmer + D * N_mt = 0
# A = 1 / R_mer
# B = 1 / R_t
# C = 1
# D = (-dp * R_t) / 2

def _membrane_stresses_thin_shell_of_revolution(R_mer: float, R_t: float, dp: float) -> tuple[Any, Any]:
    """
    Equation 8-1 from the Air Force Stress Manual.

    :param R_mer: Meridional radius of curvature
    :param R_t: Tengential radius of curvature
    :param dp: Pressure differential (inside and outside the pressure vessel
    :return: Equation 8-1 from the Air Force Stress Manual
    """
    A = 1 / R_mer
    B = 1 / R_t
    C = 1
    D = (-dp * R_t) / 2

    # Solve system of equations with least squares
    LHS = np.array([[A, B], [C, D]])
    RHS = np.array([dp, 0])
    N_mmer, N_mt = np.linalg.lstsq(LHS, RHS, rcond=None)[0]

    return N_mmer, N_mt

