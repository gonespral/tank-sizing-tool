"""
Contains formulas related to pressure.

Source(s): my mind, my brain, my heart, my soul, my body, my blood, my sweat, my tears, my life, my death, my everything
"""

import numpy as np
from databases.fos import *

def find_hydrostatic_pressure(rho: float, h: float, a: float) -> float:
    """
    Finds the hydrostatic pressure at a given height.
    :param rho: Density of fluid in kg/m^3
    :param h: Height of fluid in m
    :param a: Acceleration
    :return:
    """
    if h < 0:
        return 0
    else:
        return rho * h * a * Coeff_A_Loads * Coeff_B


def find_stress_in_dome_from_internal_pressure(p: float, r: float, t: float):
    """
    Finds the stress in a dome from an internal pressure.
    :param p: internal pressure in Pa
    :param r: radius of dome in m
    :param t: thickness of dome in m
    :return: stress in Pa
    """
    return p * Coeff_A_Pressure * Coeff_B * r / (2 * t)
