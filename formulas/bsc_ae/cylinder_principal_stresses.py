"""
These functions, used together, find the principal stresses as a function of position.

IMPORTANT NOTE: This analysis is incorrect for pressure vessels of the type being considered for Stratos V. Hence,
this module is not used.

Source(s): mechanics of materials, sad, alex's dreams
"""

import numpy as np


def stress_position(effective_internal_pressure: float, z: float, radius: float, thickness: float,
                    axialforce: float, area: float, shearforce: float, moment: float, lateralforce: float,
                    theta: float, momentofinertia: float) -> tuple:
    """
    Finds stress state in the structure as a function of position.

    :param effective_internal_pressure: effective internal pressure in the cylinder (basic pressure + hydrostatic pressure) in Pa
    :param z: z-coordinate of point of interest in m
    :param radius: cylinder radius in m
    :param thickness: cylinder thickness in m
    :param axialforce: axial point force imposed onto cylinder in N
    :param area: cross-sectional area of cylinder (only material, not space inside) in m^2
    :param shearforce: shear point force imposed onto cylinder in N
    :param moment: moment imposed onto cylinder in Nm
    :param lateralforce: distributed lateral force imposed onto cylinder wall in N/m
    :param theta: theta-coordinate of point of interest in radians
    :param momentofinertia: moment of inertia of the cross-section of the cylinder in m^4
    :return: stress state as a function of position in the cylinder (axial, hoop, shear)
    """
    axialstress = (effective_internal_pressure) * radius / (2 * thickness) - axialforce / area - (
            (shearforce * z) + moment + (lateralforce * z) * z / 2) / momentofinertia * radius * np.cos(theta)

    hoopstress = effective_internal_pressure * radius / thickness

    shearstress = -((shearforce + lateralforce * z) * (radius**2)) / momentofinertia * np.sin(theta)

    return axialstress, hoopstress, shearstress


def principal_stresses(axialstress: float, hoopstress: float, shearstress: float) -> tuple:
    """
        Finds stress state in the structure as a function of position.
        :param axialstress: longitudinal/meridional stress in the cylinder
        :param hoopstress: circumferential stress in the cylinder
        :param shearstress: circumferential stress in the cylinder
        :return: principal stresses as a function of stress state at a point (p1, p2)
        """
    principalstress1 = (axialstress + hoopstress) / 2 + np.sqrt(
        ((axialstress - hoopstress) / 2) ** 2 + shearstress ** 2)
    principalstress2 = (axialstress + hoopstress) / 2 - np.sqrt(
        ((axialstress - hoopstress) / 2) ** 2 + shearstress ** 2)
    return principalstress1, principalstress2
