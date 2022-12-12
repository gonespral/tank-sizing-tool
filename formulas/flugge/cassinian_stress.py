"""
This function finds the meridional/longitudinal and circumferential/hoop stress due to internal pressure in a cassinian dome endcap. Discontinuity effects and external loading cases are not considered.

Source: "Stresses in Shells" by Flugge
"""
from databases.fos import *

def stress_position(effective_internal_pressure: float, z: float, x: float, radius: float, n: float,
                    thickness: float) -> tuple:
    """
    Finds stress state in a two-dimensional model of a cassinian dome with affine curvature. The equation of the dome is:
    (z**2+(n**2)*(x**2))**2+2*(radius**2)*(z**2-(n**2)*(x**2))=3*(radius**4)

    :param effective_internal_pressure: effective internal pressure in the dome (basic pressure + hydrostatic pressure) in Pa
    :param z: z-coordinate of the position on the dome, ranges from -radius to +radius, the z-axis is any axis perpendicular to the x-axis in m
    :param x: x-coordinate of the position on the dome, ranges from 0 at cylinder joint to the dome height at the vertex, the x-axis is the vertical axis parallel to the cylinder length in m
    :param radius: internal (can also be median) radius of the cylinder in m
    :param n: affine transformation factor, ranges from 1 to infinity, should be kept under 1.9 to avoid compressive stresses
    :param thickness: dome thickness in m

    :return: stress state as a function of position in the cylinder (meridional, circumferential)
    """
    #Fos application
    effective_internal_pressure = Coeff_A_Pressure * effective_internal_pressure * Coeff_B


    meridional_stress = effective_internal_pressure * radius / thickness * (
                ((z ** 2) * (radius ** 2 + (n ** 2) * (x ** 2)) + (n ** 4) * (x ** 2) * (radius ** 2 - z ** 2)) ** (
                    1 / 2)) / (radius ** 2 + z ** 2 + (n ** 2) * (x ** 2))

    circumferential_stress = meridional_stress * (2 - (
                (3 * (n ** 2) * (radius ** 4) * (radius ** 2 - z ** 2 + (n ** 2) * (x ** 2))) / (
                    (radius ** 2 + z ** 2 + (n ** 2) * (x ** 2)) * (
            ((z ** 2) * (radius ** 2 + (n ** 2) * (x ** 2)) + (n ** 4) * (x ** 2) * (radius ** 2 - z ** 2))))))

    return meridional_stress, circumferential_stress
