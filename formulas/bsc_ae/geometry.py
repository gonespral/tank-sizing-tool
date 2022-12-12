"""
Contains formulas for general geometric calculations.

Source(s): my mind, my brain, my heart, my soul, my body, my blood, my sweat, my tears, my life, my death, my everything
"""

import numpy as np
from numpy.polynomial.polynomial import Polynomial
from scipy.integrate import trapz


# CUBOID #############################################################

def rectangle_I(width: float, height: float) -> float:
    """
    Calculates the area moment of inertia of a rectangular section.
    :param width: width of the rectangle
    :param height: height of the rectangle
    :return: moment of inertia of the rectangular section
    """
    return (width * (height ** 3)) / 12


def rectangle_J(width: float, height: float) -> float:
    """
    Calculates the torsional constant of a rectangle.
    :param width: width of the rectangle
    :param height: height of the rectangle
    :return: torsion constant of the rectangular section
    """
    a = max(width, height)  # Long side
    b = min(width, height)  # Short side

    A = a * (b ** 3)
    B = 1 - (b ** 4) / (12 * (a ** 4))
    C = 0.21 * (b / a)

    J = A * ((1 / 3) - C * B)
    return J


# CYLINDER ###########################################################

def cylindrical_shell_I(radius: float, thickness: float) -> float:
    """
    Calculates the moment of inertia of a cylindrical shell.
    :param radius: outer radius of the cylindrical shell
    :param thickness: thickness of the cylindrical shell
    :return: moment of inertia of the cylindrical shell
    """
    return (1 / 4) * np.pi * (radius ** 4) - (1 / 4) * np.pi * ((radius - thickness) ** 4)


def cylindrical_shell_A(radius: float, thickness: float) -> float:
    """
    Calculates the sectional area of a cylindrical shell.
    :param radius:outer radius of the cylindrical shell
    :param thickness: thickness of the cylindrical shell
    :return: sectional area of the cylindrical shell
    """
    return np.pi * (radius ** 2) - np.pi * ((radius - thickness) ** 2)


def cylinder_V(radius: float, height: float) -> float:
    """
    Calculates the volume of a cylinder.
    :param radius: outer radius of the cylinder
    :param height: height of the cylinder
    :return: cylinder volume
    """
    return np.pi * (radius ** 2) * height


def cylinder_h(volume: float, radius: float) -> float:
    """
    Calculates the height of a cylinder.
    :param volume: volume of the cylinder
    :param radius: outer radius of the cylinder
    :return: height of the cylinder
    """
    return volume / (np.pi * (radius ** 2))


def cylinder_r(volume: float, height: float) -> float:
    """
    Calculates the radius of a cylinder.
    :param volume: volume of the cylinder
    :param height: height of the cylinder
    :return: radius of the cylinder
    """
    return np.sqrt(volume / (np.pi * height))


# HEMISPHERE ########################################################

def hemisphere_V(radius: float) -> float:
    """
    Calculates the volume of a hemispherical dome.
    :param radius: radius of the hemispherical dome
    :return: volume of the hemispherical dome
    """
    # TODO: Validate this
    return (2 / 3) * np.pi * (radius ** 3)


# SEMI-ELLIPSOID ####################################################

def semi_ellipsoid_V(radius: float, height: float) -> float:
    """
    Calculates the volume of an elliptical dome. Assumes semimajor axis = semiminor axis = radius.
    :param radius: radius of the elliptical dome
    :param height: height of the elliptical dome
    :return: volume of the elliptical dome
    """

    return ((4 / 3) * np.pi * radius ** 2 * height) / 2


def semi_ellipsoid_h(volume: float, radius: float) -> float:
    """
    Calculates the height of an elliptical dome. Assumes semimajor axis = semiminor axis = radius.
    :param volume: volume of the elliptical dome
    :param radius: radius of the elliptical dome
    :return: height of the elliptical dome
    """

    return (2 * volume) / ((4 / 3) * np.pi * radius ** 2)


# SPHERICAL SECTION #################################################

def spherical_section_V(height: float, radius: float = None, meridional_radius: float = None) -> float:
    """
    Calculates the volume of a spherical section. You must supply the height plus either the radius or the meridional
    radius.
    Source: https://www.cuemath.com/measurement/volume-of-a-section-of-a-sphere/
    :param height: height of the spherical section from sectional cut
    :param radius: radius of the spherical section
    :param meridional_radius: meridional radius of the spherical section
    :return: volume of the spherical section
    """
    # TODO: Validate this
    if radius is None and meridional_radius is None:
        raise ValueError("You must supply either the radius or the meridional radius")
    elif radius is not None and meridional_radius is not None:
        raise ValueError("You must supply either the radius or the meridional radius, not both")

    if meridional_radius is not None:
        # Use meridional radius
        return (1 / 3) * np.pi * (height ** 2) * (3 * meridional_radius - height)
    elif radius is not None:
        # Use radius
        return (1 / 6) * np.pi * height * (3 * radius ** 2 - height ** 2)


def spherical_section_h(volume: float, radius: float = None, meridional_radius: float = None) -> float:
    """
    Calculates the height of a spherical section. You must supply the volume plus either the radius or the meridional
    radius.
    Source: https://www.cuemath.com/measurement/volume-of-a-section-of-a-sphere/
    :param volume: volume of the spherical section
    :param radius: radius of the spherical section
    :param meridional_radius: meridional radius of the spherical section
    :return: height of the spherical section
    """
    # TODO: Validate this
    if radius is None and meridional_radius is None:
        raise ValueError("You must supply either the radius or the meridional radius")
    elif radius is not None and meridional_radius is not None:
        raise ValueError("You must supply either the radius or the meridional radius, not both")

    if meridional_radius is not None:
        # Use meridional radius
        # Find height using np polynomial
        A = 1
        B = -3 * meridional_radius
        C = 0
        D = (3 * volume) / np.pi
        polynomial = Polynomial([D, C, B, A],
                                symbol="h")  # See https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.Polynomial.html
    elif radius is not None:
        # Use radius
        # Find height using np polynomial
        A = 1
        B = 0
        C = 3 * radius ** 2
        D = -(6 * volume) / np.pi
        polynomial = Polynomial([D, C, B, A],
                                symbol="h")  # See https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.Polynomial.html

    # Find roots
    roots = polynomial.roots()
    # Find real roots
    real_roots = roots[np.isreal(roots)]
    # Remove imaginary part
    real_roots = real_roots.real
    # Find positive real roots
    positive_real_roots = real_roots[real_roots > 0]
    # Find smallest positive real root
    smallest_positive_real_root = np.min(positive_real_roots)
    return smallest_positive_real_root


# CASSINIAN ##########################################################################

def cassinian_h(radius: float, n: float) -> float:
    """
    Calculates the height of a cassinian dome defined by (z**2+(n**2)*(x**2))**2+2*(radius**2)*(z**2-(n**2)*(x**2))=3*(radius**4).
    :param radius: radius of the cassinian oval
    :param n: affine transformation factor, ranges from 1 to 1.9 to avoid compressive stresses
    :return: height of the cassinian dome
    """
    # TODO: Validate this
    return radius / n * (3) ** (1 / 2)


def cassinian_V(radius: float, n: float) -> float:
    """
    Calculates the volume of a cassinian dome defined by (z**2+(n**2)*(x**2))**2+2*(radius**2)*(z**2-(n**2)*(x**2))=3*(radius**4).
    :param radius: radius of the cassinian oval
    :param n: affine transformation factor, ranges from 1 to 1.9 to avoid compressive stresses
    :return: volume of the cassinian dome
    """
    # TODO: Validate this
    x = np.linspace(0, radius / n * (3) ** (1 / 2), 500)
    f = (radius ** 2) * (-1 - (n ** 2) * (x / radius) ** 2 + (4 + 4 * (n ** 2) * (x / radius) ** 2) ** (1 / 2))
    return np.pi * trapz(f, x)


def cassinian_z_d_to_section_radius(z_d: float, radius: float, n: float) -> float:
    """
    Calculates the radius as a function of height of a cassinian dome defined by (z**2+(n**2)*(x**2))**2+2*(radius**2)*(z**2-(n**2)*(x**2))=3*(radius**4).
    :param radius: radius of the cassinian oval
    :param n: affine transformation factor, ranges from 1 to 1.9 to avoid compressive stresses
    :return: radius of the cassinian dome at a certain height
    """
    # TODO: Validate this
    x = z_d
    r = radius * (-1 - (n ** 2) * (x / radius) ** 2 + (4 + 4 * (n ** 2) * (x / radius) ** 2) ** (1 / 2)) ** (1 / 2)
    return r
