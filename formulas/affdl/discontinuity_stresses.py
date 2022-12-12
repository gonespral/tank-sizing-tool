"""
Contains formulas related to discontinuity stresses generated at the interface between the dome and cylindrical section
of the tank. Assumes thin wall.

Assumption: External loads don't cause any additional deformation, such that we can consider these loads alone for the
dome-cylinder connection region.

Source(s): https://engineeringlibrary.org/reference/simple-thin-pressure-vessels-air-force-stress-manual [Section 8.3.1.2.2],
 https://apps.dtic.mil/sti/pdfs/AD0759199.pdf

 In a thin pressure vessel, no stresses other than those tangential to the surface are present at points sufficiently
 removed from discontinuities in the curvature, slope, or thickness of the wall. These tangential or membrane
 stresses are constant throughout the thickness of the shell. At points near discontinuities, such as the junction of
 a cylinder and its head discontinuity, stresses must be superposed upon the membrane stresses in order to obtain the
 total stress.
"""

import os
import sys

from typing import Tuple, Any, Dict
import numpy as np

# Displacement of heads (here referred to as domes) and cylinders due to internal pressure
# w_m1 refers to displacement of the dome due to internal pressure
# w_m2 refers to displacement of the cylinder due to internal pressure

# TODO: Validate these formulas


def w_m1_hemispherical_dome(pressure: float, radius: float, youngs_modulus: float, thickness: float,
                            poisson_ratio: float) -> float:
    """
    Finds outward displacement of a hemispherical dome to internal pressure.
    :param pressure: Internal pressure.
    :param radius: Radius of the dome (section cut to edge).
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the dome.
    :param poisson_ratio: Poisson's ratio of the material.
    :return wm1: Displacement of the dome.
    """
    A = (pressure * (radius ** 2)) / (2 * youngs_modulus * thickness)
    B = (1 - poisson_ratio)
    return A * B


def w_m1_figure_of_revolution(pressure: float, radius: float, meridional_radius: float, youngs_modulus: float,
                              thickness: float, poisson_ratio: float) -> float:
    """
    Finds outward displacement of a figure of revolution to internal pressure. Note that for semi-ellipsoidal domes,
    the meridional radius varies along the dome, and as a result, so does wm1.
    :param pressure: Internal pressure.
    :param radius: Radius of the figure of revolution (section cut to edge).
    :param meridional_radius: Meridional radius of the figure of revolution.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the figure of revolution.
    :param poisson_ratio: Poisson's ratio of the material.
    :return wm1: Displacement of the figure of revolution.
    """
    A = (pressure * (radius ** 2)) / (2 * youngs_modulus * thickness)
    B = 2 - poisson_ratio - (radius / meridional_radius)
    return A * B


def w_m2_cylinder(pressure: float, radius: float, youngs_modulus: float, thickness: float,
                  poisson_ratio: float) -> float:
    """
    Finds outward displacement of a cylinder to internal pressure.
    :param pressure: Internal pressure.
    :param radius: Radius of the cylinder.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the cylinder.
    :param poisson_ratio: Poisson's ratio of the material.
    :return wm2: Displacement of the cylinder.
    """
    A = (pressure * (radius ** 2)) / (youngs_modulus * thickness)
    B = 1 - (poisson_ratio / 2)
    return A * B


# Radial and angular displacements due to edge loading. Applicable to semi-ellipsoidal, conical, dished or hemispherical
# domes.
# w_1 refers to radial displacement of the dome due to edge loading
# w_2 refers to radial displacement of the cylinder due to edge loading
# theta_1 refers to angular displacement of the dome due to edge loading
# theta_2 refers to angular displacement of the cylinder due to edge loading

def _lambda(poisson_ratio: float, radius: float, thickness: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param radius: Radius of the dome (section cut to edge).
    :param thickness: Thickness of the dome (for lambda 1) or of the cylinder (for lambda 2).
    :return: lambda_1 or lambda_2
    """
    A = 3 * (1 - poisson_ratio ** 2)
    B = (radius ** 2) * (thickness ** 2)
    return (A / B) ** (1 / 4)


# The following set of equations is created from the formulas in the source, with a function which calculates each
# coefficient:
# With: w_1 - w_2 = w_m1 - w_m2 and theta_1 - theta_2 = 0
# (A_1 - A_2) * Q_0 + (B_1 - B_2) * M_0 = w_m1 - w_m2
# (C_1 - C_2) * Q_0 + (D_1 - D_2) * M_0 = 0
# Then, numpy linalg is used to solve for Q_0 and M_0, with w_1 and w_2 are calculated from the functions above.
# Note: A lot of the following code is repeated (not reused), but for the sake of clarity, it is kept this way.

def _A_1(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_1: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the dome.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_1: lambda_1 (for the dome).
    :return: A_1
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = -(np.sin(phi_0) ** 2) / (2 * lambda_1 ** 3)
    return A * B


def _A_2(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_2: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the cylinder.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_2: lambda_2 (for the cylinder).
    :return: A_2
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = (np.sin(phi_0) ** 2) / (2 * lambda_2 ** 3)
    return A * B


def _B_1(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_1: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the dome.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_1: lambda_1 (for the dome).
    :return: B_1
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = -np.sin(phi_0) / (2 * lambda_1 ** 2)
    return A * B


def _B_2(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_2: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the cylinder.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_2: lambda_2 (for the cylinder).
    :return: B_2
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = -np.sin(phi_0) / (2 * lambda_2 ** 2)
    return A * B


def _C_1(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_1: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the dome.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_1: lambda_1 (for the dome).
    :return: C_1
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = np.sin(phi_0) / (2 * lambda_1 ** 2)
    return A * B


def _C_2(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_2: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the cylinder.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_2: lambda_2 (for the cylinder).
    :return: C_2
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = np.sin(phi_0) / (2 * lambda_2 ** 2)
    return A * B


def _D_1(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_1: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the dome.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_1: lambda_1 (for the dome).
    :return: D_1
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = -1 / lambda_1
    return A * B


def _D_2(poisson_ratio: float, youngs_modulus: float, thickness: float, phi_0: float, lambda_2: float) -> float:
    """
    Intermediate step.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness: Thickness of the cylinder.
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :param lambda_2: lambda_2 (for the cylinder).
    :return: D_2
    """
    A = (12 * (1 - poisson_ratio ** 2)) / (youngs_modulus * thickness ** 3)
    B = -1 / lambda_2
    return A * B


def _find_coefficients(poisson_ratio: float, youngs_modulus: float, thickness_dome: float, thickness_cylinder: float,
                       lambda_1: float, lambda_2: float,
                       phi_0: float) -> Tuple[float, float, float, float]:
    """
    Find the coefficients for the equations.
    :param poisson_ratio: Poisson's ratio of the material.
    :param youngs_modulus: Young's modulus of the material.
    :param thickness_dome: Thickness of the dome.
    :param thickness_cylinder: Thickness of the cylinder.
    :param lambda_1: Lambda 1 (for the dome).
    :param lambda_2: Lambda 2 (for the cylinder).
    :param phi_0: Angle between the edge loading and the cylinder axis.
    :return:
    """
    A_1 = _A_1(poisson_ratio, youngs_modulus, thickness_dome, phi_0, lambda_1)
    A_2 = _A_2(poisson_ratio, youngs_modulus, thickness_cylinder, phi_0, lambda_2)
    B_1 = _B_1(poisson_ratio, youngs_modulus, thickness_dome, phi_0, lambda_1)
    B_2 = _B_2(poisson_ratio, youngs_modulus, thickness_cylinder, phi_0, lambda_2)
    C_1 = _C_1(poisson_ratio, youngs_modulus, thickness_dome, phi_0, lambda_1)
    C_2 = _C_2(poisson_ratio, youngs_modulus, thickness_cylinder, phi_0, lambda_2)
    D_1 = _D_1(poisson_ratio, youngs_modulus, thickness_dome, phi_0, lambda_1)
    D_2 = _D_2(poisson_ratio, youngs_modulus, thickness_cylinder, phi_0, lambda_2)

    A = A_1 - A_2
    B = B_1 - B_2
    C = C_1 - C_2
    D = D_1 - D_2

    return A, B, C, D


def _find_force_moment(coeff: Tuple[float, float, float, float], w_m1: float, w_m2: float) -> tuple[float, float]:
    """
    Finds the coefficients of system of equations.
    :param w_m1: Deflection at the dome.
    :param w_m2: Deflection at the cylinder.
    :return: Q_0, M_0
    """

    LHS = np.array([[coeff[0], coeff[1]], [coeff[2], coeff[3]]])
    RHS = np.array([w_m1 - w_m2, 0])
    # Solve the system of equations with least squares (linalg.solve gives error if matrix is singular).
    Q_0, M_0 = np.linalg.lstsq(LHS, RHS, rcond=None)[0]
    return Q_0, M_0


def _hoop_normal_stress(x: float, lambda_: float, radius: float, thickness: float, R_t: float, Q_0: float,
                        M_0: float) -> float:
    """
    Calculates the hoop normal stress. (F_t)
    :param x: x-coordinate along the meridian from the junction.
    :param lambda_: lambda.
    :param radius: Radius.
    :param thickness: Thickness.
    :param R_t: Tangential radius.
    :param Q_0: Q_0 (intermediate step).
    :param M_0: M_0 (intermediate step).
    :return: Hoop normal stress.
    """
    A = (2 * lambda_ * radius) / thickness
    B = radius / R_t
    C = np.exp(-lambda_ * x)
    D = Q_0 * np.cos(np.deg2rad(lambda_ * x)) - (lambda_ * M_0 * (np.cos(np.deg2rad(lambda_ * x)) - np.sin(np.deg2rad(lambda_ * x))))
    return A * B * C * D


def _meridional_shear_stress(x: float, lambda_: float, radius: float, thickness: float, R_t: float, Q_0: float,
                             M_0: float) -> float:
    """
    Calculates the meridional shear stress. (F_smer)
    :param x: x-coordinate along the meridian from the junction.
    :param lambda_: lambda.
    :param radius: Radius.
    :param thickness: Thickness.
    :param R_t: Tangential radius.
    :param Q_0: Q_0 (intermediate step).
    :param M_0: M_0 (intermediate step).
    :return: Meridional shear stress.
    """
    A = radius / (thickness * R_t)
    B = (radius / R_t) ** (1 / 2)
    C = np.exp(-lambda_ * x)
    D = Q_0 * (np.cos(np.deg2rad(lambda_ * x)) - np.sin(np.deg2rad(lambda_ * x))) + (2 * lambda_ * M_0 * np.sin(np.deg2rad(lambda_ * x)))
    return A * B * C * D


def _maximum_meridional_bending_stress(x: float, lambda_: float, radius: float, thickness: float, R_t: float, Q_0: float,
                                       M_0: float) -> float:
    """
    Calculates the maximum meridional bending stress. Maximum is taken since this varies throughout the thickness. (F_bmer)
    :param x: x-coordinate along the meridian from the junction.
    :param lambda_: lambda.
    :param radius: Radius.
    :param thickness: Thickness.
    :param R_t: Tangential radius.
    :param Q_0: Q_0 (intermediate step).
    :param M_0: M_0 (intermediate step).
    :return: Maximum meridional shear stress.
    """
    A = (6 * radius) / ((thickness ** 2) * R_t)
    B = 1 / lambda_
    C = np.exp(-lambda_ * x)
    D = (-Q_0 * np.sin(np.deg2rad(lambda_ * x))) + (lambda_ * M_0 * (np.cos(np.deg2rad(lambda_ * x)) + np.sin(np.deg2rad(lambda_ * x))))
    return A * B * C * D


def _maximum_hoop_bending_stress(x: float, lambda_: float, radius: float, thickness: float, R_t: float, Q_0: float,
                                 M_0: float, youngs_modulus: float, poisson_ratio: float, phi_0: float) -> float:
    """
    Calculates the maximum hoop bending stress. Maximum is take since this varies throughout the thickness. (F_bt)
    :param x: x-coordinate along the meridian from the junction.
    :param lambda_: lambda.
    :param radius: Radius.
    :param thickness: Thickness.
    :param R_t: Tangential radius.
    :param Q_0: Q_0 (intermediate step).
    :param M_0: M_0 (intermediate step).
    :param youngs_modulus: Young's modulus.
    :param poisson_ratio: Poisson's ratio.
    :param phi_0: Initial angle.
    :return: Maximum hoop bending stress.
    """
    F_bmer = _maximum_meridional_bending_stress(x, lambda_, radius, thickness, R_t, Q_0, M_0)
    A = poisson_ratio * F_bmer
    B = (youngs_modulus * thickness) / 2
    C = 1 / (R_t * np.tan(np.deg2rad(phi_0)))
    D = (6 * (1 - (poisson_ratio ** 2))) / (youngs_modulus * (thickness ** 3) * (lambda_ ** 2))
    E = np.exp(-lambda_ * x)
    F = Q_0 * (np.cos(np.deg2rad(lambda_ * x)) + np.sin(np.deg2rad(lambda_ * x))) - (2 * lambda_ * M_0 * np.cos(np.deg2rad(lambda_ * x)))
    return A + (B * C * D * E * F)


def find_stress_state(pressure: float, phi_0: float, thickness_dome: float,
                      thickness_cylinder: float, poisson_ratio: float, youngs_modulus: float,
                      radius: float, x: float, dome_type: str, find_in: str, R_t: float, R_m: float = None) -> dict:
    """
    Calculates the stress state at a point on the meridian of a spherical dome. NOTE: Input and output parameters are
    in SI units.
    :param pressure: Pressure (Pa).
    :param phi_0: Angle (rad).
    :param thickness_dome: Dome thickness (m).
    :param thickness_cylinder: Cylinder thickness (m).
    :param poisson_ratio: Poisson's ratio (unitless).
    :param youngs_modulus: Young's modulus (Pa).
    :param radius: Radius (m).
    :param x: x-coordinate along the meridian from the junction (m).
    :param dome_type: Type of dome.
    :param find_in: Find in.
    :param R_t: Tangential radius (m).
    :param R_m: Meridional radius (m). Only required if find_in is 'meridional'.
    :return: Stress state at a point on the meridian of a spherical dome (Pa).
    """

    # Convert Pa to psi
    pressure = pressure / 6894.75729
    youngs_modulus = youngs_modulus / 6894.75729

    # Convert from m to in
    radius = radius * 39.3701
    thickness_dome = thickness_dome * 39.3701
    thickness_cylinder = thickness_cylinder * 39.3701
    x = x * 39.3701
    R_t = R_t * 39.3701
    if R_m is not None:
        R_m = R_m * 39.3701

    dome_types = ["hemispherical", "figure_of_revolution"]
    find_ins = ["cylinder", "dome"]

    if dome_type == "semi-ellipsoidal":
        dome_type = "figure_of_revolution"
        print("Semi-ellipsoidal dome -> figure of revolution")
    elif dome_type == "dished":
        dome_type = "figure_of_revolution"
        print("Dished dome -> figure of revolution")

    if dome_type not in dome_types:
        raise ValueError(f"Dome type must be one of {dome_types}")
    if dome_type == "figure_of_revolution" and R_m is None:
        raise ValueError("R_m must be specified for a figure of revolution dome.")
    if find_in not in find_ins:
        raise ValueError(f"Find_in must be one of {find_ins}")

    if dome_type == "hemispherical":
        w_m1 = w_m1_hemispherical_dome(radius=radius,
                                       pressure=pressure,
                                       youngs_modulus=youngs_modulus,
                                       poisson_ratio=poisson_ratio,
                                       thickness=thickness_dome)
    elif dome_type == "figure_of_revolution":
        w_m1 = w_m1_figure_of_revolution(radius=radius,
                                         meridional_radius=R_m,
                                         pressure=pressure,
                                         youngs_modulus=youngs_modulus,
                                         poisson_ratio=poisson_ratio,
                                         thickness=thickness_dome)

    w_m2 = w_m2_cylinder(radius=radius,
                         youngs_modulus=youngs_modulus,
                         poisson_ratio=poisson_ratio,
                         thickness=thickness_cylinder,
                         pressure=pressure)

    lambda_1 = _lambda(poisson_ratio=poisson_ratio,
                       radius=radius,
                       thickness=thickness_dome)

    lambda_2 = _lambda(poisson_ratio=poisson_ratio,
                       radius=radius,
                       thickness=thickness_cylinder)

    coeff = _find_coefficients(poisson_ratio=poisson_ratio,
                               youngs_modulus=youngs_modulus,
                               thickness_dome=thickness_dome,
                               thickness_cylinder=thickness_cylinder,
                               lambda_1=lambda_1,
                               lambda_2=lambda_2,
                               phi_0=phi_0)

    Q_0, M_0 = _find_force_moment(coeff=coeff,
                                  w_m1=w_m1,
                                  w_m2=w_m2)

    if find_in == "cylinder":
        lambda_ = lambda_2
        thickness = thickness_cylinder

    elif find_in == "dome":
        lambda_ = lambda_1
        thickness = thickness_dome

    hoop_normal_stress = _hoop_normal_stress(x=x,
                                             lambda_=lambda_,
                                             radius=radius,
                                             thickness=thickness,
                                             R_t=R_t,
                                             Q_0=Q_0,
                                             M_0=M_0)

    meridional_shear_stress = _meridional_shear_stress(x=x,
                                                       lambda_=lambda_,
                                                       radius=radius,
                                                       thickness=thickness,
                                                       R_t=R_t,
                                                       Q_0=Q_0,
                                                       M_0=M_0)

    maximum_meridional_bending_stress = _maximum_meridional_bending_stress(x=x, lambda_=lambda_,
                                                                           radius=radius,
                                                                           thickness=thickness,
                                                                           R_t=R_t,
                                                                           Q_0=Q_0,
                                                                           M_0=M_0)

    maximum_hoop_bending_stress = _maximum_hoop_bending_stress(x=x,
                                                               lambda_=lambda_,
                                                               radius=radius,
                                                               thickness=thickness,
                                                               R_t=R_t,
                                                               Q_0=Q_0,
                                                               M_0=M_0,
                                                               youngs_modulus=youngs_modulus,
                                                               poisson_ratio=poisson_ratio,
                                                               phi_0=phi_0)

    # Convert pressures from psi to Pa
    hoop_normal_stress = hoop_normal_stress * 6894.75729
    meridional_shear_stress = meridional_shear_stress * 6894.75729
    maximum_meridional_bending_stress = maximum_meridional_bending_stress * 6894.75729
    maximum_hoop_bending_stress = maximum_hoop_bending_stress * 6894.75729
    resultant = hoop_normal_stress + meridional_shear_stress + maximum_meridional_bending_stress + maximum_hoop_bending_stress

    # Convert moment from in-lb to N-m
    M_0 = M_0 * 0.112984829

    # Convert force from lb to N
    Q_0 = Q_0 * 4.448221615

    # Convert from in to m
    w_m1 = w_m1 * 0.0254
    w_m2 = w_m2 * 0.0254

    return {"hoop_normal_stress": hoop_normal_stress,
            "meridional_shear_stress": meridional_shear_stress,
            "maximum_meridional_bending_stress": maximum_meridional_bending_stress,
            "maximum_hoop_bending_stress": maximum_hoop_bending_stress,
            "Q_0": Q_0,
            "M_0": M_0,
            "w_m1": w_m1,
            "w_m2": w_m2,
            "lambda_1": lambda_1,
            "lambda_2": lambda_2,
            "resultant": resultant}
