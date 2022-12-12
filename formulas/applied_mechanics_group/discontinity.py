"""
DISCONTINUITY EFFECTS AT THE JUNCTIQN OF A PRESSURIZED CYLINDER AND END BULKHEAD - AMG
"""

from math import *
from mpmath import cot


# TODO: RHO IS A POISSIONS RATIO

# def second_curvature_radius(r: float, l: float, h_max: float) -> float: 
#     '''
#     Calculates the second curvature radius (the radsius of a point based on
#     a normal of a tangent line to a point on the doem)
#     :param r: cylinder radius
#     :param l: the horizontal distance between the rotation axis and the point on the dome
#     :param h_max: the maximal height of the dome, measured from its base '''
#     R = sqrt(r**2 - l**2 * ( r**2 - h_max**2 )) / h_max

#     return R

def discontinuity_forces(p, r, R, t_0, t_1, rho_1, psi):
    """
    Calculates maximal stress due to 
    :param p: inner pressure
    :param r: radius of the cylinder
    :param R: radius of the endcap
    :param t_0: thickness of the cylindrical shell
    :param t_1: thickness of the bulkhead
    :param rho_1: poission ratio of the bulkhead 
    :param psi: angle between sth and another thing 
    :return: hoop and longitudinal stresses at a partiucalr point
    """
    phi = asin(r / R) if r != R else pi / 2
    #print("\nphi: ", phi)

    # Verified
    my_lambda = (3 * (1 - rho_1 ** 2) * (R / t_1) ** 2) ** (1 / 4)
    #print("lambda: ", my_lambda)

    # Verified
    delta = my_lambda ** 2 * (r / R) * 2 * (
            (1 + (R / r) ** 0.5 * (t_1 / t_0) ** (3 / 2)) * (1 + (r / R) ** 0.5 * (t_1 / t_0) ** (5 / 2)) - (
            1 - (t_1 / t_0) ** 2) ** 2)
    #print("delta: ", delta)

    # Verified
    miu = my_lambda * (t_1 / t_0) ** 2 * (1 - (r / R) ** 2) ** (1 / 2) * (
            1 + (R / r) ** (1 / 2) * (t_0 / t_1) ** (1 / 2))

    #print("miu: ", miu)

    # Verified
    v = my_lambda ** 2 * (t_1 / t_0) ** 2 * (1 - (r / R) ** 2) ** (1 / 2) * (
            1 + 2 * (R / r) ** (1 / 2) * (t_0 / t_1) ** (1 / 2) + (t_1 / t_0) ** 2)
    #print("v: ", v)

    M_0 = -0.5 * p * r ** 2 * miu / delta
    M_1 = -M_0

    N_0 = -0.5 * p * r * v / delta
    N_1 = -N_0

    M_psi = exp(- my_lambda * psi) * (
            sqrt(2) * M_1 * sin(my_lambda * psi + pi / 4) - R * N_1 / my_lambda * sin(phi) * sin(my_lambda * psi))

    M_theta = rho_1 * M_psi

    Q_psi = exp(- my_lambda * psi) * (2 * my_lambda * M_1 / R * sin(my_lambda * psi)) - sqrt(2) * N_1 * sin(phi) * sin(
        my_lambda * psi - pi / 4)

    N_psi = -Q_psi * cot(phi - psi)

    N_theta = - 2 * my_lambda * exp(- my_lambda * psi) * (
            N_1 * sin(phi) * cos(my_lambda * psi) + sqrt(2) * my_lambda * M_1 / R * (sin(my_lambda * psi) - pi / 4))

    s_meridional = 6 * M_psi / t_1 ** 2 + p * R / (2 * t_1) + N_psi / t_1

    s_hoop = 6 * M_theta / t_1 ** 2 + p * R / (2 * t_1) + N_theta / t_1
    s_max = sqrt(s_meridional ** 2 + s_hoop ** 2 - s_meridional * s_hoop)

    return s_hoop, s_meridional, s_max


# def discontinuity_max_stress_state(r: float, l: float, h_max: float, p: float, t_0: float , t_1: float,
# rho_1:float) -> tuple: R = second_curvature_radius(r, l, h_max)

if __name__ in "__main__":
    p = 60E5
    r = 0.18
    R = 0.181
    t_0 = 0.004
    t_1 = 0.003
    rho_1 = 0.33
    psi = radians(20)
    stress = p * R / (2 * t_1)
    s_max = discontinuity_forces(p, r, R, t_0, t_1, rho_1, psi)

    print(s_max / 1E6)
    print(stress / 1E6)
