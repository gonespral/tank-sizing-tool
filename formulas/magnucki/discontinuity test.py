import numpy as np
import matplotlib.pyplot as plt

def discontinuity_stress(p0: float, t1: float, t2: float, a: float, b: float, E: float, nu: float,
                         theta: float) -> tuple:
    """
    Finds stress state in at a singular point in the ellipse.

    :param p0: effective internal pressure in the dome (basic pressure + hydrostatic pressure) in Pa
    :param a: major axis of the ellipse, in this case equal to the radius (probably internal, but no practical difference from median) in m
    :param b: minor axis of the ellipse, in this case equal to the height of the semi-ellipse (again, probably internal, but no practical difference from median) in m
    :param t1: constant dome thickness in m
    :param t2: cylinder thickness in m
    :param E: Young's modulus of both ellipse and cylinder
    :param nu: Poisson's ratio
    :param theta: angle between a line tangent to the surface at the point in question and the horizontal, for clearer explanation check the paper
    :return: stress state as a function of position in the cylinder (meridional, circumferential)
    """

    x1 = t2 / t1  # ratio of wall thicknesses
    x2 = a / t2  # characteristic ratio of cylinder
    beta = b / a  # relative convexity of dome
    xi = beta ** (-2) - 1  # dimensionless parameter
    cv = (3 * (1 - nu ** 2)) ** (1 / 4)


    alpha1_11 = 2 / E * cv * (a / t1) ** (
                3 / 2)  # random factors that go into the solving matrix for the discontinuity shear and moment
    alpha1_12 = -2 / E / t1 * (cv ** 2) * (a / t1)
    alpha1_22 = 4 / E / (t1 ** 2) * (cv ** 3) * (a / t1) ** (1 / 2)
    alpha2_11 = -2 / E * cv * (a / t2) ** (3 / 2)
    alpha2_12 = -2 / E / t2 * (cv ** 2) * (a / t2)
    alpha2_22 = -4 / E / (t2 ** 2) * (cv ** 3) * (a / t2) ** (
                1 / 2)  # paper says t1 for the first t2, but it might be an error

    u0_star = (2 - nu - (a / b) ** 2) * p0 * (
                a ** 2) / 2 / E / t1  # base displacements, go into the displacement vector
    w0_star = (2 - nu) * p0 * (a ** 2) / 2 / E / t2

    A = np.array([[alpha2_11 - alpha1_11, alpha2_12 - alpha1_12],
                  [alpha2_12 - alpha1_12,
                   alpha2_22 - alpha1_22]])  # solving matrix, signifies displacement compatibility
    displacement = np.array([u0_star - w0_star, 0])

    solution = np.linalg.solve(A, displacement)  # solution vector, contains shear force and moment
    Q0 = solution[0]  # discontinuity shear force
    M0 = solution[1]  # discontinuity moment


    sigmaM = 6 * M0 / t1 ** 2
    sigmaQ = 2 * cv * (x1 * x2) ** (1 / 2) * Q0 / t1

    x = np.linspace(theta, np.pi / 2, 100)
    f = 1 / (1 + xi * (np.sin(x)) ** 2) ** (5 / 4)
    psi = cv * ((x1 * x2 / beta) ** (1 / 2)) * np.trapz(f, x)

    fthetaBD = sigmaM*(np.cos(psi) + np.sin(psi)) - 3 / (cv ** 2) * sigmaQ * np.sin(psi)
    sigmatheta1 = x1 * x2 / 2 / beta / (1 + xi * (np.sin(theta)) ** 2) ** (1 / 2) - fthetaBD * np.exp(-psi)
    sigmatheta2 = x1 * x2 / 2 / beta / (1 + xi * (np.sin(theta)) ** 2) ** (1 / 2) + fthetaBD * np.exp(-psi)
    sigmatheta = max(sigmatheta1, sigmatheta2)  # dimensionless meridional factor

    fphiBD1 = sigmaM * (((cv ** 2) / 3 - nu) * np.cos(psi) - (cv ** 2) / 3 + nu) * np.cos(psi) - sigmatheta * (
                np.cos(psi) - 3 * nu / (cv ** 2) * np.sin(psi))
    fphiBD2 = sigmaM * (((cv ** 2) / 3 + nu) * np.cos(psi) - (cv ** 2) / 3 - nu) * np.cos(psi) - sigmatheta * (
                np.cos(psi) + 3 * nu / (cv ** 2) * np.sin(psi))
    sigmaphi1 = x1 * x2 / 2 / beta * (1 - xi * (np.sin(theta)) ** 2) / (1 + xi * (np.sin(theta)) ** 2) ** (
                1 / 2) - fphiBD1 * np.exp(-psi)
    sigmaphi2 = x1 * x2 / 2 / beta * (1 - xi * (np.sin(theta)) ** 2) / (1 + xi * (np.sin(theta)) ** 2) ** (
                1 / 2) - fphiBD2 * np.exp(-psi)
    sigmaphi = max(sigmaphi1, sigmaphi2)  # dimensionless circumferential factor

    meridional_stress = sigmatheta * p0
    circumferential_stress = sigmaphi * p0

    return meridional_stress, circumferential_stress


a=0.178
b=0.11
t1=0.005
t2=0.006
E=70*10**9
p0=45
nu=	0.33

stress_array=[]

for theta in np.linspace(0, np.pi/2, 100):
    stress_array.append(discontinuity_stress(p0, t1, t2, a, b, E, nu,theta))

plt.plot(stress_array[0], label='meridional stress')
plt.plot(stress_array[1], label='circumferential stress')
plt.xlabel('theta'); plt.ylabel('stress in Pa')
plt.legend()
plt.show()