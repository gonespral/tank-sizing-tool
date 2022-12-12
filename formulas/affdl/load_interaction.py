"""
Contains formulas related to interaction equations between different loadings.

Assumption: External loads don't cause any additional deformation, such that we can consider these loads alone for the
dome-cylinder connection region.

Source(s): https://engineeringlibrary.org/reference/simple-thin-pressure-vessels-air-force-stress-manual [Section 8.3.1.5.4],
 https://apps.dtic.mil/sti/pdfs/AD0759199.pdf

 In a thin pressure vessel, no stresses other than those tangential to the surface are present at points sufficiently
 removed from discontinuities in the curvature, slope, or thickness of the wall. These tangential or membrane
 stresses are constant throughout the thickness of the shell. At points near discontinuities, such as the junction of
 a cylinder and its head discontinuity, stresses must be superposed upon the membrane stresses in order to obtain the
 total stress.
"""

# Terms used for pressurized cylinders are defined as follows:
#   - R_b: Bending moment ratio
#   - R_c: Compressive load ratio
#   - R_t: Tensile load ratio
#   - R_s: transverse shear load ratio
#   - R_st: Torsional moment ratio

allowed_variation = 0.1


############### Interaction equations ###############

def axial_comp__pure_bending(R_b: float, R_c: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with axial compression and pure bending.
    :param R_b: Bending moment ratio
    :param R_c: Compressive load ratio
    :return: True if in buckling region with axial compression and pure bending, False if not
    """
    A = R_b + R_c
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def axial_comp__pure_bending__transverse_shear(R_b: float, R_c: float, R_s: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with axial compression, pure bending and transverse shear.
    :param R_b: Bending moment ratio
    :param R_c: Compressive load ratio
    :param R_s: Transverse shear load ratio
    :return: True if in pure bending region with transverse shear, False if not
    """
    A = R_c * (R_b ** 3 + R_s ** 3) ** (1 / 3)
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def pure_bending__transverse_shear(R_b: float, R_s: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with pure bending and transverse shear.
    :param R_b: Bending moment ratio
    :param R_s: Transverse shear load ratio
    :return: True if in buckling region with pure bending and transverse shear, False if not
    """
    A = R_s ** 3 + R_b ** 3
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def axial_comp__torsion(R_c: float, R_st: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with axial compression and torsion.
    :param R_c: Compressive load ratio
    :param R_t: Torsional moment ratio
    :return: True if in buckling region with pure axial and torsion, False if not
    """
    A = R_c + R_st ** 2
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def axial_tens__torsion(R_st: float, R_t: float) -> float:
    """
    Determine if the cylinder is in the buckling region with axial tension and torsion.
    :param R_st: Torsional moment ratio
    :param R_t: Tensile load ratio
    :return: True if in buckling region with pure axial and torsion, False if not
    """
    if not R_t < 0.8:
        raise ValueError("Tensile load ratio must be less than 0.8")
    A = R_st ** 3 - R_t
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def pure_bending__torsion(R_b: float, R_st: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with pure bending and torsion.
    :param R_b: Bending moment ratio
    :param R_st: Torsional moment ratio
    :return: True if in buckling region with pure bending and torsion, False if not
    """
    A = R_b ** 1.5 + R_st ** 2
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def pure_bending__torsion__transverse_shear(R_b: float, R_s: float, R_st: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with pure bending, torsion and transverse shear.
    :param R_b: Bending moment ratio
    :param R_s: Transverse shear load ratio
    :param R_st: Torsional moment ratio
    :return: True if in buckling region with pure bending, torsion and transverse shear, False if not
    """
    p = 2.25  # 1.5 <= p <= 3.0
    q = 2.5  # 2.0 <= q <= 3.0
    A = R_b ** p + (R_s + R_st) ** q
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def axial_comp__pure_bending__transverse_shear__torsion(R_c: float, R_st: float, R_s: float, R_b: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with axial compression, pure bending, transverse shear and torsion.
    :param R_c: Compressive load ratio
    :param R_st: Torsional moment ratio
    :param R_s: Transverse shear load ratio
    :param R_b: Bending moment ratio
    :return: True if in buckling region with pure axial, pure bending, transverse shear and torsion, False if not
    """
    A = R_c + R_st ** 2 + (R_s ** 3 + R_b ** 3) ** (1 / 3)
    if A < 1 - allowed_variation:
        return False
    else:
        return True


def axial_load__pure_bending__torsion(R_s: float, R_b: float, R_st: float) -> bool:
    """
    Determine if the cylinder is in the buckling region with axial compression, pure bending and torsion.
    :param R_s: Transverse shear load ratio
    :param R_b: Bending moment ratio
    :param R_st: Torsional moment ratio
    :return: True if in buckling region with pure axial, pure bending and torsion, False if not
    """
    # TODO: unsure about this one
    pass
