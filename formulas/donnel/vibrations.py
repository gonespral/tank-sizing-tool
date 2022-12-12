"""
[File description]

Source(s): [sources]
"""
# TODO: Fill in file header
# TODO: Add sources per function

from math import sqrt, pi


# Disclaimer: This might not be correct
def find_1st_bending_freq_cylindrical_shell(r_o: float, L: float, t: float, v: float, E: float, rho: float) -> int:
    """
    Calculates the 1st bending frequency of the 
    :param r_o: outer radius of the cylindrical shell
    :param L: length of the cylindrical section
    :param t: thickness of the cylindrical shell
    :param v: Poisson's ratio
    :param E: modulus of elasticity
    :param rho: density
    :param Nt: circumferential load per unit length, often results from differential pressure p THIS IS INCORRECT
    :return: 1st bending frequency
    """
    axial = 100
    circumferential = 100
    # Case 8: Beam Bending, i=1, Modes of SS, Natural Frequency of Plates and shells p.239
   
    # lambda_ij = (j * pi * r_o / L) ** 2 * sqrt((1 - v ** 2) / 2)

    

    # f_ij_pressurized = sqrt( f_ij**2 + 1 / ( 4 * pi**2 ) * ( Nx / (rho *t) * (j * pi / L)**2 + Nt / (rho *t) * (j *
    # pi / L)**2 ) )

    # f_ij_pressurized = sqrt(f_ij ** 2 + 1 / (4 * pi ** 2) * (Nt / (rho * t) * (j * pi / L) ** 2))

    #Case 9: Simply supported cylindrical shell
    freq_lst =[]
    for i in range (2, circumferential):
        for j in range (1, axial):
            
            lambda_ij = (( (1-v**2) *( j * pi * r_o / L)**4 + t**2 / (12 * r_o**2) * (i**2 +(j * pi * r_o / L)**2)**4) )**(1/2) / (i**2 + (j* pi * r_o / L)**2 )
    
            f_ij = lambda_ij / (2 * pi * r_o) * sqrt(E / (rho * (1 - v ** 2)))
            freq_lst.append(f_ij)

    return min(freq_lst)


# test  =  find_1st_bending_freq_cylindrical_shell(r_o = 0.175, L = 0.2, t =0.004, v =0.33, E= 69E9, rho = 2850)
# print(test)