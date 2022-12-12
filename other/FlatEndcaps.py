# from databases import materials
from math import pi 
'''
Source: https://www.engineersedge.com/material_science/circular_plate_edges_clamped_13639.htm
'''
def flat(p , r, t , D):
    stress_max = 3 * p * r**2 / ( 4 * t**2 )
    y_max = p * r ** 4 / ( 64 * D )
    return stress_max, y_max

def flexural_rigidity(E, v , t):
    return E * t**3 / ( 12 * ( 1 - v**2 ) )

if __name__ in "__main__":
    yield_ = 80E6
    E = 69E9
    v = 0.33
    stress_max = 250E6 
    y_max = 0.1
    t = 1E-3
    p = 45E6
    rho = 2700
    r = 0.366/2
    while stress_max >= yield_ / 1.567 or y_max >= 0.005:
        D = flexural_rigidity(E, v, t)
        stress_max, y_max = flat(p, r, t, D)
        t+=0.001

    Mass = t * pi * r**2 * rho
    print(f"Thicknes: {t}, Mass: {Mass}, Mass deflection: {y_max}")
