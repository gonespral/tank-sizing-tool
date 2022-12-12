from math import ln
def isPV(P_1, V, gamma_, P_2):
    """
    Finds if a vessel is a pressure vessel from ECSS-E-ST-32-02C
    :param P_1: Internal pressure in Pa
    :param P_2: external pressure in Pa
    :param V: volume in m^3
    :param gamma_: specific heat ratio
    :return: boolean
    """
    E = P_1 * V / ( gamma_ - 1 ) * ( 1 - ( P_2 / P_1 ) ** ( gamma_ -1 / gamma_ ) )  

    if P_1 >= 0.69E6 or E>= 19310:
        return True
    else:
        return False
    
# if __name__ in "__main__":
#     # P_1 = 2.5E5
#     # P_2 = 1E5
#     # V = 65/1000 
#     # gamma_ = 1.395 
#     P_1 = 2.0E5 
#     P_2 = 1E5 
#     V = 61.5/1000 
#     gamma_ = 1.18 
#     P_1 = 3.5E5
#     P_2 = 1E5
#     V = 50/1000 
#     gamma_ = 1.67
#     print(isPV(P_1, V, gamma_, P_2))

if __name__ in "__main__":
    UTS = 10.66 -19.42*ln(1-60/109)
    print(f'Magnitude: {UTS}')