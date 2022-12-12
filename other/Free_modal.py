import numpy as np 

def free(i, h, v, R, L, rho):
    num = i**2 * ( i**2 - 1 )**2 * h**2 * (1 + 24 * ( 1- v ) * R**2 / ( i**2 * L**2))
    den = (i**2 + 1) * (12 * R**2) * (1+ 12 * R**2) / (i**2 * (i**2 +1) * L**2)
    
    lambda_ij = (num / den)**(1/2)

    f_ij = lambda_ij / (2 * np.pi * R) * np.sqrt(E / (rho * (1 - v ** 2)))
            # freq_lst.append(f_ij)
    return f_ij

if __name__ in "__main__":
    i = 2 
    R = 0.1435 #m
    h = 0.005 #m
    v = 0.33
    L = 6.5
    E = 70E9
    rho = 2700

    freq  =free(i,h,v,R,L,rho)

    print(f"Freq: {round(freq,2)}Hz")
