# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 10:28:26 2023

@author: Luca
"""

# Example thermal model of a tank from LOx to skirt interface

from CoolProp.CoolProp import PropsSI
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

# gravitational acceleration
g = 9.81 # [m/s^2]
# tank pressure
P_tank = 45e5 # [Pa]
# height of fluid column
h_fluid_column = 0.355
# number of fluid layers
n_layer_fluid = 10
# number of skirt layers
n_layer_skirt = 10
# skirt thickness
t_skirt = 5e-3
# tank diameter
d_tank = 0.355
# height of the skirt
h_skirt = 0.3
# number of dome layers
n_layer_dome = 5
# thickness of dome
t_dome = 6e-3
# ambient temperature
T_ambient = 273

# alu properties, I use constant values for now
rho_alu = 2700 # assumed constant, 2700 at room temperature and 2731 at 90 K
cp_alu = 422 # [J/(kg*K)], depends on temperature, like 422 for 90K and about 820 for room temperature
k_alu = 30 # [W/(mK)] depends on temperature, 69 at 293K and 30 at 90K


# ---------- first calculate constant parameters
# fluid column
T_inf_lox = 90  # [K] # temperature of LOx bulk
# length of conductive path in LOx
L_lox = np.full(n_layer_fluid, h_fluid_column/n_layer_fluid)
# cross sectional area of fluid column
A_lox = np.full(n_layer_fluid, np.pi*(d_tank/2)**2)

# dome
# area of dome surface
A_dome = np.full(n_layer_dome, np.pi*(d_tank/2)**2 )# circular plate, todo: make it represent the dome shape (ellipsoid or cassinian or smth)
# length of conductive path through dome
L_dome = np.full(n_layer_dome, t_dome/n_layer_dome)
# mass of dome nodes
m_dome = A_dome*L_dome*rho_alu

# skirt
# area of skirt ring
A_skirt = np.full(n_layer_skirt, np.pi*(((d_tank+2*t_skirt)/2)**2-(d_tank/2)**2))
# length of conductive path thorugh skirt
L_skirt = np.full(n_layer_skirt, h_skirt/n_layer_skirt)
# mass of skirt nodes
m_skirt = A_skirt*L_skirt*rho_alu

# set up simulation time
t_start = 0
t_stop = 3600
# number of discrete time steps
delta_t = 0.01
n_t = int((t_stop-t_start)/delta_t)
t_sim = np.arange(t_start,t_stop,delta_t)

# set up temperature vector with initial temperatures
n_layer = n_layer_fluid+n_layer_dome+n_layer_skirt
T = np.zeros((n_layer,n_t))
T_lox = np.full(n_layer_fluid,T_inf_lox)
T_dome = np.full(n_layer_dome,T_ambient)
T_skirt = np.full(n_layer_skirt, T_ambient)
T_inital = np.concatenate((T_lox, T_dome, T_skirt))
T[:,0] = T_inital

# simulation loop
for i in tqdm(range(len(t_sim)-1)):
    # get fluid properties for current temperature
    # conductivity
    k_lox = PropsSI('L', 'T', T[0:n_layer_fluid,i], 'P', P_tank, 'O2')
    # dynamic viscosity
    mu_lox = PropsSI('V', 'T', T[0:n_layer_fluid,i], 'P', P_tank, 'O2')
    # density 
    rho_lox = PropsSI('D', 'T', T[0:n_layer_fluid,i], 'P', P_tank, 'O2')
    # specific heat
    cp_lox = PropsSI('C', 'T', T[0:n_layer_fluid,i], 'P', P_tank, 'O2')
    # Prandtl number
    Pr = PropsSI('Prandtl', 'T', T[0:n_layer_fluid,i], 'P', P_tank, 'O2')
    # kinematic viscosity
    nu_lox = mu_lox[n_layer_fluid-1]/rho_lox[n_layer_fluid-1]
    # thermal expansion coefficient for LOx
    beta_lox = 9.53e-3
    # Grashof number
    Gr = g*beta_lox*(T[n_layer_fluid,i]-T[n_layer_fluid-1,i])*d_tank**3/nu_lox**2
    # Rayleigh number
    Ra = Gr*Pr[n_layer_fluid-1]
    # Nusselt number
    Nu = 0.15*Ra**(1/3)
    
    # mass of fluid layers
    m_lox = A_lox*L_lox*rho_lox
    # heat capacity of fluid layers
    c_lox = cp_lox*m_lox
    
    # heat capacity of dome nodes
    c_dome = m_dome*cp_alu
    # heat capacity of skirt nodes
    c_skirt = m_skirt*cp_alu
    
    # heat capacity vector
    c = np.concatenate((c_lox,c_dome,c_skirt))
    
    # conductive couplings
    GL_lox = (k_lox*A_lox/L_lox)[:-1]
    GL_dome = k_alu*A_dome/L_dome
    GL_skirt = (k_alu*A_skirt/L_skirt)[:-1]
    
    # convection from fluid to dome
    h_conv = Nu*k_lox[n_layer_fluid-1]/d_tank
    # convective coupling
    GC = np.full(1,h_conv*A_lox[n_layer_fluid-1])
    
    # coupling matrix
    # this is very badly set up with dumb variable names
    # but it works for now
    couplings = np.concatenate((GL_lox,GC,GL_dome,GL_skirt))
    A = np.eye(n_layer,k=1)*np.concatenate(([0],couplings))
    B = np.eye(n_layer,k=-1)*np.concatenate((couplings,[0]))
    d = np.concatenate(([0],couplings,[0]))
    e = np.zeros(len(d)-1)
    for k in range(len(d)-1):
        e[k] = -d[k+1]-d[k]
    C = np.eye(n_layer)*e
    # D is the complete coupling matrix
    D = A+B+C
    
    # calculate temperature change
    delta_T = delta_t/c*np.matmul(D,T[:,i])
    # calculate temperature of next time step
    T[:,i+1] = T[:,i]+delta_T
    
# plotting the node at the end of the skirt
plt.plot(t_sim,T[-1])
plt.grid()
plt.xlabel("Time after LOx filling [s]")
plt.ylabel("Temperature at the interface [K]")
plt.show()
    
    
    
    
    
    