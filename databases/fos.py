"""
FoS for pressurized metallic structures according to ECSS-E-ST-32-02C and ECSS-E-ST-32-10C

DO NOT CHANGE WITHOUT PERMISSION.
"""

'''
Due to the fact that the leak is hazardous for LOx, the vessel is desinged with non LBB behaviour (p.36)
The safe life of the vessel has to be demonsrated (by analysis or test or both): pre-flawed mettalic items,
leak tighness and no rupture after 4 x service life (Actually it should be 5)
'''

# Minimum values of factors of safety for internal pressure
'''
FOSY_p: yield design factor of safety for internal pressure
FOSU_p: ultimate design factor of safety for internal pressure
PF: proof factor
'''
FOSY_p = 1.3  # FAA , last updated: 14.03.2023
FOSU_p = 1.5 # FAA , last updated: 14.03.2023
ProofFactor = 1.25  #FAA, last updated: 14.03.2023
BurstFactor = 1.5 #FAA, last updated: 14.03.2023 
# -----------------------------------------------------------------------------------------------------
# Factors of safety for external loads
'''
LL: limit load
DLL: design limit load = DF * LL
DUL: design ultimate load = DLL * coresspodning FOSU * corresponding K_LD > FOSU_p
DYL: design yield load = DLL * coresspodning FOSY * corresponding K_LD > FOSY_p
'''

# FoS to adjust for project stage and confidance in the model
'''
K_M: model factor, justified based on previous missions
K_P: project factor, justfied based on previoius missions
DF: design factor
'''
K_M_Thrust = 1.1 #14.03.2023
K_M_Loads = 1.2   # Estimated guess for the preliminary tool, last updated: 14.03.2023
K_M_Pressure = 1.0  # Final level of subsystem, last updated: 14.03.2023
K_P = 1.05
Coeff_A_Thrust = K_M_Thrust * K_P
Coeff_A_Loads = K_M_Loads * K_P
Coeff_A_Pressure = K_M_Pressure * K_P
# -----------------------------------------------------------------------------------------------------
# Constant safety factors for tanks subsytem, metallic strucutres
'''
Veryfication of the subsystem is by analysis only

FOSU: ultimate design factor of safety
FOSY: yield design factor of safety
K_MP: margin policy factor, valid solely for the tank subsytems 
'''
K_MP = 1.15  # Ariane 5, 17.02.2023
FOSY = 1.3 # FAA for test, last updated: 06.04.2023, was 1.25 Jan wants it 
FOSU = 1.5 # FAA for test, last updated: 14.03.2023
Coeff_B = FOSY * K_MP 

Coeff_C = FOSU * K_MP

# -----------------------------------------------------------------------------------------------------
# Safety factors for joints, inserts and connections ANALYSIS ONLY
'''
FOSY_cij: yield design factor of safety for joints, inserts and connections
FOSU_cij: ultimate design factor of safety for joints, inserts and connections
'''
FOSY_cij = 2.0  # ECSS for analysis, last updated: 09.01.2023
FOSU_cij = 2.0  # ECSS for analysis, last updated: 09.01.2023
# -----------------------------------------------------------------------------------------------------
# Safety factors for buckling ANALYSIS ONLY
'''
FOSY_b_local: yield design factor of safety for local buckling
FOSU_b_local: ultimate design factor of safety for local buckling
FOSY_b_global: yield design factor of safety for global buckling
FOSU_b_global: ultimate design factor of safety for global buckling
'''
FOSY_b_local = 2.0  # ECSS for test, last updated: 21.02.2023
FOSU_b_local = 1.25  # ECSS for analysis, last updated: 21.02.2023
FOSY_b_global = 2.0 # ECSS for analysis, last updated: 21.02.2023
FOSU_b_global = 1.25  # ECSS for analysis, last updated: 21.02.2023
# -----------------------------------------------------------------------------------------------------
# Safety factor for local components NOT FINISHED
'''
K_LD: local design factor, this factor accounts for specific uncertainties linked to the analysis difficulties or to the lack
of reliable dimensioning methodology or criteria where significant stress gradients occur (e.g. geometric singularities, fitting, welding,
riveting, bonding, holes, inserts and, for composite, lay-up drop out, sandwich core thickness change, variation of ply consolidation
as a result of drape over corners)
'''
# Each time you add new component a new local factor has to be added
K_LD = 1.2  # Guestimated, last updated: 09.01.2023
# -----------------------------------------------------------------------------------------------------
