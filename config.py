# Configuration file. Change these parameters to generate different configurations.

# All units are in SI, unless otherwise specified.

import numpy as np
import valispace
from databases import fluids

# Initialize Valispace
valispace = valispace.API(url='https://projectstratos.valispace.com/', username='', password='')

# ############################# General parameters ##########################
# These parameters can be used in all sections

thicknesses = np.arange(0.003, 0.006, 0.0001)

inner_radii = [0.280 / 2]

# ############################ Fluid parameters ############################
fluid = "lox"  # See fluids database

# Fuel: 13466
# Oxidizer: 13469
fluid_mass = valispace.get_vali(13466)["value"] if fluid == "ethanol" else valispace.get_vali(13469)["value"]
fluid_volume = fluid_mass / fluids.fluids[fluid]["density"]

# Fuel: 13470
# Oxidizer: 13471
internal_pressure = 100000 * valispace.get_vali(13470)["value"] if fluid == "ethanol" else 100000 * valispace.get_vali(13471)["value"]

# Ullage_min_fuel = 13465
# Ullage_min_oxidizer = 14257
# Ullage_max : ullage_min + 0.05
ullage_range = [valispace.get_vali(13465)["value"] if fluid == "ethanol" else valispace.get_vali(14257)["value"], 0]
if fluid == "lox": 
    ullage_range[0] += valispace.get_vali(14254)["value"]
ullage_range[1] = ullage_range[0] + 0.05
# 
# ############################ Dome iterables ##############################
dome_inner_radii = inner_radii

dome_thicknesses = thicknesses

dome_materials = ["6082-T6 aluminium"]  # See materials database

dome_types = ["semi-ellipsoidal"]  # See sections/dome.py

# -- Cassinian dome parameters --
dome_cassinian_ns = np.arange(1.5, 2.0, 0.1)

# -- Semi-ellipsoidal dome parameters --
dome_heights = np.arange(0.05, max(dome_inner_radii) + max(dome_thicknesses), 0.01)

# ############################# Cylinder iterables #########################
cylinder_inner_radii = inner_radii

cylinder_thicknesses = thicknesses

cylinder_heights = np.arange(0.1, 0.8, 0.005)

cylinder_materials = ["6082-T6 aluminium"]  # See materials database

# ############################# Skirt iterables #############################
# Skirt thickness is determined by the thickness of the dome/cylinder -> t_skirt = t_cylinder - t_dome
# Skirt radius is determined by the outer radius of the dome/cylinder -> r_skirt = r_dome + t_dome

skirt_extensions = [0]  # Extension of skirt from top of dome (extra space for feed)

skirt_materials = ["6082-T6 aluminium"]  # See materials database. Material is the same for stiffeners and skirt.

# -- Ring-stringer parameters --
skirt_ring_thickness = 0.000001
skirt_ring_height = 0.000001
skirt_ring_num = 10  # Circumferentially - Make sure this not too low, assumptions will not work.

skirt_ring_inside_surface = True  # True for inside, False for outside

skirt_stringer_thickness = 0.002
skirt_stringer_height = 0.002
skirt_stringer_num = 8  # Axially

stringer_volume = skirt_stringer_thickness * skirt_stringer_height * skirt_stringer_num
ring_volume = skirt_ring_thickness * skirt_ring_height * skirt_ring_num
stiffening_volume = ring_volume + stringer_volume

skirt_stringer_inside_surface = True  # True for inside, False for outside