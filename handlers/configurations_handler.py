"""
Functions for manipulating and generating configurations sets.
"""

from sections import cylinder
from sections import dome
from sections import skirt
from sections import configuration

from databases import materials
from databases import fluids

from formulas.nasa_sp_8007 import elastic_constants
from formulas.bsc_ae import geometry

import numpy as np
import pickle as pkl
import datetime
from tqdm import trange
import itertools


def generate(configurations: list, config) -> list:
    """
    Generates all configurations using tqdm progress bar.
    :param configurations: Configurations list
    :param config: Config file
    :return: Populated configurations list
    """
    # Get total num of configurations to be generated, from a set of lists, using itertools
    total_num_of_configurations = len(list(itertools.product(config.dome_types, config.dome_inner_radii,
                                                             config.dome_thicknesses, config.dome_heights,
                                                             config.dome_materials, config.dome_cassinian_ns,
                                                             config.cylinder_inner_radii, config.cylinder_thicknesses,
                                                             config.cylinder_heights, config.cylinder_materials,
                                                             config.skirt_extensions, config.skirt_materials)))

    print(f"Generating {total_num_of_configurations} configurations...")

    with trange(total_num_of_configurations) as pbar:
        for dome_inner_radius in config.dome_inner_radii:
            for dome_thickness in config.dome_thicknesses:
                for dome_height in config.dome_heights:
                    for dome_material in config.dome_materials:
                        for dome_type in config.dome_types:
                            for dome_n in config.dome_cassinian_ns:
                                for cylinder_inner_radius in config.cylinder_inner_radii:
                                    for cylinder_thickness in config.cylinder_thicknesses:
                                        for cylinder_height in config.cylinder_heights:
                                            for cylinder_material in config.cylinder_materials:
                                                for skirt_extension in config.skirt_extensions:
                                                    for skirt_material in config.skirt_materials:
                                                        try:
                                                            pbar.update(1)

                                                            # Convert inner radius to outer radius
                                                            dome_outer_radius = dome_inner_radius + dome_thickness
                                                            cylinder_outer_radius = cylinder_inner_radius + cylinder_thickness

                                                            # Check that dome and cylinder are overlapping
                                                            if dome_outer_radius - dome_thickness > cylinder_outer_radius \
                                                                or dome_outer_radius < cylinder_outer_radius - cylinder_thickness:
                                                                raise ValueError("Dome and cylinder are not "
                                                                                 "overlapping.")

                                                            if dome_thickness >= cylinder_thickness:
                                                                raise ValueError("Dome is thicker than the cylinder.")

                                                            # Derived parameters
                                                            skirt_thickness = cylinder_thickness - dome_thickness
                                                            skirt_outer_radius = dome_outer_radius + skirt_thickness

                                                            # Calculate stiffener elastic constants for
                                                            # ring-stringer method.
                                                            stiffener_stringer_spacing = (2 * np.pi * (skirt_outer_radius - skirt_thickness / 2)) / config.skirt_stringer_num
                                                            stiffener_ring_spacing = (skirt_extension + dome_height) / config.skirt_ring_num
                                                            stiffener_stringer_cross_sectional_area = config.skirt_stringer_height * config.skirt_stringer_thickness
                                                            stiffener_ring_cross_sectional_area = config.skirt_ring_height * config.skirt_ring_thickness
                                                            stiffener_ring_I = geometry.rectangle_I(config.skirt_ring_thickness, config.skirt_ring_height)
                                                            stiffener_stringer_I = geometry.rectangle_I(config.skirt_stringer_thickness, config.skirt_stringer_height)
                                                            stiffener_ring_J = geometry.rectangle_J(config.skirt_ring_thickness, config.skirt_ring_height)
                                                            stiffener_stringer_J = geometry.rectangle_J(config.skirt_stringer_thickness, config.skirt_stringer_height)
                                                            stiffener_ring_z = -(config.skirt_ring_height / 2) if config.skirt_ring_inside_surface else (config.skirt_ring_height / 2)
                                                            stiffener_stringer_z = -(config.skirt_stringer_height / 2) if config.skirt_stringer_inside_surface else (config.skirt_stringer_height / 2)
                                                            skirt_material_E = materials.materials[skirt_material]["youngs_modulus"]
                                                            skirt_material_G = materials.materials[skirt_material]["shear_modulus"]
                                                            skirt_material_v = materials.materials[skirt_material]["poisson_ratio"]
                                                            ec = elastic_constants.isotropic_cylinders_with_rings_and_stringers(
                                                                E=skirt_material_E,
                                                                E_s=skirt_material_E,
                                                                E_r=skirt_material_E,
                                                                G_s=skirt_material_G,
                                                                G_r=skirt_material_G,
                                                                A_s=stiffener_stringer_cross_sectional_area,
                                                                A_r=stiffener_ring_cross_sectional_area,
                                                                I_s=stiffener_stringer_I,
                                                                I_r=stiffener_ring_I,
                                                                J_s=stiffener_stringer_J,
                                                                J_r=stiffener_ring_J,
                                                                b_s=stiffener_stringer_spacing,
                                                                b_r=stiffener_ring_spacing,
                                                                z_s=stiffener_stringer_z,
                                                                z_r=stiffener_ring_z,
                                                                v=skirt_material_v,
                                                                t=skirt_thickness)

                                                            # Create a forward dome object
                                                            f_dome_object = dome.Dome(type_=dome_type,
                                                                                      outer_radius=dome_outer_radius,
                                                                                      thickness=dome_thickness,
                                                                                      height=dome_height if dome_type == "semi-ellipsoidal" else None,
                                                                                      material=materials.materials[dome_material],
                                                                                      n=dome_n if dome_type == "cassinian" else None)

                                                            # Create an aft dome object
                                                            a_dome_object = dome.Dome(type_=dome_type,
                                                                                      outer_radius=dome_outer_radius,
                                                                                      thickness=dome_thickness,
                                                                                      height=dome_height if dome_type == "semi-ellipsoidal" else None,
                                                                                      material=materials.materials[dome_material],
                                                                                      n=dome_n if dome_type == "cassinian" else None)

                                                            # Create a forward skirt object
                                                            f_skirt_object = skirt.Skirt(outer_radius=skirt_outer_radius,
                                                                                         thickness=skirt_thickness,
                                                                                         extension=skirt_extension,
                                                                                         material=materials.materials[skirt_material],
                                                                                         stiffening_elastic_constants=ec,
                                                                                         stiffening_volume=config.stiffening_volume)

                                                            # Create an aft skirt object
                                                            a_skirt_object = skirt.Skirt(outer_radius=skirt_outer_radius,
                                                                                         thickness=skirt_thickness,
                                                                                         extension=skirt_extension,
                                                                                         material=materials.materials[skirt_material],
                                                                                         stiffening_elastic_constants=ec,
                                                                                         stiffening_volume=config.stiffening_volume)

                                                            # Create a cylinder object
                                                            cylinder_object = cylinder.Cylinder(outer_radius=cylinder_outer_radius,
                                                                                                thickness=cylinder_thickness,
                                                                                                height=cylinder_height,
                                                                                                material=materials.materials[cylinder_material])

                                                            # Create a configuration object
                                                            configuration_object = configuration.Configuration(forward_dome=f_dome_object,
                                                                                                               aft_dome=a_dome_object,
                                                                                                               forward_skirt=f_skirt_object,
                                                                                                               aft_skirt=a_skirt_object,
                                                                                                               cylinder=cylinder_object,
                                                                                                               fluid=fluids.fluids[config.fluid],
                                                                                                               fluid_volume=config.fluid_volume,
                                                                                                               internal_pressure=config.internal_pressure)

                                                            configurations.append(configuration_object)

                                                        except Exception as e:
                                                            # print(f"Exception: {e}. Skipping configuration...")
                                                            continue

    if len(configurations) == 0:
        raise Exception("No configurations generated. Check your parameters.")
    elif len(configurations) != total_num_of_configurations:
        print(f"Generated {len(configurations)} configurations out of {total_num_of_configurations}.")
    else:
        print(f"Generated {len(configurations)} configurations")
    return configurations


def sort_by_mass(configurations: list) -> list:
    """
    Sorts all configurations by the total weight of the tank.
    :param configurations: Configurations list object
    :return: sorted_configurations list
    """
    print("Sorting configurations by mass")
    sorted_configurations = sorted(configurations, key=lambda configuration: configuration.mass)
    return sorted_configurations

def sort_by_height(configurations: list) -> list:
    """
    Sorts all configurations by the height of the tank.
    :param configurations: Configurations list object
    :return: sorted_configurations list
    """
    print("Sorting configurations by height")
    sorted_configurations = sorted(configurations, key=lambda configuration: configuration.height)
    return sorted_configurations


def load_from_file(path: str) -> tuple:
    """
    Loads configurations from a pickle file, containing (configurations, config) tuple.
    :param path: Path to the file
    :return: (configurations list, config object)
    """
    print("Loading configurations from file")
    with open(path, "rb") as f:
        configurations = pkl.load(f)
    return configurations


def save_to_file(configurations: tuple, path: str = "saves/"):
    """
    Saves configurations to a pickle file, containing (configurations, config) tuple.
    :param configurations: Configurations tuple object (passed_configurations, failed_configurations)
    :param path: Path to the file
    """
    filename = f"{path}configurations_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.pkl"
    print(f"Saving configurations to file {filename}. This may take a while.")
    with open(filename, "wb") as f:
        pkl.dump(configurations, f, protocol=pkl.HIGHEST_PROTOCOL)
    print("Done")
