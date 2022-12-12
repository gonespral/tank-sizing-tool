"""
Validation handler.
"""
import numpy as np
from tqdm import trange

from formulas.nasa_sp_8007 import buckling
from formulas.nasa_sp_8007 import bending
from formulas.nasa_sp_8007 import hoop_stress
from databases.fos import *
from formulas.flugge import yield_criterion

from formulas.affdl import load_interaction as li

from formulas.donnel import vibrations as vib

def max_mass(configurations: list, lim: int) -> tuple:
    """
    Sets a maximum mass for the tank configuration.
    :param configurations: Configurations list object
    :param lim: Maximum mass of the tank
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:
            if configuration.mass > lim:
                configuration.validator_parameters["max_mass"] = {"parameter": "configuration mass",
                                                                  "value": configuration.mass, "status": "failed"}
                failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["max_mass"] = {"parameter": "configuration mass",
                                                                  "value": configuration.mass, "status": "passed"}
                passed_configurations.append(configuration)
            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to mass being greater than {lim} kg."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


# ############################# CYLINDER ########################################

def cylinder_hoop_collapse(configurations: list) -> tuple:
    """
    Validates all configurations for yielding due to the hoop stress
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:
            stress = hoop_stress.hoop_stress(p=configuration.internal_pressure,
                                             r=configuration.cylinder.outer_radius,
                                             t=configuration.cylinder.thickness)

            if stress >= configuration.cylinder.material['yield_stress']:
                configuration.validator_parameters["cylinder_hoop_collapse"] = {"parameter": "hoop stress",
                                                                                "value": stress,
                                                                                "status": "failed"}
                failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["cylinder_hoop_collapse"] = {"parameter": "hoop stress",
                                                                                "value": stress,
                                                                                "status": "passed"}
                passed_configurations.append(configuration)
            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to hoop stress collapse."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


def cylinder_axial_buckling(configurations: list, loads: dict) -> tuple:
    """
    Validates all configurations for axial buckling of the cylinder material.
    :param loads: Loading dictionary
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:

            critical_buckling_load = buckling.critical_cylinder_buckling(p=configuration.internal_pressure,
                                                                         r=configuration.cylinder.outer_radius,
                                                                         t=configuration.cylinder.thickness,
                                                                         l=configuration.cylinder.height,
                                                                         E=configuration.cylinder.material[
                                                                             'youngs_modulus'],
                                                                         v=configuration.cylinder.material[
                                                                             'poisson_ratio'])
            if abs(loads["axial"]) >= critical_buckling_load:
                configuration.validator_parameters["cylinder_axial_buckling"] = {
                    "parameter": "critical axial buckling load", "value": critical_buckling_load, "status": "failed"}
                failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["cylinder_axial_buckling"] = {
                    "parameter": "critical axial buckling load", "value": critical_buckling_load, "status": "passed"}
                passed_configurations.append(configuration)
            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to axial buckling."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


def cylinder_bending_buckling(configurations: list, loads: dict) -> tuple:
    """
    Validates all configurations for cylinder bending.
    :param loads: Loading dictionary
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:

            critical_buckling_load = bending.critical_cylinder_bending(p=configuration.internal_pressure,
                                                                       r=configuration.cylinder.outer_radius,
                                                                       t=configuration.cylinder.thickness,
                                                                       E=configuration.cylinder.material[
                                                                           'youngs_modulus'],
                                                                       v=configuration.cylinder.material[
                                                                           'poisson_ratio'])
            if abs(loads["moment"]) >= critical_buckling_load:
                configuration.validator_parameters["cylinder_bending_buckling"] = {
                    "parameter": "critical bending buckling load", "value": critical_buckling_load, "status": "failed"}
                failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["cylinder_bending_buckling"] = {
                    "parameter": "critical bending buckling load", "value": critical_buckling_load, "status": "passed"}
                passed_configurations.append(configuration)
            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to bending buckling."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


def cylinder_load_interaction(configurations: list, loads: dict) -> tuple:
    """
    Validates all configurations for cylinder loads interaction.
    :param loads: Loading dictionary
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    # Assign loads. If zero, set them to None, so that the correct interaction equation can then be used.
    applied_compressive_load = loads["axial"] if loads["axial"] < 0 else None
    applied_bending_moment = loads["moment"] if loads["moment"] != 0 else None

    if applied_compressive_load and applied_bending_moment:

        with trange(len(configurations)) as pbar:
            for configuration in configurations:

                critical_compressive_load = buckling.critical_cylinder_buckling(p=configuration.internal_pressure,
                                                                                r=configuration.cylinder.outer_radius,
                                                                                t=configuration.cylinder.thickness,
                                                                                l=configuration.cylinder.height,
                                                                                E=configuration.cylinder.material[
                                                                                    "youngs_modulus"],
                                                                                v=configuration.cylinder.material[
                                                                                    "poisson_ratio"])

                critical_bending_moment = bending.critical_cylinder_bending(p=configuration.internal_pressure,
                                                                            r=configuration.cylinder.outer_radius,
                                                                            t=configuration.cylinder.thickness,
                                                                            E=configuration.cylinder.material[
                                                                                "youngs_modulus"],
                                                                            v=configuration.cylinder.material[
                                                                                "poisson_ratio"])

                if applied_compressive_load and applied_bending_moment:
                    # Axial compression and bending moment interaction
                    R_c = applied_compressive_load / critical_compressive_load
                    R_b = applied_bending_moment / critical_bending_moment

                    if li.axial_comp__pure_bending(R_c, R_b):
                        configuration.validator_parameters["cylinder_load_interaction"] = {
                            "parameter": "axial compression and bending moment interaction", "value": "failed",
                            "status": "failed"}
                        failed_configurations.append(configuration)
                    else:
                        configuration.validator_parameters["cylinder_load_interaction"] = {
                            "parameter": "axial compression and bending moment interaction", "value": "passed",
                            "status": "passed"}
                        passed_configurations.append(configuration)

                pbar.update(1)

    else:
        print("No axial compression and bending moment interaction to validate.")

    print(f"Failed {len(failed_configurations)} configurations due to loads interaction."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations

def cylinder_stiffness(configurations: list) -> tuple:
    """
    Validates all configurations for cylinder 1st bending frequency.
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:                                                         
        
            freq = vib.find_1st_bending_freq_cylindrical_shell( r_o=configuration.cylinder.outer_radius,
                                                                L=configuration.cylinder.height,
                                                                t=configuration.cylinder.thickness,
                                                                v= configuration.cylinder.material["poisson_ratio"],
                                                                E = configuration.cylinder.material["youngs_modulus"],
                                                                rho = configuration.cylinder.material["density"])
            if freq < 100 * Coeff_A_Loads * Coeff_C: 
                failed_configurations.append(configuration)
            else:
                passed_configurations.append(configuration)

            pbar.update(1)


    print(f"Failed {len(failed_configurations)} configurations due to 1st bending mode."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations

# ########################## SKIRT #######################################

def skirt_axial_buckling(configurations: list, loads: dict) -> tuple:
    """
    Validates all configurations for axial buckling of the skirts.
    :param loads: Loading dictionary
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:
            # Forward skirt
            ec = configuration.forward_skirt.stiffening_elastic_constants
            L = configuration.forward_skirt.height
            r = (configuration.forward_skirt.outer_radius + configuration.forward_skirt.inner_radius) / 2
            t = configuration.forward_skirt.thickness
            v = configuration.forward_skirt.material['poisson_ratio']
            critical_axial_load = buckling.critical_stiffened_cylinder_buckling(*ec, L, r, t, v)

            if abs(loads["axial"]) >= critical_axial_load:
                configuration.validator_parameters["forward_skirt_axial_buckling"] = {
                    "parameter": "critical axial buckling load", "value": critical_axial_load, "status": "failed"}
                failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["forward_skirt_axial_buckling"] = {
                    "parameter": "critical axial buckling load", "value": critical_axial_load, "status": "passed"}
                passed_configurations.append(configuration)

            # Aft skirt
            ec = configuration.aft_skirt.stiffening_elastic_constants
            L = configuration.aft_skirt.height
            r = (configuration.aft_skirt.outer_radius + configuration.aft_skirt.inner_radius) / 2
            t = configuration.aft_skirt.thickness
            v = configuration.aft_skirt.material['poisson_ratio']
            critical_axial_load = buckling.critical_stiffened_cylinder_buckling(*ec, L, r, t, v)

            if abs(loads["axial"]) >= critical_axial_load:
                configuration.validator_parameters["aft_skirt_axial_buckling"] = {
                    "parameter": "critical axial buckling load", "value": critical_axial_load, "status": "failed"}
                if configuration not in failed_configurations: failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["aft_skirt_axial_buckling"] = {
                    "parameter": "critical axial buckling load", "value": critical_axial_load, "status": "passed"}
                if configuration not in passed_configurations: passed_configurations.append(configuration)

            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to skirt axial buckling."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


def skirt_bending_buckling(configurations: list, loads: dict) -> tuple:
    """
    Validates all configurations for bending buckling of the skirts.
    :param loads: Loading dictionary
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:
            # Forward skirt
            ec = configuration.forward_skirt.stiffening_elastic_constants
            L = configuration.forward_skirt.height
            r = (configuration.forward_skirt.outer_radius + configuration.forward_skirt.inner_radius) / 2
            t = configuration.forward_skirt.thickness
            v = configuration.forward_skirt.material['poisson_ratio']
            critical_bending_load = bending.critical_stiffened_cylinder_bending(*ec, L, r, t, v)

            if abs(loads["moment"]) >= critical_bending_load:
                configuration.validator_parameters["forward_skirt_bending_buckling"] = {
                    "parameter": "critical bending buckling load", "value": critical_bending_load, "status": "failed"}
                failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["forward_skirt_bending_buckling"] = {
                    "parameter": "critical bending buckling load", "value": critical_bending_load, "status": "passed"}
                passed_configurations.append(configuration)

            # Aft skirt
            ec = configuration.aft_skirt.stiffening_elastic_constants
            L = configuration.aft_skirt.height
            r = (configuration.aft_skirt.outer_radius + configuration.aft_skirt.inner_radius) / 2
            t = configuration.aft_skirt.thickness
            v = configuration.aft_skirt.material['poisson_ratio']
            critical_bending_load = bending.critical_stiffened_cylinder_bending(*ec, L, r, t, v)

            if abs(loads["moment"]) >= critical_bending_load:
                configuration.validator_parameters["aft_skirt_bending_buckling"] = {
                    "parameter": "critical bending buckling load", "value": critical_bending_load, "status": "failed"}
                if configuration not in failed_configurations: failed_configurations.append(configuration)
            else:
                configuration.validator_parameters["aft_skirt_bending_buckling"] = {
                    "parameter": "critical bending buckling load", "value": critical_bending_load, "status": "passed"}
                if configuration not in passed_configurations: passed_configurations.append(configuration)

            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to skirt bending buckling."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


def skirt_load_interaction(configurations: list, loads: dict) -> tuple:
    """
    Validates all configurations for skirt loads interaction.
    :param loads: Loading dictionary
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    # Assign loads. If zero, set them to None, so that the correct interaction equation can then be used.
    applied_compressive_load = loads["axial"] if loads["axial"] < 0 else None
    applied_bending_moment = loads["moment"] if loads["moment"] != 0 else None

    if applied_compressive_load and applied_bending_moment:

        with trange(len(configurations)) as pbar:
            for configuration in configurations:

                for skirt in [configuration.forward_skirt, configuration.aft_skirt]:

                    critical_compressive_load = buckling.critical_cylinder_buckling(p=configuration.internal_pressure,
                                                                                    r=skirt.outer_radius,
                                                                                    t=skirt.thickness,
                                                                                    l=skirt.height,
                                                                                    E=skirt.material["youngs_modulus"],
                                                                                    v=skirt.material["poisson_ratio"])

                    critical_bending_moment = bending.critical_cylinder_bending(p=configuration.internal_pressure,
                                                                                r=skirt.outer_radius,
                                                                                t=skirt.thickness,
                                                                                E=skirt.material["youngs_modulus"],
                                                                                v=skirt.material["poisson_ratio"])

                    # Axial compression and bending moment interaction
                    R_c = applied_compressive_load / critical_compressive_load
                    R_b = applied_bending_moment / critical_bending_moment

                    if li.axial_comp__pure_bending(R_c, R_b):
                        configuration.validator_parameters["skirt_load_interaction"] = {
                            "parameter": "axial compression and bending moment interaction", "value": "failed",
                            "status": "failed"}
                        if configuration not in failed_configurations: failed_configurations.append(configuration)
                    else:
                        configuration.validator_parameters["skirt_load_interaction"] = {
                            "parameter": "axial compression and bending moment interaction", "value": "passed",
                            "status": "passed"}
                        if configuration not in passed_configurations: passed_configurations.append(configuration)

                pbar.update(1)

    else:
        print("No axial compression and bending moment interaction to validate.")

    print(f"Failed {len(failed_configurations)} configurations due to loads interaction."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


# ################### DOME ####################

def dome_yielding(configurations: list, loads: dict) -> tuple:
    """
    Validates all configurations for dome yielding.
    :param configurations: Configurations list object
    :param loads: Loads dictionary
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:

        for configuration in configurations:
            failed = False

            for dome in [configuration.forward_dome, configuration.aft_dome]:

                # Iterate over points on dome
                for z_d in np.linspace(0, dome.height, 100):

                    # Get stress state
                    s_mer, s_circ = dome.stress_state(z_d, loads)

                    # Check if yielding occurs
                    if not yield_criterion.yield_check(s_mer, s_circ, dome.material["yield_stress"]):
                        # Yielding occurs - configuration fails
                        configuration.validator_parameters["dome_yielding"] = {
                            "parameter": "dome yielding", "value": "failed", "status": "failed"}
                        if configuration not in failed_configurations: failed_configurations.append(configuration)
                        failed = True
                        break

                if failed: break

            if not failed:
                configuration.validator_parameters["dome_yielding"] = {
                    "parameter": "dome yielding", "value": "passed", "status": "passed"}
                if configuration not in passed_configurations: passed_configurations.append(configuration)

            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to dome yielding."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations


# ################### CONFIGURATION VALIDATION ####################

def ullage_check(configurations: list, config: object, ullage_lb: float, ullage_ub: float) -> tuple:
    """
    Validates all configurations for ullage.
    :param ullage_ub: Upper bound of ullage as a fraction of inner volume
    :param ullage_lb: Lower bound of ullage as a fraction of inner volume
    :param config: Configuration file object
    :param configurations: configurations list object
    :return: passed_configurations list, failed_configurations list
    """
    print("Validating configurations")

    failed_configurations = []
    passed_configurations = []

    with trange(len(configurations)) as pbar:
        for configuration in configurations:

            inner_volume = configuration.inner_volume
            inner_volume_ub = config.fluid_volume * (1 + ullage_ub)
            inner_volume_lb = config.fluid_volume * (1 + ullage_lb)

            if inner_volume_lb <= inner_volume <= inner_volume_ub:
                configuration.validator_parameters["ullage"] = {
                    "parameter": "inner volume", "value": inner_volume, "status": "passed"}
                passed_configurations.append(configuration)

            else:
                #print(f"Ullage failed for configuration {configuration.id} with inner volume {inner_volume}. Volume bounds are "
                #       f"{inner_volume_lb} and {inner_volume_ub}.")
                configuration.validator_parameters["ullage"] = {
                    "parameter": "inner volume", "value": inner_volume, "status": "failed"}
                failed_configurations.append(configuration)

            pbar.update(1)

    print(f"Failed {len(failed_configurations)} configurations due to ullage."
          f" {len(passed_configurations)} configurations passed.")

    return passed_configurations, failed_configurations
