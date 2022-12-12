"""
Functions for plotting a visual representation of a tank.
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_tank(configuration):
    """
    Plots a tank based on the configuration object passed in.
    """
    # Extract cylinder properties from configuration object
    cyl_or = configuration.cylinder.outer_radius
    cyl_ir = configuration.cylinder.inner_radius
    cyl_h = configuration.cylinder.height

    # Extract aft skirt properties from configuration object
    aft_skirt_or = configuration.aft_skirt.outer_radius
    aft_skirt_ir = configuration.aft_skirt.inner_radius
    aft_skirt_h = configuration.aft_skirt.height

    # Extract forward skirt properties from configuration object
    fwd_skirt_or = configuration.forward_skirt.outer_radius
    fwd_skirt_ir = configuration.forward_skirt.inner_radius
    fwd_skirt_h = configuration.forward_skirt.height

    # Create figure and axis, make scale equal
    fig, ax = plt.subplots()
    ax.set_aspect("equal")

    # Plot four corners of the inner radius of the cylinder
    ax.plot([-cyl_ir, cyl_ir], [0, 0], color="blue", label="Cylinder")
    ax.plot([-cyl_ir, cyl_ir], [cyl_h, cyl_h], color="blue", label="Cylinder")
    ax.plot([-cyl_ir, -cyl_ir], [0, cyl_h], color="blue", label="Cylinder")
    ax.plot([cyl_ir, cyl_ir], [0, cyl_h], color="blue", label="Cylinder")

    # Plot four corners of the outer radius of the cylinder
    ax.plot([-cyl_or, cyl_or], [0, 0], color="blue", label="Cylinder")
    ax.plot([-cyl_or, cyl_or], [cyl_h, cyl_h], color="blue", label="Cylinder")
    ax.plot([-cyl_or, -cyl_or], [0, cyl_h], color="blue", label="Cylinder")
    ax.plot([cyl_or, cyl_or], [0, cyl_h], color="blue", label="Cylinder")

    # Plot four corners of the inner radius of the aft skirt
    ax.plot([-aft_skirt_ir, aft_skirt_ir], [0, 0], color="red", label="Aft Skirt")
    ax.plot([-aft_skirt_ir, aft_skirt_ir], [-aft_skirt_h, -aft_skirt_h], color="red", label="Aft Skirt")
    ax.plot([-aft_skirt_ir, -aft_skirt_ir], [0, -aft_skirt_h], color="red", label="Aft Skirt")
    ax.plot([aft_skirt_ir, aft_skirt_ir], [0, -aft_skirt_h], color="red", label="Aft Skirt")

    # Plot four corners of the outer radius of the aft skirt
    ax.plot([-aft_skirt_or, aft_skirt_or], [0, 0], color="red", label="Aft Skirt")
    ax.plot([-aft_skirt_or, aft_skirt_or], [-aft_skirt_h, -aft_skirt_h], color="red", label="Aft Skirt")
    ax.plot([-aft_skirt_or, -aft_skirt_or], [0, -aft_skirt_h], color="red", label="Aft Skirt")
    ax.plot([aft_skirt_or, aft_skirt_or], [0, -aft_skirt_h], color="red", label="Aft Skirt")

    # Plot four corners of the inner radius of the forward skirt
    ax.plot([-fwd_skirt_ir, fwd_skirt_ir], [cyl_h, cyl_h], color="red", label="Forward Skirt")
    ax.plot([-fwd_skirt_ir, fwd_skirt_ir], [cyl_h + fwd_skirt_h, cyl_h + fwd_skirt_h], color="red", label="Forward Skirt")
    ax.plot([-fwd_skirt_ir, -fwd_skirt_ir], [cyl_h, cyl_h + fwd_skirt_h], color="red", label="Forward Skirt")
    ax.plot([fwd_skirt_ir, fwd_skirt_ir], [cyl_h, cyl_h + fwd_skirt_h], color="red", label="Forward Skirt")

    # Plot four corners of the outer radius of the forward skirt
    ax.plot([-fwd_skirt_or, fwd_skirt_or], [cyl_h, cyl_h], color="red", label="Forward Skirt")
    ax.plot([-fwd_skirt_or, fwd_skirt_or], [cyl_h + fwd_skirt_h, cyl_h + fwd_skirt_h], color="red", label="Forward Skirt")
    ax.plot([-fwd_skirt_or, -fwd_skirt_or], [cyl_h, cyl_h + fwd_skirt_h], color="red", label="Forward Skirt")
    ax.plot([fwd_skirt_or, fwd_skirt_or], [cyl_h, cyl_h + fwd_skirt_h], color="red", label="Forward Skirt")

    # Extract dome properties from configuration object
    for i, dome in enumerate([configuration.aft_dome, configuration.forward_dome]):

        if dome.type_ == "cassinian":
            # Get 100 points on the dome between -r and r
            dome_x_or = np.linspace(-dome.outer_radius, dome.outer_radius, 100)
            dome_x_ir = np.linspace(-dome.inner_radius, dome.inner_radius, 100)
            # Calculate the set of y values for the dome
            # Cassinian formula: y = sqrt((2 sqrt(a^2 n^4 (a^2 - x^2)))/n^4 + a^2/n^2 - x^2/n^2) where a = inner or outer
            # radius, n = dome n, x = x value
            dome_y_or = np.sqrt(
                (2 * np.sqrt(dome.outer_radius ** 2 * dome.n ** 4 * (dome.outer_radius ** 2 - dome_x_or ** 2))) / dome.n ** 4
                + dome.outer_radius ** 2 / dome.n ** 2
                - dome_x_or ** 2 / dome.n ** 2
            )
            dome_y_ir = np.sqrt(
                (2 * np.sqrt(dome.inner_radius ** 2 * dome.n ** 4 * (dome.inner_radius ** 2 - dome_x_ir ** 2))) / dome.n ** 4
                + dome.inner_radius ** 2 / dome.n ** 2
                - dome_x_ir ** 2 / dome.n ** 2
            )
            # Translate dome and flip if necessary
            if i == 0:
                # Aft dome
                dome_y_ir = -dome_y_ir
                dome_y_or = -dome_y_or
            elif i == 1:
                # Forward dome
                dome_y_ir = dome_y_ir + cyl_h
                dome_y_or = dome_y_or + cyl_h

            # Plot the dome
            ax.plot(dome_x_ir, dome_y_ir, color="green", label="Dome")
            ax.plot(dome_x_or, dome_y_or, color="green", label="Dome")

        elif dome.type_ == "semi-ellipsoidal":
            # Get 100 points on the dome between -r and r
            dome_x_or = np.linspace(-dome.outer_radius, dome.outer_radius, 100)
            dome_x_ir = np.linspace(-dome.inner_radius, dome.inner_radius, 100)

            # Calculate the set of y values for the dome
            # Ellipse formula (for semi-ellipsoidal dome): a!=0, y = (b sqrt(a^2 - x^2))/a, b!=0 where a = outer radius,
            # b = dome height, x = x value
            dome_y_or = (dome.height * np.sqrt(dome.outer_radius ** 2 - dome_x_or ** 2)) / dome.outer_radius
            dome_y_ir = ((dome.height - dome.thickness) * np.sqrt(dome.inner_radius ** 2 - dome_x_ir ** 2)) / dome.inner_radius

            # Translate dome and flip if necessary
            if i == 0:
                # Aft dome
                dome_y_ir = -dome_y_ir
                dome_y_or = -dome_y_or
            elif i == 1:
                # Forward dome
                dome_y_ir = dome_y_ir + cyl_h
                dome_y_or = dome_y_or + cyl_h

            # Plot the dome
            ax.plot(dome_x_ir, dome_y_ir, color="green", label="Dome")
            ax.plot(dome_x_or, dome_y_or, color="green", label="Dome")

        elif dome.type_ == "hemispherical":
            # Get 100 points on the dome between -r and r
            dome_x_or = np.linspace(-dome.outer_radius, dome.outer_radius, 100)
            dome_x_ir = np.linspace(-dome.inner_radius, dome.inner_radius, 100)

            # Calculate the set of y values for the dome
            # Hemispherical formula: y = sqrt(r^2 - x^2) where r = outer radius, x = x value
            dome_y_or = np.sqrt(dome.outer_radius ** 2 - dome_x_or ** 2)
            dome_y_ir = np.sqrt(dome.inner_radius ** 2 - dome_x_ir ** 2)

            # Translate dome and flip if necessary
            if i == 0:
                # Aft dome
                dome_y_ir = -dome_y_ir
                dome_y_or = -dome_y_or
            elif i == 1:
                # Forward dome
                dome_y_ir = dome_y_ir + cyl_h
                dome_y_or = dome_y_or + cyl_h

            # Plot the dome
            ax.plot(dome_x_ir, dome_y_ir, color="green", label="Dome")
            ax.plot(dome_x_or, dome_y_or, color="green", label="Dome")

    return fig, ax





