"""
Dome section class.
"""

from formulas.bsc_ae import pressure, geometry
from formulas.flugge import cassinian_stress, ellipse_stress

import numpy as np

dome_types = ["hemispherical", "semi-ellipsoidal", "cassinian"]


class Dome:
    def __init__(self,
                 type_: str,
                 outer_radius: float,
                 thickness: float,
                 material: dict,
                 height: float = None,
                 n: float = None):
        """
        Dome object, containing all relevant parameters. Dome coordinate system is defined with the origin at the center,
        in the plane of the connection with the cylinder. The z_d axis moves along the dome's axis in the direction of
        decreasing radius. Hence, for the aft dome, the coordinate system points in the direction opposite to the
        global coordinate system for the tank configuration.

        :param type_: Type of dome.
        :param outer_radius: in m
        :param thickness: in m
        :param material: dictionary object from materials database (databases/materials.py)
        :param height: in m
        :param n: Affine parameter for Cassinian dome. If not provided, will be calculated from the dome's geometry.
        """

        self.type_ = type_
        self.outer_radius = outer_radius
        self.thickness = thickness
        self.material = material
        self.height = height
        self.n = n
        self.parent_configuration = None

        if not self.is_thin_walled: raise ValueError(
            f"Dome is not thin walled (OD: {self.outer_diameter}, H: {self.height}), t: {self.thickness})")

        if self.height / self.outer_radius > 0.707 and self.type_ == "semi-ellipsoidal": raise ValueError(
            f"Semi ellipsoidal ratio exceeds critical value ({self.height / self.outer_radius} > 0.707).")

    @property
    def type_(self):
        return self._type_

    @type_.setter
    def type_(self, value):
        if value not in dome_types:
            raise ValueError(f"Dome type not supported. Supported types are: {dome_types}")
        else:
            self._type_ = value

    @property
    def inner_radius(self) -> float:
        return self.outer_radius - self.thickness

    @property
    def inner_diameter(self) -> float:
        return 2 * self.inner_radius

    @property
    def outer_diameter(self) -> float:
        return 2 * self.outer_radius

    @property
    def height(self) -> float:
        if self.type_ == "cassinian":
            return geometry.cassinian_h(radius=(self.outer_radius + self.inner_radius) / 2,
                                        n=self.n)
        else:
            return self._height

    @height.setter
    def height(self, value):
        # Set height to outer radius for hemispherical dome.
        if self.type_ == "hemispherical" and value is not None and value != self.outer_radius:
            raise ValueError("Height of hemispherical dome must be equal to outer radius.")
        elif self.type_ == "hemispherical" and value is None:
            self._height = self.outer_radius
        elif self.type_ == "semi-ellipsoidal" and value is None:
            raise ValueError("Height must be provided for semi-ellipsoidal dome.")
        elif self.type_ == "cassinian" and value is not None:
            raise ValueError("Height must not be provided for Cassinian dome.")
        else:
            self._height = value

    @property
    def n(self) -> float:
        return self._n

    @n.setter
    def n(self, value):
        # Set n only for Cassinian dome.
        if value is None and self.type_ == "cassinian":
            raise ValueError("n must be provided for Cassinian dome.")
        else:
            self._n = value

    @property
    def inner_volume(self) -> float:
        if self.type_ == "hemispherical":
            return geometry.hemisphere_V(self.inner_radius)
        elif self.type_ == "semi-ellipsoidal":
            return geometry.semi_ellipsoid_V(self.inner_radius, self.height)
        elif self.type_ == "cassinian":
            return geometry.cassinian_V(self.inner_radius, self.n)

    @property
    def outer_volume(self) -> float:
        if self.type_ == "hemispherical":
            return geometry.hemisphere_V(self.outer_radius)
        elif self.type_ == "semi-ellipsoidal":
            return geometry.semi_ellipsoid_V(self.outer_radius, self.height)
        elif self.type_ == "cassinian":
            return geometry.cassinian_V(self.outer_radius, self.n)

    @property
    def mass(self) -> float:
        material_volume = self.outer_volume - self.inner_volume
        return material_volume * self.material["density"]

    @property
    def is_thin_walled(self) -> bool:
        return (self.thickness <= 0.1 * self.outer_diameter) and (self.thickness <= 0.1 * self.height)

    def z_d_to_section_radius(self, z_d: float) -> float:
        """
    `   Finds the sectional radius of the dome section at a given height, using the local dome coordinate system (z_d).
        :param z_d: in m
        :return: radius in m
        """
        if self.type_ == "hemispherical":
            return np.sqrt(self.outer_radius ** 2 - z_d ** 2)
        elif self.type_ == "semi-ellipsoidal":
            return self.outer_radius * np.sqrt(1 - (z_d / self.height) ** 2)
        elif self.type_ == "cassinian":
            return geometry.cassinian_z_d_to_section_radius(z_d=z_d,
                                                            n=self.n,
                                                            radius=(self.inner_radius + self.outer_radius) / 2)

    def stress_state(self, z_d: float, loads: dict) -> tuple:
        """
        Finds the stress state of the dome at a given height, using the local dome coordinate system (z_d).
        :param z_d: in m
        :return: (sigma_meridional, sigma_circumferential)
        """
        # Hydrostatic pressure for aft dome
        if self.parent_configuration.aft_dome == self:
            hydrostatic_pressure = pressure.find_hydrostatic_pressure(rho=self.parent_configuration.fluid["density"],
                                                                      h=self.parent_configuration.height,
                                                                      a=loads["a_axial"])
        else:
            hydrostatic_pressure = 0

        if self.type_ == "hemispherical":
            s_press = pressure.find_stress_in_dome_from_internal_pressure(p=self.parent_configuration.internal_pressure + hydrostatic_pressure,
                                                                          r=(self.inner_radius + self.outer_radius) / 2,
                                                                          t=self.thickness)
            return s_press, s_press

        elif self.type_ == "semi-ellipsoidal":
            local_radius = self.z_d_to_section_radius(z_d=z_d)
            phi = np.arctan(local_radius / z_d)
            s_mer_press, s_circ_press = ellipse_stress.stress_position(
                effective_internal_pressure=self.parent_configuration.internal_pressure + hydrostatic_pressure,
                phi=phi,
                major_axis=(self.inner_radius + self.outer_radius) / 2,
                minor_axis=self.height,
                thickness=self.thickness)
            return s_mer_press, s_circ_press

        elif self.type_ == "cassinian":
            x = self.z_d_to_section_radius(z_d=z_d)
            s_mer_press, s_circ_press = cassinian_stress.stress_position(
                effective_internal_pressure=self.parent_configuration.internal_pressure + hydrostatic_pressure,
                z=x,
                x=z_d,
                radius=(self.inner_radius + self.outer_radius) / 2,
                n=self.n,
                thickness=self.thickness)
            return s_mer_press, s_circ_press
