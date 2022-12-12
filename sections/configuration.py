"""
Configuration class. Inherits parameters from the sections, and calculates the configuration parameters.
"""

from formulas.bsc_ae import geometry

import uuid
import warnings


class Configuration:
    def __init__(self,
                 forward_dome,
                 aft_dome,
                 forward_skirt,
                 aft_skirt,
                 cylinder,
                 internal_pressure,
                 fluid,
                 fluid_volume):
        """
        Configuration object, containing all relevant parameters. Configuration coordinate system (z) is defined with the
        origin at the aft tank dome's center. The z axis moves axially along the tank.
        :param forward_dome: Dome object
        :param aft_dome: Dome object
        :param forward_skirt: Skirt object
        :param aft_skirt: Skirt object
        :param cylinder: Cylinder object
        :param internal_pressure: in Pa
        :param fluid: dictionary object from fluids database (databases/fluids.py)
        :param fluid_volume: in m^3
        """

        self.id = uuid.uuid4()
        # The line above creates a unique ID for the configuration. This is used to identify the configuration in the
        # set.
        self.forward_dome = forward_dome
        self.aft_dome = aft_dome
        self.forward_skirt = forward_skirt
        self.aft_skirt = aft_skirt
        self.cylinder = cylinder
        self.fluid = fluid
        self.fluid_volume = fluid_volume
        self.internal_pressure = internal_pressure
        self.validator_parameters = {}  # Stores parameters from the validators. {"validator_name": {"parameter":
        # parameter_name, "value": value, "status": "passed/failed"}}

    @property
    def forward_dome(self):
        return self._forward_dome

    @forward_dome.setter
    def forward_dome(self, value):
        # Check if forward and aft domes are the same. If so, raise error. Else, set the forward dome.
        if (not hasattr(self, "_aft_dome")) or (hasattr(self, "_aft_dome") and self.aft_dome != value):
            self._forward_dome = value
            value.parent_configuration = self
        else:
            raise ValueError("Forward and aft dome objects cannot be the same.")

    @property
    def aft_dome(self):
        return self._aft_dome

    @aft_dome.setter
    def aft_dome(self, value):
        # Check if forward and aft domes are the same. If so, raise error. Else, set the aft dome.
        if (not hasattr(self, "_forward_dome")) or (hasattr(self, "_forward_dome") and self.forward_dome != value):
            self._aft_dome = value
            value.parent_configuration = self
        else:
            raise ValueError("Forward and aft dome objects cannot be the same.")

    @property
    def forward_skirt(self):
        return self._forward_skirt

    @forward_skirt.setter
    def forward_skirt(self, value):
        self._forward_skirt = value
        value.parent_configuration = self

    @property
    def aft_skirt(self):
        return self._aft_skirt

    @aft_skirt.setter
    def aft_skirt(self, value):
        self._aft_skirt = value
        value.parent_configuration = self

    @property
    def cylinder(self):
        return self._cylinder

    @cylinder.setter
    def cylinder(self, value):
        self._cylinder = value
        value.parent_configuration = self

    @property
    def height(self) -> float:
        # Ignores skirts and welds for now
        return self.cylinder.height + self.forward_skirt.height + self.aft_skirt.height

    @property
    def mass(self) -> float:
        # Ignores skirts and welds for now
        return self.forward_dome.mass + self.cylinder.mass + self.aft_dome.mass + self.forward_skirt.mass + \
                  self.aft_skirt.mass

    @property
    def inner_volume(self) -> float:
        # Ignores skirts and welds for now
        return self.forward_dome.inner_volume + self.cylinder.inner_volume + self.aft_dome.inner_volume

    @property
    def outer_volume(self) -> float:
        # Ignores skirts and welds for now
        return self.forward_dome.outer_volume + self.cylinder.outer_volume + self.aft_dome.outer_volume

    def effective_internal_pressure(self, z: float, acceleration: float) -> float:
        """
        Finds the effective internal pressure at a given height. Takes values from domes and cylinder,
        :param z: Height in m
        :return:
        """
        if z < 0:
            raise ValueError("Height cannot be negative.")
        elif z < self.aft_dome.height:
            return self.aft_dome.effective_internal_pressure(z_d=self.aft_dome.height - z, acceleration=acceleration)
        elif z < self.aft_dome.height + self.cylinder.height:
            return self.cylinder.effective_internal_pressure(z_c=z - self.aft_dome.height,
                                                             acceleration=acceleration)
        elif z <= self.aft_dome.height + self.cylinder.height + self.forward_dome.height:
            return self.forward_dome.effective_internal_pressure(z_d=z - self.aft_dome.height - self.cylinder.height,
                                                                 acceleration=acceleration)
        else:
            raise ValueError(f"Height {z} is outside of tank.")
