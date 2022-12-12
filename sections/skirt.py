# Skirt section class.

from formulas.bsc_ae import geometry


class Skirt:
    def __init__(self,
                 outer_radius: float,
                 thickness: float,
                 extension: float,
                 material: dict,
                 stiffening_elastic_constants: tuple,
                 stiffening_volume: float = 0
                 ):
        """
        Skirt object, containing all relevant parameters. Skirt coordinate system is defined with the origin at the
        connection between the skirt and the cylinder. The z_s axis moves along the skirt's axis in the direction
        from connection to edge.

        :param outer_radius: Outer radius of the skirt
        :param thickness: Thickness of the skirt
        :param extension: Extension of the skirt past the dome
        :param material: Material dictionary, containing the following keys: 'E', 'v', 'G', 'rho'
        :param stiffening_elastic_constants: Tuple containing the elastic constants for the stiffening method used, in the
        form (E_x, E_y, E_xy, G_xy, C_x, C_y, C_xy, K_xy, D_x, D_y, D_xy, F_x, F_y, F_xy, H_x, H_y, H_xy, M_x, M_y, M_xy)
        :param stiffening_volume: Volume of the stiffeners
        """

        # These are the parameters that are passed to the class:
        self.outer_radius = outer_radius
        self.thickness = thickness
        self.extension = extension
        self.material = material
        self._stiffening_elastic_constants = stiffening_elastic_constants
        self.stiffening_volume = stiffening_volume
        self.parent_configuration = None

    @property
    def stiffening_elastic_constants(self):
        return self._stiffening_elastic_constants

    @stiffening_elastic_constants.setter
    def stiffening_elastic_constants(self, value):
        self._stiffening_elastic_constants = value

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
    def inner_volume(self) -> float:
        return geometry.cylinder_V(self.inner_radius, self.height)

    @property
    def outer_volume(self) -> float:
        return geometry.cylinder_V(self.outer_radius, self.height)

    @property
    def section_Ixx(self) -> float:
        # TODO: Adjust for stiffeners
        return geometry.cylindrical_shell_I(self.outer_radius, self.thickness)

    @property
    def section_Iyy(self) -> float:
        # TODO: Adjust for stiffeners
        return self.section_Ixx

    @property
    def mass(self) -> float:
        # TODO: Adjust for stiffeners
        material_volume = self.outer_volume - self.inner_volume
        material_volume += self.stiffening_volume
        return self.material["density"] * material_volume

    @property
    def sectional_area(self) -> float:
        # TODO: Adjust for stiffeners
        return geometry.cylindrical_shell_A(self.outer_radius, self.thickness)

    @property
    def height(self) -> float:
        # Check if this is aft skirt
        if self.parent_configuration.aft_skirt is self:
            return self.parent_configuration.aft_dome.height + self.extension
        elif self.parent_configuration.forward_skirt is self:
            return self.parent_configuration.forward_dome.height + self.extension
