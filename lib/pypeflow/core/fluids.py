"""
## Definitions of fluids used in piping networks
"""
from typing import Dict, Type
from CoolProp.CoolProp import PropsSI
from lib import quantities as qty


class Fluid:
    """
    Base class that defines the available properties for fluids that are needed in piping network calculations.
    Fluid properties are retrieved using the third-party package [CoolProp](http://www.coolprop.org/).
    """
    fluid = None

    def __init__(self, T: float, p_gauge: float = 0.0):
        """
        Create Fluid object.

        **Parameters:**

        - `T`: *float*<br>
        The temperature of the fluid in °C.
        - `p_gauge`: *float*<br>
        The gauge pressure of the fluid in Pa (optional, default is 0 Pa).
        """
        T_abs = 273.15 + T
        P_abs = 101325.0 + p_gauge
        self._density = PropsSI('D', 'T', T_abs, 'P', P_abs, self.fluid)            # [kg/m^3]
        self._dynamic_viscosity = PropsSI('V', 'T', T_abs, 'P', P_abs, self.fluid)  # [Pa.s]
        self._kinematic_viscosity = self._dynamic_viscosity / self._density         # [m^2/s]
        self._specific_heat = PropsSI('C', 'T', T_abs, 'P', P_abs, self.fluid)      # J/(kg.K)

    @property
    def density(self) -> qty.MassDensity:
        """
        Get the mass density (*quantities.MassDensity*) of the fluid.

        """
        return qty.MassDensity(self._density)

    @property
    def kinematic_viscosity(self) -> qty.KinematicViscosity:
        """
        Get the kinematic viscosity (*quantities.KinematicViscosity*) of the fluid.

        """
        return qty.KinematicViscosity(self._kinematic_viscosity)

    @property
    def dynamic_viscosity(self) -> qty.DynamicViscosity:
        """
        Get the dynamic viscosity (*quantities.DynamicViscosity*) of the fluid.

        """
        return qty.DynamicViscosity(self._dynamic_viscosity)

    @property
    def specific_heat(self) -> qty.SpecificHeat:
        """Get the specific heat capacity (*tsq.DynamicViscosity*) of the fluid."""
        return qty.SpecificHeat(self._specific_heat)


class Water(Fluid):
    """Child class of Fluid that defines water."""
    fluid = 'Water'


class Air(Fluid):
    """Child class of Fluid that defines air."""
    fluid = 'Air'


FLUIDS: Dict[str, Type[Fluid]] = {
    'water': Water,
    'air': Air
}
"""Dictionary that holds the available types of Fluid"""
