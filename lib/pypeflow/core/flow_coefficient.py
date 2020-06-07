"""
## Convert flow coefficients of fittings and valves
"""
import math
from lib import quantities as qty
from lib.pypeflow.core.fluids import Fluid, Water


class FlowCoefficient:
    """Class that groups class methods to convert between flow coefficient units."""

    rho: float = Water(15).density('kg/m^3')  # water density @ 15 °C

    @classmethod
    def Av_to_Kv(cls, Av: float) -> float:
        """
        Convert Av value (flow rate in m^3/s and pressure in Pa) (*float*) to Kv value (flow rate in m^3/h, pressure in
        bar and with density of water at 15 °C) (*float*).

        """
        Kv = Av * 3.6e5 * math.sqrt(10) / math.sqrt(cls.rho)
        return Kv

    @classmethod
    def Kv_to_Av(cls, Kv: float) -> float:
        """
        Convert Kv (*float*) to Av value (*float*).

        """
        Av = Kv * math.sqrt(cls.rho) / (3.6e5 * math.sqrt(10))
        return Av

    @classmethod
    def calc_Kv(cls, V: qty.VolumeFlowRate, dp: qty.Pressure, fluid: Fluid = Water(15.0)) -> float:
        """
        Calculate flow coefficient Kv (*float*) of a piping element if flow rate (*quantities.VolumeFlowRate*) and
        pressure drop (*quantities.Pressure*) are known.

        """
        V_base = V('m^3/h')
        dp_base = dp('bar')
        sg = fluid.density('kg/m^3') / cls.rho
        Kv = V_base / math.sqrt(dp_base / sg)
        return Kv

    @staticmethod
    def Kv_to_R(Kv: float):
        """
        Convert flow coefficient to hydraulic resistance
        """
        return 1 / (Kv ** 2)

    @staticmethod
    def R_to_Kv(R: float):
        """
        Convert hydraulic resistance to flow coefficient
        """
        return math.sqrt(1 / R)
