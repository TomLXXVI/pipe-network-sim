from typing import Type, Union, Optional
import math
from CoolProp.CoolProp import PropsSI
from lib.quantities.base import Quantity


def create_from_str(type_: Type[Quantity], value: str, unit: Optional[str] = None,
                    fall_back_value: Optional[float] = 0.0) -> Optional[Quantity]:
    try:
        value = float(value)
    except ValueError:
        if fall_back_value is not None:
            value = fall_back_value
        else:
            return None
    return type_(value, unit)


class Length(Quantity):
    base_unit = 'm'
    units = {'m': 1e0, 'mm': 1e3}


class Area(Quantity):
    base_unit = 'm^2'
    units = {'m^2': 1e0, 'mm^2': 1e6}


class Velocity(Quantity):
    base_unit = 'm/s'
    units = {'m/s': 1e0, 'km/h': 3.6e0}


class Pressure(Quantity):
    # mass density of water @ 10 °C and standard atmospheric pressure
    rho = PropsSI('D', 'T', 273.15 + 10.0, 'P', 101325.0, 'Water')  # kg/m^3
    g = 9.81  # m/s^2
    base_unit = 'Pa'
    units = {'Pa': 1e0, 'kPa': 1e-3, 'bar': 1e-5, 'MPa': 1e-6, 'm': 1.0 / (rho * g)}

    @classmethod
    def convert(cls, src_value: float, src_unit: str, des_unit: str):
        if src_unit == "m":
            base_value = src_value * cls.rho * cls.g
        else:
            src_cf = cls.units.get(src_unit)
            base_value = src_value / src_cf
        if des_unit == "m":
            des_value = base_value / (cls.rho * cls.g)
        else:
            des_cf = cls.units.get(des_unit)
            des_value = base_value * des_cf
        return des_value


class Angle(Quantity):
    base_unit = 'rad'
    units = {'rad': 1e0, 'deg': math.pi / 180.0}


class KinematicViscosity(Quantity):
    base_unit = 'm^2/s'
    units = {'m^2/s': 1.0}


class DynamicViscosity(Quantity):
    base_unit = 'Pa.s'
    units = {'Pa.s': 1.0}


class MassDensity(Quantity):
    base_unit = 'kg/m^3'
    units = {'kg/m^3': 1.0}


class Power(Quantity):
    base_unit = 'W'
    units = {'W': 1.0, 'kW': 1.0e-3}


class SpecificHeat(Quantity):
    base_unit = 'J/(kg.K)'
    units = {'J/(kg.K)': 1.0, 'kJ/(kg.K)': 1.0e-3}


class Temperature(Quantity):
    base_unit = '°C'
    units = {'°C': 0.0, 'K': 273.15}

    @classmethod
    def convert(cls, src_value: float, src_unit: str, des_unit: str):
        if src_unit == '°C' and des_unit == 'K':
            return src_value + cls.units['K']
        elif src_unit == 'K' and des_unit == '°C':
            return src_value - cls.units['K']
        else:
            return src_value


class VolumeFlowRate(Quantity):
    base_unit = 'm^3/s'
    units = {'m^3/s': 1e0, 'L/min': 6e4, 'm^3/h': 3.6e3, 'L/s': 1e3}

    def to_mass(self, T=Temperature(10.0), P_atm=Pressure(101325.0), fluid='Water') -> 'MassFlowRate':
        rho = PropsSI('D', 'T', T('K'), 'P', P_atm('Pa'), fluid)
        return MassFlowRate(self() * rho)

    @classmethod
    def from_mass(cls, m: 'MassFlowRate', T=Temperature(10.0), P_atm=Pressure(101325.0),
                  fluid='Water') -> 'VolumeFlowRate':
        rho = PropsSI('D', 'T', T('K'), 'P', P_atm('Pa'), fluid)
        return cls(m() / rho)


class MassFlowRate(Quantity):
    base_unit = 'kg/s'
    units = {'kg/s': 1.0}

    def to_volume(self, T=Temperature(10.0), P_atm=Pressure(101325.0), fluid='Water') -> VolumeFlowRate:
        rho = PropsSI('D', 'T', T('K'), 'P', P_atm('Pa'), fluid)
        return VolumeFlowRate(self() / rho)

    @classmethod
    def from_volume(cls, V: VolumeFlowRate, T=Temperature(10.0), P_atm=Pressure(101325.0),
                    fluid='Water') -> 'MassFlowRate':
        rho = PropsSI('D', 'T', T('K'), 'P', P_atm('Pa'), fluid)
        return cls(rho * V())
