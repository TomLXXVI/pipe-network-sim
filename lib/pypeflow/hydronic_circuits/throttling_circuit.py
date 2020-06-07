"""
## DESIGNING A THROTTLING CIRCUIT

Design procedure:
1.  Determine design flow rate, diameters of the sections in the branch and calculate pressure losses in the branch
    without control valve and balancing valve. To determine the design flow rate use:
        - `T_sup_des` to set supply water temperature at design conditions
        - `T_ret_des` to set return water temperature at design conditions
        - `Q_des` to set the required heating capacity at design conditions
    If the design flow rate is already known, it can be set directly using `V_des`.
2.  Determine the Kvs of the balancing valve. Use:
        - `calc_bal_valve` to calculate a preliminary Kvs
        - `size_bal_valve` to set a commercially available Kvs
3.  Determine the Kvs of the control valve. Use:
        - `dP_br_wov` to set the pressure loss in the branch without balancing valve and without control valve
        - `calc_ctrl_valve` to calculate a preliminary Kvs
        - `size_ctrl_valve` to set a commercially available Kvs
4.  Balance the branch. Use:
        - `set_bal_valve` to get the calculated Kvr setting for the balancing valve
5.  Check control valve authority. Use:
        - `authority` to get the control valve authority after balancing.
"""

from typing import Optional
from lib import quantities as qty
from lib.pypeflow.core import Water
from lib.pypeflow.core import BalancingValve, ControlValve


class ThrottlingCircuit:

    def __init__(self, name: str):
        self.name = name
        self._Q_des: Optional[qty.Power] = None
        self._T_sup: Optional[qty.Temperature] = None
        self._T_ret: Optional[qty.Temperature] = None
        self._V_des: Optional[qty.VolumeFlowRate] = None
        self.water: Optional[Water] = None
        self._dP_br_wov: Optional[qty.Pressure] = None
        self._dP_cv: Optional[qty.Pressure] = None
        self._dP_bv: Optional[qty.Pressure] = None
        self._bal_valve = BalancingValve()
        self._ctrl_valve = ControlValve()

    @property
    def T_sup_des(self) -> qty.Temperature:
        """Get/set supply water temperature at design conditions."""
        return self._T_sup

    @T_sup_des.setter
    def T_sup_des(self, T_sup: qty.Temperature):
        self._T_sup = T_sup
        T_avg = (T_sup() + self._T_ret()) / 2.0
        self.water = Water(T_avg)

    @property
    def T_ret_des(self) -> qty.Temperature:
        """Get/set return water temperature at design conditions."""
        return self._T_ret

    @T_ret_des.setter
    def T_ret_des(self, T_ret: qty.Temperature):
        self._T_ret = T_ret
        T_avg = (self._T_sup() + T_ret()) / 2.0
        self.water = Water(T_avg)

    @property
    def Q_des(self) -> qty.Power:
        """Get/set heating capacity at design conditions."""
        return self._Q_des

    @Q_des.setter
    def Q_des(self, Q: qty.Power):
        self._Q_des = Q

    @property
    def V_des(self) -> qty.VolumeFlowRate:
        """Get/set required design flow rate."""
        rho = self.water.density()
        c = self.water.specific_heat()
        Q_des = self._Q_des()
        dT_des = self._T_sup() - self._T_ret()
        self._V_des = Q_des / (rho * c * dT_des)
        return qty.VolumeFlowRate(self._V_des)

    @V_des.setter
    def V_des(self, V_des: qty.VolumeFlowRate):
        self._V_des = V_des

    def calc_bal_valve(self, dp_100: qty.Pressure = qty.Pressure(3.0, 'kPa')) -> float:
        """Calculate preliminary Kvs of fully open balancing valve for a given pressure drop."""
        self._bal_valve = BalancingValve.create(
            fluid=self.water,
            flow_rate=self.V_des,
            dp_100=dp_100
        )
        self._dP_bv = dp_100
        return self._bal_valve.Kvs

    def size_bal_valve(self, Kvs: float):
        """Set commercially available Kvs of balancing valve."""
        self._bal_valve.Kvs = Kvs

    @property
    def dP_br_wov(self) -> qty.Pressure:
        """
        Get/set pressure drop in the branch without balancing valve and without control valve at design flow rate.
        """
        return self._dP_br_wov

    @dP_br_wov.setter
    def dP_br_wov(self, dP_br_wov: qty.Pressure):
        self._dP_br_wov = dP_br_wov

    def calc_ctrl_valve(self, auth_des=0.5) -> float:
        """
        Calculate preliminary Kvs of control valve to attain a certain valve authority at design conditions.
        Before calling this method, the Kvs of the balancing valve must have been specified and the pressure drop
        in the branch without balancing valve and control valve must have been set.
        """
        self._ctrl_valve = ControlValve.create(
            fluid=self.water,
            flow_rate=self.V_des,
            target_authority=auth_des,
            dp_crit_path=qty.Pressure(self._dP_br_wov() + self._dP_bv())
        )
        return self._ctrl_valve.Kvs

    def size_ctrl_valve(self, Kvs: float):
        """Set commercially available Kvs of control valve."""
        self._ctrl_valve.Kvs = Kvs
        self._dP_cv = self._ctrl_valve.pressure_drop

    def set_bal_valve(self, dP_br_avail: qty.Pressure) -> float:
        """
        Determine the balancing valve's Kvr setting if the available feed pressure across the branch is known.
        Call this method after the control valve's commercially available Kvs has been set.
        """
        dP_br_des = self._dP_br_wov() + self._dP_bv() + self._dP_cv()
        dP_sur = qty.Pressure(dP_br_avail() - dP_br_des)
        self._bal_valve.set_pressure_excess(dP_sur)
        self._dP_bv = self._bal_valve.pressure_drop
        return self._bal_valve.Kvr

    @property
    def authority(self) -> float:
        """
        Get the control valve's authority.
        """
        dP_cv_0 = self._dP_br_wov() + self._dP_bv() + self._dP_cv()
        return self._ctrl_valve.authority(qty.Pressure(dP_cv_0))

    @property
    def dP_br(self) -> qty.Pressure:
        """
        Get pressure drop across the branch with balancing valve and control valve included.
        """
        dP_br = self._dP_br_wov() + self._dP_bv() + self._dP_cv()
        return qty.Pressure(dP_br)

    @property
    def dP_bv(self) -> qty.Pressure:
        """
        Get pressure drop across the balancing valve.
        """
        return self._dP_bv

    @property
    def dP_cv(self) -> qty.Pressure:
        """
        Get pressure drop across the control valve.
        """
        return self._dP_cv
