from lib import quantities as qty
from pypeflow.hydronic_circuits.throttling_circuit import ThrottlingCircuit, BalancingValve


class DivergingCircuit(ThrottlingCircuit):

    def __init__(self, name):
        super().__init__(name)
        self._bal_valve_byp = BalancingValve()

    def set_bypass(self, Kvs: float) -> float:
        """
        Get the Kvr setting of the balancing valve in the bypass.
        """
        self._bal_valve_byp = BalancingValve.create_from_str(
            fluid=self.water,
            flow_rate=self.V_des,
            dp_100=qty.Pressure(3.0, 'kPa')
        )
        self._bal_valve_byp.Kvs = Kvs
        dp_sur = self._dP_br_wov() - self._bal_valve_byp.pressure_drop()
        self._bal_valve_byp.set_pressure_excess(qty.Pressure(dp_sur))
        return self._bal_valve_byp.Kvr
