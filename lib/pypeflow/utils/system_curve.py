"""
## Calculate and draw the system curve of a flow path in a piping network

"""
from typing import Dict, Tuple, List
import numpy as np
from lib import quantities as qty
from lib.nummath.graphing2 import LineGraph
from lib.pypeflow.analysis.network import FlowPath


class SystemCurve:

    def __init__(self, R: float, R_units: Dict[str, str] = None, desired_units: Dict[str, str] = None):
        """
        Create *SystemCurve* object.

        **Parameters:**

        - `R`: (*float*) = (equivalent) hydraulic resistance of flow path
        - `R_units`: (*Dict[str, str]*) = the measuring units associated with the hydraulic resistance. Keys:
            + 'flow_rate'
            + 'pressure'
        - `desired_units`: (*Dict[str, str]*) = the desired measuring units in which to express the system curve
        """
        self._R = R
        self._V_unit: str = R_units['flow_rate'] if R_units else 'm^3/h'
        self._p_unit: str = R_units['pressure'] if R_units else 'bar'
        self._desired_units = desired_units if desired_units else {'flow_rate': 'm^3/s', 'pressure': 'Pa'}
        self._dp_stat: float = 0.0
        self._dp_elev: float = 0.0

    def set_static_head(self, p_stat: qty.Pressure):
        """Set static head (*quantities.Pressure*) of flow path."""
        self._dp_stat: float = p_stat(self._p_unit)

    def set_elevation_head(self, p_elev: qty.Pressure):
        """Set elevation head (*quantities.Pressure*) of flow path."""
        self._dp_elev: float = p_elev(self._p_unit)

    def create_system_curve(self, V_start: qty.VolumeFlowRate, V_end: qty.VolumeFlowRate, num: int = 50) \
            -> Tuple[List[float], List[float]]:
        """
        Calculate the system curve between an initial and final flow rate.

        **Parameters:**

        - `V_start`: (*quantities.VolumeFlowRate*) = first flow rate to put on the system curve
        - `V_end`: (*quantities.VolumeFlowRate*) = last flow rate to put on the system curve
        - `num`: (*int*) = number of points on the system curve (default = 50)

        **Returns:**
        Tuple with 1st element a list of the flow rates and 2nd element a list of the corresponding
        pressures, both expressed in the desired measuring units set at instantiation of the *SystemCurve*-object.

        """
        V_i = V_start(self._V_unit)
        V_f = V_end(self._V_unit)
        V_arr = np.linspace(V_i, V_f, num, endpoint=True)
        p_arr = self._R * V_arr ** 2 + self._dp_stat + self._dp_elev
        V_qty = [qty.VolumeFlowRate(V, self._V_unit) for V in V_arr]
        p_qty = [qty.Pressure(p, self._p_unit) for p in p_arr]
        V_sys = [V(self._desired_units['flow_rate']) for V in V_qty]
        p_sys = [p(self._desired_units['pressure']) for p in p_qty]
        return V_sys, p_sys

    def draw_system_curve(self, V_start: qty.VolumeFlowRate, V_end: qty.VolumeFlowRate, **kwargs) -> LineGraph:
        """
        Draw the calculated system curve.

        **Parameters:**

        - `V_start`: (*quantities.VolumeFlowRate*) = first flow rate to put on the system curve
        - `V_end`: (*quantities.VolumeFlowRate*) = last flow rate to put on the system curve
        - `kwargs`: optional keyword arguments
            + `fig_size`: (*Tuple[float, float]*) = the width and height of the figure in inches
            + `dpi`: (*int*) = dots per inch of the figure
            + `num`: (*int*) = number of calculated points to draw
            + `V_step`: (*quantities.VolumeFlowRate*) = step between ticks on the flow rate axis of the diagram
            + `V_max`: (*quantities.VolumeFlowRate*) = the maximum flow rate shown on the axis
            + `p_step`: (*quantities.Pressure*) = step between ticks on the pressure axis of the diagram
            + `p_max`: (*quantities.Pressure*) = maximum pressure shown on the axis

        **Returns:** (*nummath.graphing2.LineGraph*)<br>
        Call show() on the returned *LineGraph* object to show the diagram.
        """
        fig_size: Tuple[int, int] = kwargs.get('fig_size', (6, 4))
        dpi: int = kwargs.get('dpi', 96)
        num: int = kwargs.get('num', 50)
        V_step: qty.VolumeFlowRate = kwargs.get('V_step')
        V_max: qty.VolumeFlowRate = kwargs.get('V_max')
        p_step: qty.Pressure = kwargs.get('p_step')
        p_max: qty.Pressure = kwargs.get('p_max')
        V, p = self.create_system_curve(V_start, V_end, num)
        graph = LineGraph(fig_size=fig_size, dpi=dpi)
        graph.add_dataset(name="system curve", x1_data=V, y1_data=p)
        graph.x1.set_title(f'flow rate [{self._desired_units["flow_rate"]}]')
        if V_max is not None and V_step is not None:
            graph.x1.scale(
                lim_down=0.0,
                lim_up=V_max(self._desired_units['flow_rate']),
                step_size=V_step(self._desired_units['flow_rate'])
            )
        graph.y1.set_title(f'pressure [{self._desired_units["pressure"]}]')
        if p_max is not None and p_step is not None:
            graph.y1.scale(
                lim_down=0.0,
                lim_up=p_max(self._desired_units['pressure']),
                step_size=p_step(self._desired_units['pressure'])
            )
        return graph


def calculate_system_curves(flow_paths: List[FlowPath], V_wp: qty.VolumeFlowRate, dp_wp: qty.Pressure, **kwargs) \
        -> List[Tuple[List[float], List[float]]]:
    """
    Calculate the system curves of a list of `flow_paths` in a piping network with one common pump.
    To calculate the hydraulic resistance of the flow paths the working point of the pump must be known, specified
    by `V_wp` and `dp_wp`.

    **kwargs:**

    - `"unit_flow_rate"`: desired unit for flow rate
    - `"unit_pressure"`: desired unit for pressure
    - `"V_start"`: first flow rate to put on the system curve
    - `"V_end"`: last flow rate to put on the system curve
    """
    sys_curves = []
    for p in flow_paths:
        dp_loss = p.head_loss('bar')
        dp_loss += dp_wp('bar')
        R = dp_loss / V_wp('m^3/h') ** 2
        curve = SystemCurve(
            R=R,
            desired_units={
                'flow_rate': kwargs.get('unit_flow_rate', 'm^3/h'),
                'pressure': kwargs.get('unit_pressure', 'bar')
            }
        )
        curve.set_static_head(p.static_head)
        curve.set_elevation_head(p.elevation_head)
        V, dp = curve.create_system_curve(
            V_start=kwargs.get('V_start', qty.VolumeFlowRate(0.0, 'm^3/h')),
            V_end=kwargs.get('V_end', qty.VolumeFlowRate(5.0, 'm^3/h'))
        )
        sys_curves.append((V, dp))
    return sys_curves


def draw_curves(
        sys_curves: List[Tuple[List[float], List[float]]],
        pump_curves: List[Tuple[List[float], List[float]]],
        units: Tuple[str, str],
        **kwargs
):
    """
    Draw system curves and pump curves.

    kwargs
    ------
    + `fig_size`: (*Tuple[float, float]*) = the width and height of the figure in inches
    + `dpi`: (*int*) = dots per inch of the figure
    + `V_step`: (*float*) = step between ticks on the flow rate axis of the diagram
    + `V_max`: (*float*) = the maximum flow rate shown on the axis
    + `p_step`: (*float*) = step between ticks on the pressure axis of the diagram
    + `p_max`: (*float*) = maximum pressure shown on the axis
    + `working_point`: (*Tuple[float, float]*) = working point of the pump (shown as a red dot on the diagram)
    + `pump_curve_labels`: (*Tuple[str,...]) = legend labels for pump curves
    + `sys_curve_labels`: (*Tuple[str,...]) = legend labels for system curves
    """
    graph = LineGraph(
        fig_size=kwargs.get('fig_size', (10, 6)),
        dpi=kwargs.get('dpi', 96),
        figure_constructs=kwargs.get('figure_constructs')
    )
    sys_curve_labels = kwargs.get('sys_curve_labels')
    for i in range(len(sys_curves)):
        graph.add_dataset(
            name=f'{sys_curve_labels[i]}' if sys_curve_labels else f'path {i}',
            x1_data=sys_curves[i][0],
            y1_data=sys_curves[i][1]
        )
    pump_curve_labels = kwargs.get('pump_curve_labels')
    for i in range(len(pump_curves)):
        graph.add_dataset(
            name=f'{pump_curve_labels[i]}' if pump_curve_labels else f'pump {i}',
            x1_data=pump_curves[i][0],
            y1_data=pump_curves[i][1]
        )
    working_point = kwargs.get('working_point')
    if working_point:
        graph.add_dataset(
            name="working point",
            x1_data=working_point[0],
            y1_data=working_point[1],
            layout={'marker': 'o', 'linestyle': 'None', 'color': 'red'}
        )
    graph.x1.set_title(f"flow rate [{units[0]}]")
    graph.y1.set_title(f"pressure [{units[1]}]")
    V_max = kwargs.get('V_max')
    V_step = kwargs.get('V_step')
    p_max = kwargs.get('p_max')
    p_step = kwargs.get('p_step')
    if V_max is not None and V_step is not None:
        graph.x1.scale(lim_down=0.0, lim_up=V_max, step_size=V_step)
    if p_max is not None and p_step is not None:
        graph.y1.scale(lim_down=0.0, lim_up=p_max, step_size=p_step)
    graph.add_legend()
    return graph
