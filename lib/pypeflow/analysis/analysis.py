"""
##  User interface for doing network flow analysis using the Hardy Cross method
"""
from typing import Dict, Optional, Tuple, List
import csv
import pandas as pd
from lib import quantities as qty
from lib.pypeflow.analysis.network import Network
from lib.pypeflow.core.fluids import FLUIDS
from lib.pypeflow.core.pipe_schedules import PIPE_SCHEDULES


class Analyzer:
    """Class that encapsulates the user interface methods for network flow analysis"""
    network: Network = Network()
    """Reference to the *Network* object"""
    units: Dict[str, str] = {
        'length': 'm',
        'diameter': 'mm',
        'flow_rate': 'L/s',
        'pressure': 'bar',
        'velocity': 'm/s'
    }
    """Dictionary that holds the measuring units in which the quantity values are expressed"""

    @classmethod
    def set_units(cls, units: Dict[str, str]):
        """
        Set the measuring SI-units of the quantities that will be used as input and that will be returned as output.

        Parameter `units` is a dictionary (*Dict[str, str]*) that can contain the following keys with corresponding
        measuring units assigned to them:

        - *'length'* (default value = *'m'*)
        - *'diameter'* (default value = *'mm'*)
        - *'flow_rate'* (default value = *'L/s'*)
        - *'pressure'* (default value = *'bar'*)
        - *'velocity'* (default value = *'m/s'*)

        """
        cls.units.update(units)

    @classmethod
    def create_network(cls, **kwargs):
        """
        Create the *Network* object.

        **kwargs:**

        - `start_node_id`: (*str*) = start node id of the network
        - `end_node_id`: (*str*) = end node id of the network
        - `fluid`: (*str*) = the fluid that flows in the network (default = *'water'*)
        - `fluid_temperature`: (*float*) = the fluid temperature [°C]
        - `pipe_schedule`: (*str*) = the pipe schedule of the sections (default = *'pipe_schedule_40'*)

        """
        start_node_id: str = kwargs.get('start_node_id', '')
        end_node_id: str = kwargs.get('end_node_id', '')
        fluid_str: str = kwargs.get('fluid', 'water')
        fluid_temperature: float = kwargs.get('fluid_temperature', 10.0)
        sch_str: str = kwargs.get('pipe_schedule', 'pipe_schedule_40')

        fluid = cls._create_fluid(fluid_str, fluid_temperature)
        pipe_schedule = cls._create_pipe_schedule(sch_str)
        cls.network = Network.create(
            start_node_id=start_node_id,
            end_node_id=end_node_id,
            fluid=fluid,
            pipe_schedule=pipe_schedule
        )

    @classmethod
    def _create_fluid(cls, fluid: str, temperature: float):
        try:
            fl = FLUIDS[fluid.lower()]
        except KeyError:
            raise KeyError(f'Fluid {fluid} unknown.')
        else:
            return fl(temperature)

    @classmethod
    def _create_pipe_schedule(cls, pipe_schedule: str):
        try:
            sch = PIPE_SCHEDULES[pipe_schedule.lower()]
        except KeyError:
            raise KeyError(f'Pipe schedule {pipe_schedule} unknown.')
        else:
            return sch

    @classmethod
    def configure_network(cls, file_path: str):
        """
        Configure network via a network configuration file (.csv-file). Parameter `file_path` (*str*) is the file path
        to this configuration file. The configuration data is organised in a table. Each row contains the configuration
        data of a pipe section in the network. Each row has the following fields (columns) in the order as
        mentioned here:

        0. id of the loop to which the section belongs
        1. id of the section in the network
        2. id of the start node of the section
        3. id of the end node of the section
        4. nominal diameter of the section
        5. length of the section
        6. sum of the resistance coefficients of fittings/valves in the section
        7. pump coefficient a0 in the equation dp_pump = a0 + a1 * V + a2 * V ** 2 (leave empty if no pump is present
        in the section)
        8. pump coefficient a1
        9. pump coefficient a2
        10. fixed pressure difference between start and end node of the section (only in case of pseudo section, leave
        empty if the section is not a pseudo section)
        11. flow rate through the section (leave empty in case of a pseudo section)

        Fixed pressure differences and flow rates in sections must carry a sign with reference to the positive loop
        sense (by convention clockwise sense).
        """
        with open(file_path) as f:
            reader = csv.reader(f)
            for i, row in enumerate(reader):
                if i == 0:
                    continue
                else:
                    cls.network.add_section(
                        loop_id=row[0],
                        section_id=row[1],
                        start_node_id=row[2],
                        end_node_id=row[3],
                        nominal_diameter=qty.create_from_str(qty.Length, row[4], cls.units['diameter'], 0.0),
                        length=qty.create_from_str(qty.Length, row[5], cls.units['length'], 0.0),
                        zeta=cls._set_zeta(row[6], 0.0),
                        pump_curve=cls._set_pump_curve(row[7], row[8], row[9]),
                        dp_fixed=qty.create_from_str(qty.Pressure, row[10], cls.units['pressure'], None),
                        flow_rate=qty.create_from_str(qty.VolumeFlowRate, row[11], cls.units['flow_rate'], 0.0)
                    )

    @staticmethod
    def _set_pump_curve(a0: str, a1: str, a2: str) -> Optional[Tuple[float, float, float]]:
        try:
            a0 = float(a0)
            a1 = float(a1)
            a2 = float(a2)
        except ValueError:
            return None
        else:
            return a0, a1, a2

    @staticmethod
    def _set_zeta(value: str, fallback_value: float = 0.0) -> Optional[float]:
        try:
            value = float(value)
        except ValueError:
            return fallback_value
        else:
            return value

    @classmethod
    def check_flow_balance(cls, V_ext: Dict[str, Dict[str, Optional[List[qty.VolumeFlowRate]]]] = None) \
            -> Dict[str, qty.VolumeFlowRate]:
        """
        Check the flow balances at the nodes of the network.
        External flow rates that enter or leave at network nodes are added through `V_ext`. It is a dictionary
        with the following syntax:
        ```V_ext = {<node_id>: {"in": [<qty.VolumeFlowRate>,...], "out": [<qty.VolumeFlowRate>,...]}}```
        The method returns a dictionary `checks` like {<node_id>: <qty.VolumeFlowRate>, ...}
        The value of given `<node_id>`-key is the difference between the incoming flow rate and the outgoing flow rate.
        This net flow rate should be zero in order to respect the physical law of continuity. If positive, incoming
        flow rate to the node is greater than flow rate that is leaving the node. On the other hand, if negative,
        outgoing flow rate is greater than incoming flow rate.
        """
        return cls.network.check_flow_balance(V_ext)

    @classmethod
    def solve(cls, error: float = 1.0e-3, i_max: int = 30):
        """
        Solve the piping network for flow rates and pressure drops.

        **Parameters:**

        - `error`: (*float*) = allowable deviation from zero for the pressure drop around each loop
        - `i_max`: (*int*) = the maximum number of iterations to find a solution within the given error tolerance

        If no solution within the given fault tolerance is found after maximum number of iterations an *OverflowError*
        exception is raised.

        """
        cls.network.solve(error, i_max)

    @classmethod
    def get_network(cls) -> pd.DataFrame:
        """Return the solved network as a Pandas DataFrame."""
        keys = [
            'loop_id',
            'section_id',
            'start_node_id',
            'end_node_id',
            f'length [{cls.units["length"]}]',
            f'diameter [{cls.units["diameter"]}]',
            'zeta',
            f'flow_rate [{cls.units["flow_rate"]}]',
            f'velocity [{cls.units["velocity"]}]',
            f'pressure_drop [{cls.units["pressure"]}]'
        ]
        d = {k: [] for k in keys}
        for loop in cls.network.loops:
            for section in loop.sections.values():
                d[keys[0]].append(loop.id)
                d[keys[1]].append(section.id)
                d[keys[2]].append(section.start_node.id)
                d[keys[3]].append(section.end_node.id)
                d[keys[4]].append(section.length(cls.units['length'], 3))
                d[keys[5]].append(section.nominal_diameter(cls.units['diameter'], 3))
                d[keys[6]].append(section.zeta)
                d[keys[7]].append(section.sign * section.flow_rate(cls.units['flow_rate'], 3))
                d[keys[8]].append(section.sign * section.velocity(cls.units['velocity'], 3))
                d[keys[9]].append(section.pressure_drop(cls.units['pressure'], 3))
        return pd.DataFrame(d)

    @classmethod
    def configure_network_from_df(cls, df: pd.DataFrame, clear=False):
        """
        Configure network from a Pandas DataFrame object. If clear is set to True the internal network objects will
        be cleared first; this may be needed in case the network needs to be reconfigured.
        """
        if clear is True:
            cls.network.clear()
        for _, row in df.iterrows():
            row = cls._check_row(row)
            cls.network.add_section(
                loop_id=row[0],
                section_id=row[1],
                start_node_id=row[2],
                end_node_id=row[3],
                nominal_diameter=qty.Length(row[4], cls.units['diameter']),
                length=qty.Length(row[5], cls.units['length']),
                zeta=row[6],
                pump_curve=row[7],
                dp_fixed=row[8],
                flow_rate=qty.VolumeFlowRate(row[9], cls.units['flow_rate'])
            )

    @classmethod
    def _check_row(cls, row):
        row_ = [row.iloc[i] for i in (0, 1, 2, 3)]
        for i in (4, 5, 6):
            if pd.isna(row.iloc[i]):
                row_.append(0.0)
            else:
                row_.append(row.iloc[i])
        if (pd.isna(row.iloc[7:10])).any():
            row_.append(None)
        else:
            row_.append((row.iloc[7], row.iloc[8], row.iloc[9]))
        if pd.isna(row.iloc[10]):
            row_.append(None)
        else:
            row_.append(qty.Pressure(row.iloc[10], cls.units['pressure']))
        if pd.isna(row.iloc[11]):
            row_.append(0.0)
        else:
            row_.append(row.iloc[11])
        return row_

    @classmethod
    def get_paths(cls) -> pd.DataFrame:
        """
        Get the flow paths in the solved network, returned as a Pandas DataFrame.
        For each flow path is returned:

        - velocity head
        - elevation head
        - dynamic head
        - static head

        """
        keys = [
            'path',
            f'dp,vel [{cls.units["pressure"]}]',
            f'dp,elev [{cls.units["pressure"]}]',
            f'dp,loss (net) [{cls.units["pressure"]}]',
            f'dp,stat [{cls.units["pressure"]}]',
        ]
        d = {k: [] for k in keys}
        for path in cls.network.paths:
            if len(path) == 1 and path[0] is cls.network.network_section:  # skip the pseudo equivalent network section
                pass
            else:
                d[keys[0]].append(repr(path))
                d[keys[1]].append(path.velocity_head(cls.units['pressure'], 3))
                d[keys[2]].append(path.elevation_head(cls.units['pressure'], 3))
                d[keys[3]].append(path.head_loss(cls.units['pressure'], 3))
                d[keys[4]].append(path.static_head(cls.units['pressure'], 3))
        return pd.DataFrame(d)

    @classmethod
    def set_zetas(cls, section_ids: List[str], zetas: List[float]):
        """
        Change the zeta value of the given sections.

        **Parameters:**

        - `section_ids`: (*List[str]*) list with the id's of the selected sections
        - `zetas`: (*List[float]*) list with the corresponding zeta values of the selected sections

        """
        selection = [cls.network.sections[selected_id] for selected_id in section_ids]
        for i, lst in enumerate(selection):
            for section in lst:
                section.zeta = zetas[i]
