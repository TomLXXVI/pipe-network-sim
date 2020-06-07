from typing import List, Tuple
import pandas as pd
import ipywidgets as iw
import matplotlib.figure as mpl_figure
import matplotlib.axes as mpl_axes
import matplotlib.pyplot as mpl_plt
import ipysheet
from lib.pypeflow.analysis import Analyzer
from lib.pypeflow.utils.system_curve import calculate_system_curves, draw_curves
from lib.pypeflow.utils.pump_curve import PumpCurve, calculate_pump_curve
from lib.pypeflow.core.resistance_coefficient import ResistanceCoefficient
from lib.pypeflow.core.pipe_schedules import PipeSchedule40
import lib.quantities as qty
import lib.nummath.graphing2 as graphing


class PipingNetworkSimulator:
    """
    Simulate water flow in a drinking water installation of an apartment building.
    """
    # widgets
    run_btn: iw.Button = None
    output: iw.Output = None
    fig: mpl_figure.Figure = None
    ax: mpl_axes.Axes = None
    table: ipysheet.Sheet = None
    sliders: List[iw.FloatSlider] = None
    slider_panel: iw.VBox = None
    dashboard: iw.VBox = None
    # data objects
    df_building: pd.DataFrame = None
    df_floors: pd.DataFrame = None
    graph: graphing.LineGraph = None
    pump_curve: PumpCurve = None
    sys_curves_axes: List[Tuple[List[float], List[float]]] = None
    pump_curve_axes: Tuple[List[float], List[float]] = None
    V_wp: qty.VolumeFlowRate = None
    dp_wp: qty.Pressure = None
    ratio: float = None
    total_flow_rate: float = None

    @classmethod
    def init(cls):
        cls._init_booster_pump()
        cls._init_analyzer('./input_files/hardy_gebouw.csv')
        cls.ratio = 100.0 / 7.245
        return cls._create_dashboard()

    @classmethod
    def _create_dashboard(cls):
        title = iw.HTML(value="<h1>Piping Network Simulator</h1>")
        cls._init_sliders()
        cls._init_run_btn()
        cls.output = iw.Output()
        container1 = iw.VBox([cls.run_btn, cls.output])
        container2 = iw.HBox([cls.slider_panel, container1])
        cls._init_table()
        cls._init_plot()
        cls.dashboard = iw.VBox([title, container2, cls.table, cls.fig.canvas])
        return cls.dashboard

    @classmethod
    def _init_sliders(cls):
        cls.sliders = []
        for i in range(8):
            slider_params = {
                'value': 6.0,
                'min': 1.0,
                'max': 100.0,
                'step': 1.0,
                'description': f'V{i + 1} [%]',
                'disabled': False,
                'continuous_update': False,
                'orientation': 'horizontal',
                'readout': True,
                'readout_format': '.1f'
            }
            cls.sliders.append(iw.FloatSlider(**slider_params))
        cls.slider_panel = iw.VBox(cls.sliders)

    @classmethod
    def _init_run_btn(cls):
        cls.run_btn = iw.Button(description='run')
        cls.run_btn.on_click(cls._run)

    @classmethod
    def _init_table(cls):
        cls.table = ipysheet.sheet(
            rows=cls.df_floors.shape[0],
            columns=cls.df_floors.shape[1],
            column_headers=cls.df_floors.columns.to_list()
        )
        cls.table.cells = [ipysheet.column(c, cls.df_floors[column].tolist()) for c, column in enumerate(cls.df_floors)]
        cls.table.layout.width = '1000px'

    @classmethod
    def _init_plot(cls):
        mpl_plt.ioff()
        cls.fig = mpl_plt.figure(figsize=(10, 6), dpi=96, tight_layout={'pad': 1})
        cls.ax = cls.fig.add_subplot(1, 1, 1)
        cls._setup_plot()

    @classmethod
    def _setup_plot(cls):
        cls.graph = draw_curves(
            sys_curves=cls.sys_curves_axes,
            pump_curves=[cls.pump_curve_axes],
            units=('L/s', 'bar'),
            V_max=3.0,
            V_step=0.1,
            p_max=6.0,
            p_step=0.5,
            working_point=(cls.V_wp('L/s'), cls.dp_wp('bar')),
            figure_constructs=(cls.fig, cls.ax)
        )
        cls.graph.draw(grid_on=True)
        mpl_plt.tight_layout()

    @classmethod
    def _run(cls, button):
        # event handler for run button click
        button.disabled = True
        # choose random resistance coefficients for the equivalent valves on each floor
        cls._print_message("Setting new Kv's")
        cls._choose_zetas()
        # re-analyze the piping network with the new resistance coefficients
        cls._print_message("Analyzing...")
        Analyzer.configure_network_from_df(cls.df_building, clear=True)
        try:
            cls._analyze()
        except (OverflowError, ValueError) as err:
            cls._print_message(err.args[0])
        else:
            # update dashboard
            cls._print_message("Updating...")
            cls._update_plot()
            cls._update_table()
            cls._print_message(f"Update finished. Total flow rate is {cls.total_flow_rate:.2f} L/s.")
        finally:
            button.disabled = False

    @classmethod
    def _update_table(cls):
        cls.table.cells = [ipysheet.column(c, cls.df_floors[column].tolist()) for c, column in enumerate(cls.df_floors)]

    @classmethod
    def _update_plot(cls):
        mpl_plt.cla()
        cls._setup_plot()
        cls.fig.canvas.draw()
        cls.fig.canvas.flush_events()

    @classmethod
    def _init_analyzer(cls, fp):
        Analyzer.set_units({
            'length': 'm',
            'diameter': 'mm',
            'flow_rate': 'L/s',
            'pressure': 'bar',
            'velocity': 'm/s'
        })
        Analyzer.create_network(
            start_node_id='n1',
            end_node_id='n0',
            fluid='water',
            fluid_temperature=10.0,
            pipe_schedule='pipe_schedule_40'
        )
        Analyzer.configure_network(fp)
        with open(fp) as fh:
            cls.df_building = pd.read_csv(fh)
        cls._analyze()

    @classmethod
    def _analyze(cls):
        # analyze piping network and get the new data for showing in the dashboard
        try:
            Analyzer.solve(error=1.0e-3, i_max=5000)
        except (OverflowError, ValueError) as err:
            raise err
        else:
            df = Analyzer.get_network()
            # only select the sections with equivalent "floor valves to show in the table"
            cls.df_floors = df.iloc[[r for r in range(1, 32, 4)], :]
            cls.total_flow_rate = cls.df_floors['flow_rate [L/s]'].sum()
            # calculate working point of booster pump
            cls.V_wp = Analyzer.network.flow_rate
            cls.dp_wp = cls.pump_curve.pump_pressure(cls.V_wp)
            # calculate flow path curves
            cls.sys_curves_axes = calculate_system_curves(
                flow_paths=Analyzer.network.paths,
                V_wp=Analyzer.network.flow_rate,
                dp_wp=cls.dp_wp,
                unit_flow_rate='L/s',
                unit_pressure='bar',
                V_end=qty.VolumeFlowRate(3.0, 'L/s')
            )

    @classmethod
    def _init_booster_pump(cls):
        cls.pump_curve = PumpCurve(desired_units={'flow_rate': 'L/s', 'pressure': 'bar'})
        cls.pump_curve.set_coefficients(
            (5.845551552456625, -0.58652795301776, -0.4598136345779855),
            {'flow_rate': 'L/s', 'pressure': 'bar'}
        )
        cls.pump_curve_axes = calculate_pump_curve(
            coeff_values=cls.pump_curve.get_coefficients(),
            coeff_units={'flow_rate': 'L/s', 'pressure': 'bar'},
            unit_flow_rate='L/s',
            unit_pressure='bar',
            V_end=qty.VolumeFlowRate(3.0, 'L/s')
        )

    @classmethod
    def _choose_zetas(cls):
        zeta1_pos = [i for i in range(1, 32, 4)]
        zeta2_pos = [i for i in range(7, 32, 4)]
        slider_values = [slider.value for slider in cls.sliders]
        Kv_values = [value / cls.ratio for value in slider_values]
        dn = qty.Length(40, 'mm')
        di = PipeSchedule40.inside_diameter(dn)
        zeta_values = [ResistanceCoefficient.from_Kv(Kv, di) for Kv in Kv_values]
        cls.df_building.loc[zeta1_pos, 'zeta'] = zeta_values
        cls.df_building.loc[zeta2_pos, 'zeta'] = zeta_values[:-1]

    @classmethod
    def _print_message(cls, message: str):
        with cls.output:
            print(message)
            cls.output.clear_output(wait=True)
