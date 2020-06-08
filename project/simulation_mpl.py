from typing import List, Tuple, Any
import pandas as pd
import ipywidgets as iw
import matplotlib.figure as mpl_figure
import matplotlib.axes as mpl_axes
import matplotlib.pyplot as mpl_plt
import ipysheet
from IPython.display import display, HTML
from lib.pypeflow.analysis import Analyzer
from lib.pypeflow.utils.system_curve import calculate_system_curves, draw_curves
from lib.pypeflow.utils.pump_curve import PumpCurve, calculate_pump_curve
from lib.pypeflow.core.resistance_coefficient import ResistanceCoefficient
from lib.pypeflow.core.pipe_schedules import PipeSchedule40
import lib.quantities as qty
import lib.nummath.graphing2 as graphing

FLOW_RATE_MAX = 3 * 1.162  # L/s (max. flow rate per storey @ 3 bar )
MAX_ITERATIONS = 10000
SINGLE_PUMP = (584555.1552456624, -58652795.301775984, -45981363457.798546)
PARALLEL_PUMPS = (584555.1552456621, -29326397.650887705, -11495340864.449678)
SERIES_PUMPS = (1169110.3104913242, -117305590.60355082, -91962726915.59743)


class PipingNetworkSimulator:
    """
    Simulate water flow in a drinking water installation of an apartment building.
    """
    # widgets
    run_btn: iw.Button = None
    output: iw.Output = None
    radio_buttons: iw.RadioButtons = None
    fig: mpl_figure.Figure = None
    ax: mpl_axes.Axes = None
    table: ipysheet.Sheet = None
    sliders: List[iw.FloatSlider] = None
    float_valve_pos: List[iw.FloatText] = None
    float_flow_rates: List[iw.FloatText] = None
    slider_panel: iw.HBox = None
    dashboard: iw.VBox = None
    # data objects
    df_building: pd.DataFrame = None
    df_floors: pd.DataFrame = None
    graph: graphing.LineGraph = None
    current_pump_curve: PumpCurve = None
    single_pump_curve: PumpCurve = None
    parallel_pump_curve: PumpCurve = None
    series_pump_curve: PumpCurve = None
    current_pump_curve_axes: Tuple[List[float], List[float]] = None
    single_pump_curve_axes: Tuple[List[float], List[float]] = None
    parallel_pump_curve_axes: Tuple[List[float], List[float]] = None
    series_pump_curve_axes: Tuple[List[float], List[float]] = None
    sys_curves_axes: List[Tuple[List[float], List[float]]] = None
    V_wp: qty.VolumeFlowRate = None
    dp_wp: qty.Pressure = None
    ratio: float = None
    total_flow_rate: float = None

    @classmethod
    def init(cls):
        cls._init_single_pump()
        cls._init_parallel_pumps()
        cls._init_series_pumps()
        cls._init_analyzer('./input_files/hardy_gebouw.csv')
        cls.ratio = 100.0 / 7.245
        return cls._create_dashboard()

    @classmethod
    def _create_dashboard(cls):
        title = iw.HTML(value="<h1>Piping Network Simulator</h1>")
        cls._init_sliders()
        cls._init_radio_buttons()
        cls._init_run_btn()
        cls.output = iw.Output()
        cls.output.layout.height='100px'
        container1 = iw.VBox([cls.radio_buttons, cls.run_btn, cls.output])
        container1.layout.justify_content = 'center'
        container2 = iw.HBox([cls.slider_panel, container1])
        cls._init_table()
        cls._init_plot()
        cls.dashboard = iw.VBox([title, container2, cls.table, cls.fig.canvas])
        cls.dashboard.layout.justify_content = 'space-around'
        return cls.dashboard

    @classmethod
    def _init_radio_buttons(cls):
        cls.radio_buttons = iw.RadioButtons(
            options=['single pump', '2 pumps in parallel', '2 pumps in series'],
            index=0,
            disabled=False
        )

    @classmethod
    def _init_sliders(cls):
        cls.sliders = []
        cls.float_valve_pos = []
        cls.float_flow_rates = []
        container_valves = [iw.HTML('<b>Valve opening position</b>')]
        container_flow_rates = [iw.HTML('<b>Expected flow rate</b>')]
        for i in range(8):
            slider_params = {
                'value': 6.0,
                'min': 1.0,
                'max': 100.0,
                'step': 1.0,
                'description': f'L{i + 1}',
                'disabled': False,
                'continuous_update': False,
                'orientation': 'horizontal',
                'readout': False,
                'readout_format': '.1f'
            }
            float_valve_pos_params = {
                'value': round(6.0, 3)
            }
            float_flow_rate_params = {
                'value': round(6.0 * FLOW_RATE_MAX / 100.0, 3),
                'disabled': True
            }
            fs = iw.FloatSlider(**slider_params)
            ft_vp = iw.FloatText(**float_valve_pos_params)
            ft_vp.layout.width = '50px'
            iw.link((fs, 'value'), (ft_vp, 'value'))
            ft_fr = iw.FloatText(**float_flow_rate_params)
            ft_fr.layout.width = '75px'
            cls.sliders.append(fs)
            cls.float_valve_pos.append(ft_vp)
            cls.float_flow_rates.append(ft_fr)
            h_box = iw.HBox([fs, ft_vp, iw.Label('%')])
            container_valves.append(h_box)
            container_flow_rates.append(iw.HBox([ft_fr, iw.Label('L/s')]))
        cls.slider_panel = iw.HBox([iw.VBox(container_valves), iw.VBox(container_flow_rates)])
        for slider, float_text in zip(cls.sliders, cls.float_flow_rates):
            slider.observe(cls._attach_event_handler(float_text), names='value')

    @classmethod
    def _attach_event_handler(cls, float_text):
        def handle_slider_change(change):
            nonlocal float_text
            float_text.value = round(change['new'] * FLOW_RATE_MAX / 100.0, 3)
        return handle_slider_change

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
        cls.table.cells = [ipysheet.column(c, cls.df_floors[column].tolist())
                           for c, column in enumerate(cls.df_floors)]
        # cls.table.layout.width = '1000px'

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
            pump_curves=[cls.single_pump_curve_axes, cls.parallel_pump_curve_axes, cls.series_pump_curve_axes],
            units=('L/s', 'bar'),
            V_max=6.0,
            V_step=0.5,
            p_max=10.0,
            p_step=0.5,
            working_point=(cls.V_wp('L/s'), cls.dp_wp('bar')),
            figure_constructs=(cls.fig, cls.ax),
            sys_curve_labels=('L1', 'eq. network', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
            pump_curve_labels=('single pump', '2 pumps in parallel', '2 pumps in series')
        )
        cls.graph.draw(grid_on=True)
        mpl_plt.tight_layout()

    @classmethod
    def _run(cls, button):
        # event handler for run button click
        button.disabled = True
        # choose random resistance coefficients for the equivalent valves on each floor
        cls._print_message(HTML('<p style="color:green;">Setting new Kv coefficients</p>'))
        cls._choose_zetas()
        cls._set_pump(cls.radio_buttons.index)
        # re-analyze the piping network with the new resistance coefficients
        cls._print_message(HTML('<p style="color:green;">Analyzing...</p>'))
        Analyzer.configure_network_from_df(cls.df_building, clear=True)
        try:
            cls._analyze()
        except OverflowError:
            cls._print_message(HTML(
                f'<p style="color:orange;">Could not find a solution after {MAX_ITERATIONS}.<br>'
                'Please try again with other values.</p>'
            ))
        except ValueError:
            cls._print_message(HTML(
                f'<p style="color:orange;">Could not solve the network with the current combination of valve openings.'
                f'<br>Please try again with other values.</p>'
            ))
        else:
            # update dashboard
            cls._print_message(HTML('<p style="color:green;">Updating...</p>'))
            cls._update_plot()
            cls._update_table()
            cls._print_message(HTML(
                f'<p style="color:green;">Update finished. Total flow rate is {cls.total_flow_rate:.2f} L/s.</p>'
            ))
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
        cls._set_pump(0)
        cls._analyze()

    @classmethod
    def _analyze(cls):
        # analyze piping network and get the new data for showing in the dashboard
        try:
            Analyzer.solve(error=1.0e-2, i_max=MAX_ITERATIONS)
        except (OverflowError, ValueError) as err:
            raise err
        else:
            df = Analyzer.get_network()
            # only select the sections with equivalent "floor valves to show in the table"
            cls.df_floors = df.iloc[[r for r in range(1, 32, 4)], :]
            cls.total_flow_rate = cls.df_floors['flow_rate [L/s]'].sum()
            # calculate working point of booster pump
            cls.V_wp = Analyzer.network.flow_rate
            cls.dp_wp = cls.current_pump_curve.pump_pressure(cls.V_wp)
            # calculate flow path curves
            cls.sys_curves_axes = calculate_system_curves(
                flow_paths=Analyzer.network.paths,
                V_wp=Analyzer.network.flow_rate,
                dp_wp=cls.dp_wp,
                unit_flow_rate='L/s',
                unit_pressure='bar',
                V_end=qty.VolumeFlowRate(6.0, 'L/s')
            )

    @classmethod
    def _init_single_pump(cls):
        cls.single_pump_curve = PumpCurve(desired_units={'flow_rate': 'L/s', 'pressure': 'bar'})
        cls.single_pump_curve.set_coefficients(
            (5.845551552456625, -0.58652795301776, -0.4598136345779855),
            {'flow_rate': 'L/s', 'pressure': 'bar'}
        )
        cls.single_pump_curve_axes = calculate_pump_curve(
            coeff_values=cls.single_pump_curve.get_coefficients(),
            coeff_units={'flow_rate': 'L/s', 'pressure': 'bar'},
            unit_flow_rate='L/s',
            unit_pressure='bar',
            V_end=qty.VolumeFlowRate(6.0, 'L/s')
        )

    @classmethod
    def _init_parallel_pumps(cls):
        cls.parallel_pump_curve = PumpCurve(desired_units={'flow_rate': 'L/s', 'pressure': 'bar'})
        cls.parallel_pump_curve.set_coefficients(
            (5.845551552456622, -0.2932639765088771, -0.11495340864449681),
            {'flow_rate': 'L/s', 'pressure': 'bar'}
        )
        cls.parallel_pump_curve_axes = calculate_pump_curve(
            coeff_values=cls.parallel_pump_curve.get_coefficients(),
            coeff_units={'flow_rate': 'L/s', 'pressure': 'bar'},
            unit_flow_rate='L/s',
            unit_pressure='bar',
            V_end=qty.VolumeFlowRate(6.0, 'L/s')
        )

    @classmethod
    def _init_series_pumps(cls):
        cls.series_pump_curve = PumpCurve(desired_units={'flow_rate': 'L/s', 'pressure': 'bar'})
        cls.series_pump_curve.set_coefficients(
            (11.691103104913244, -1.1730559060355084, -0.9196272691559745),
            {'flow_rate': 'L/s', 'pressure': 'bar'}
        )
        cls.series_pump_curve_axes = calculate_pump_curve(
            coeff_values=cls.series_pump_curve.get_coefficients(),
            coeff_units={'flow_rate': 'L/s', 'pressure': 'bar'},
            unit_flow_rate='L/s',
            unit_pressure='bar',
            V_end=qty.VolumeFlowRate(3.0, 'L/s')
        )

    @classmethod
    def _set_single_pump(cls):
        cls.current_pump_curve = cls.single_pump_curve
        cls.current_pump_curve_axes = cls.single_pump_curve_axes

    @classmethod
    def _set_parallel_pumps(cls):
        cls.current_pump_curve = cls.parallel_pump_curve
        cls.current_pump_curve_axes = cls.parallel_pump_curve_axes

    @classmethod
    def _set_series_pumps(cls):
        cls.current_pump_curve = cls.series_pump_curve
        cls.current_pump_curve_axes = cls.single_pump_curve_axes

    @classmethod
    def _set_pump(cls, index):
        if index == 0:
            cls.df_building.loc[0, ('a0', 'a1', 'a2')] = SINGLE_PUMP
            cls._set_single_pump()
        elif index == 1:
            cls.df_building.loc[0, ('a0', 'a1', 'a2')] = PARALLEL_PUMPS
            cls._set_parallel_pumps()
        else:
            cls.df_building.loc[0, ('a0', 'a1', 'a2')] = SERIES_PUMPS
            cls._set_series_pumps()

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
    def _print_message(cls, message: HTML):
        with cls.output:
            display(message)
            cls.output.clear_output(wait=True)
