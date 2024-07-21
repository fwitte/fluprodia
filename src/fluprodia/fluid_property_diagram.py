# -*- coding: utf-8

"""Module for fluid property diagram creation.

This file is part of project fluprodia (github.com/fwitte/fluprodia). It's
copyrighted by the contributors recorded in the version control history of the
file, available from its original location
src/fluprodia/fluid_property_diagram.py

SPDX-License-Identifier: MIT
"""
import json
import os

import CoolProp as CP
import numpy as np

from fluprodia import __version__


def _beautiful_unit_string(unit):
    r"""Convert unit fractions to latex.

    Parameters
    ----------
    unit : str
        Value of unit for input, e.g. :code:`m^3/kg`.

    Returns
    -------
    unit : str
        Value of unit for output, e.g. :code:`$\frac{m^3}{kg}$`.
    """
    if '/' in unit:
        numerator = unit.split('/')[0]
        denominator = unit.split('/')[1]
        unit = '$\\frac{' + numerator + '}{' + denominator + '}$'

    return unit


def _isolines_log(val_min, val_max):
    """Generate default logarithmic isolines.

    Parameters
    ----------
    val_min : float
        Minimum value for isoline range.

    val_max : float
        Maximum value for isoline range.

    Returns
    -------
    arr : ndarray
        numpy array with logarithmically spaced values starting from the
        minimum value going to the maximum value in steps of :code:`1ek`,
        :code:`2ek` and :code:`5ek`.
    """
    arr = [val_min]
    digits = int(np.floor(np.log10(val_min)))
    while arr[-1] < val_max:
        arr += [1 * 10 ** digits]
        arr += [2 * 10 ** digits]
        arr += [5 * 10 ** digits]
        digits += 1

    arr = np.unique(np.asarray(arr + [val_max]))
    return arr[(arr >= val_min) & (arr <= val_max)]


class FluidPropertyDiagram:
    u"""Short summary.

    Parameters
    ----------
    fluid : str
        Fluid for diagram.

    width : float
        Width of all diagrams (default value: :code:`width=16.0`).

    height : float
        Height of all diagrams (default value: :code:`height=10.0`).

    Example
    -------
    This is a small example of how to create a fluid property dataset for water
    and export to Ts-, hs- and logph-diagram.

    >>> from fluprodia import FluidPropertyDiagram
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> diagram = FluidPropertyDiagram('water')

    After object creation it is possible to specify isolines. There are deault
    isolines available, but these might not suit your requirements. We will
    define temperature and enthalpy isolines for this case. Before that, we
    need to specify the units if we do not want to use SI units. For available
    units see the
    :py:meth:`fluprodia.fluid_property_diagram.FluidPropertyDiagram.set_unit_system`
    documentation.

    >>> diagram.set_unit_system(T='°C', p='MPa', s='kJ/kgK', h='kJ/kg')
    >>> iso_T = np.arange(50, 701, 50)
    >>> iso_h = np.arange(0, 3601, 200)
    >>> diagram.set_isolines(T=iso_T, h=iso_h)

    Now we can calculate the diagram data and create the diagram, the data
    should be plotted to. If you want to plot additional data to your diagram,
    first call the `draw_isolines` method, and then plot your data.

    >>> diagram.calc_isolines()

    After that it is possible to specify the view range of the plot and draw
    the isolines of a specific type of diagram. Last step is to export your
    diagram. Any file format supported by matplotlib is possible.

    >>> fig, ax = plt.subplots(1)
    >>> diagram.draw_isolines(diagram_type='Ts', fig=fig, ax=ax, x_min=0, x_max=8, y_min=0, y_max=700)
    >>> plt.tight_layout()
    >>> fig.savefig('Ts_Diagramm.pdf')

    If we want to create a different diagram, e.g. hs-diagram, it is not
    necessary to recalculate the isolines. Instead, create a new figure and draw
    the isolines for a different diagram.

    >>> fig, ax = plt.subplots(1, figsize=(8, 5))
    >>> diagram.draw_isolines(diagram_type='hs', fig=fig, ax=ax, x_min=0, x_max=8, y_min=0, y_max=3600)
    >>> plt.tight_layout()
    >>> fig.savefig('hs_Diagramm.pdf')

    >>> fig, ax = plt.subplots(1, figsize=(8, 5))
    >>> diagram.draw_isolines(diagram_type='logph', fig=fig, ax=ax, x_min=0, x_max=3600, y_min=1e-2, y_max=5e2)
    >>> plt.tight_layout()
    >>> fig.savefig('logph_Diagramm.pdf')

    It is also possible to specify/modify the isolines to plot. For example,
    the lines of constant specific enthalpy should be plotted in red color
    and with a linewidth of 2. Also, the lines of constant specific volumen
    should not be plotted at all. For more information see the
    :py:meth:`fluprodia.fluid_property_diagram.FluidPropertyDiagram.draw_isolines`
    method.

    >>> fig, ax = plt.subplots(1)
    >>> diagram.draw_isolines(
    ...     diagram_type='Ts', fig=fig, ax=ax,
    ...     isoline_data={
    ...         'h': {
    ...            'values': iso_h,
    ...            'style': {'linewidth': 2, 'color': '#ff0000'}
    ...         }, 'v': {'values': np.array([])}},
    ...     x_min=0, x_max=8, y_min=0, y_max=700
    ... )
    >>> plt.tight_layout()
    >>> fig.savefig('Ts_Diagramm.pdf')

    It is also possible to export the data of a diagram:

    >>> diagram.to_json("./tmp/water.json")

    And, if you have your data ready to use, you can import the data again like
    this:

    >>> pressure_isolines_before_export = diagram.pressure["isolines"]
    >>> pressure_unit_before_export = diagram.units["p"]
    >>> diagram = FluidPropertyDiagram.from_json("tmp/water.json")
    >>> diagram.fluid
    'water'
    >>> diagram.units["p"] == pressure_unit_before_export
    True
    >>> all(diagram.pressure["isolines"] == pressure_isolines_before_export)
    True

    """

    def __init__(self, fluid):
        u"""Create a FluidPropertyDiagram object.

        Parameters
        ----------
        fluid : str
            Fluid for diagram.

        width : float
            Width of all diagrams (default value: :code:`width=16.0`).

        height : float
            Height of all diagrams (default value: :code:`height=10.0`).
        """
        self.fluid = fluid
        self.state = CP.AbstractState('HEOS', self.fluid)

        self.converters = {}
        self.converters['p'] = {
            'Pa': 1, 'hPa': 1e2, 'mbar': 1e2, 'psi': 6894.7572931783,
            'kPa': 1e3, 'bar': 1e5, 'MPa': 1e6
        }
        self.converters['T'] = {
            'K': [0, 1], '°C': [273.15, 1], '°F': [459.67, 5 / 9]
        }
        self.converters['s'] = {'J/kgK': 1, 'kJ/kgK': 1e3, 'MJ/kgK': 1e6}
        self.converters['h'] = {'J/kg': 1, 'kJ/kg': 1e3, 'MJ/kg': 1e6}
        self.converters['v'] = {'m^3/kg': 1, 'l/kg': 1e-3}
        self.converters['Q'] = {'-': 1, '%': 0.01}

        self.units = {
            'p': 'Pa',
            's': 'J/kgK',
            'h': 'J/kg',
            'v': 'm^3/kg',
            'Q': '-',
            'T': 'K',
        }

        self.properties = {
            'p': 'pressure',
            'v': 'volume',
            'T': 'temperature',
            'h': 'enthalpy',
            's': 'entropy',
            'Q': 'quality'
        }

        self.supported_diagrams = {
            'Ts': {
                'x_property': 's',
                'y_property': 'T',
                'x_scale': 'linear',
                'y_scale': 'linear'
            },
            'hs': {
                'x_property': 's',
                'y_property': 'h',
                'x_scale': 'linear',
                'y_scale': 'linear'
            },
            'logph': {
                'x_property': 'h',
                'y_property': 'p',
                'x_scale': 'linear',
                'y_scale': 'log'
            },
            'Th': {
                'x_property': 'h',
                'y_property': 'T',
                'x_scale': 'linear',
                'y_scale': 'linear'
            },
            'plogv': {
                'x_property': 'v',
                'y_property': 'p',
                'x_scale': 'log',
                'y_scale': 'linear'
            }
        }

        self._setup_functions_and_inputs()
        self._setup_datastructures()
        self._setup_isoline_defaults()
        self._setup_default_line_layout()
        self._setup_default_label_positioning()

    @classmethod
    def from_json(cls, path):
        with open(path, "r") as f:
            data = json.load(f)

        metadata = data["META"].copy()
        del data["META"]
        instance = cls(metadata["fluid"])
        instance.set_unit_system(**metadata["units"])
        for key in data:
            isoprop = getattr(instance, key)
            isoprop["isolines"] = np.array([float(value) for value in data[key].keys()])
            for isoline in data[key]:
                datapoints = {}
                for prop, values in data[key][isoline].items():
                    datapoints[prop] = np.array(values)

                isoprop[float(isoline)] = datapoints

        return instance

    def _setup_functions_and_inputs(self):
        """Setup lookup tables for isoline functions and CoolProp inputs."""
        self.single_isoline_functions = {
            'p': self._single_isobaric,
            'v': self._single_isochoric,
            'T': self._single_isothermal,
            'h': self._single_isenthalpic,
            's': self._single_isentropic
        }
        self.CoolProp_inputs = {
            'p': CP.iP,
            'v': CP.iDmass,
            'T': CP.iT,
            'h': CP.iHmass,
            's': CP.iSmass,
            'Q': CP.iQ
        }
        self.CoolProp_results = {
            'p': self.state.p,
            'v': self.state.rhomass,
            'T': self.state.T,
            'h': self.state.hmass,
            's': self.state.smass,
            'Q': self.state.Q
        }

    def _setup_datastructures(self):
        """Set up datastructures for all isolines."""
        self.pressure = {'isolines': np.array([])}
        self.entropy = {'isolines': np.array([])}
        self.temperature = {'isolines': np.array([])}
        self.enthalpy = {'isolines': np.array([])}
        self.volume = {'isolines': np.array([])}
        self.quality = {'isolines': np.array([])}

    def _setup_default_line_layout(self):
        """Definition of the default isoline layout."""
        self.pressure['style'] = {
            'linestyle': '-.',
            'color': '#363636',
            'linewidth': 0.5
        }
        self.volume['style'] = {
            'linestyle': 'dotted',
            'color': '#363636',
            'linewidth': 0.5
        }
        self.quality['style'] = {
            'linestyle': '-',
            'color': '#363636',
            'linewidth': 0.5
        }
        self.entropy['style'] = {
            'linestyle': 'solid',
            'color': '#d1d1d1',
            'linewidth': 0.5
        }
        self.enthalpy['style'] = {
            'linestyle': '--',
            'color': '#363636',
            'linewidth': 0.5
        }
        self.temperature['style'] = {
            'linestyle': '--',
            'color': '#363636',
            'linewidth': 0.5
        }

    def _setup_default_label_positioning(self):
        """Definition of the default label positioning."""
        self.pressure['label_position'] = 0.85
        self.volume['label_position'] = 0.7
        self.quality['label_position'] = 0.225
        self.entropy['label_position'] = 0.95
        self.enthalpy['label_position'] = 0.85
        self.temperature['label_position'] = 0.95

    def set_isolines(self, **kwargs):
        """Set the isolines.

        Parameters
        ----------
        p : ndarray
            Isolines for pressure.

        T : ndarray
            Isolines for temperature.

        Q : ndarray
            Isolines for vapor mass fraction.

        s : ndarray
            Isolines for specific entropy.

        h : ndarray
            Isolines for specific enthalpy.

        v : ndarray
            Isolines for specific volume.
        """
        keys = ['p', 'T', 'Q', 's', 'h', 'v']
        for key in kwargs:
            if key in keys:
                obj = getattr(self, self.properties[key])
                obj['isolines'] = self.convert_to_SI(kwargs[key], key).round(8)
            else:
                msg = (
                    f'The specified isoline \'{key}\' is not available. '
                    f'Select from: {", ".join(keys)}.'
                )
                print(msg)

    def _setup_isoline_defaults(self):
        """Calculate the default values for the isolines."""
        self.p_trip = self.state.trivial_keyed_output(CP.iP_triple)
        self.p_max = self.state.trivial_keyed_output(CP.iP_max)
        self.T_trip = self.state.trivial_keyed_output(CP.iT_triple)
        self.T_max = self.state.trivial_keyed_output(CP.iT_max)

        self.p_crit = self.state.trivial_keyed_output(CP.iP_critical)
        self.T_crit = self.state.trivial_keyed_output(CP.iT_critical)
        self.v_crit = 1 / self.state.trivial_keyed_output(CP.irhomass_critical)

        self.state.update(CP.PQ_INPUTS, (self.p_crit + self.p_trip) / 2, 1)
        self.v_intermediate = 1 / self.state.rhomass()

        self.state.update(CP.PT_INPUTS, self.p_trip, self.T_max)
        self.v_max = 1 / self.state.rhomass()
        self.s_max = self.state.smass()
        self.h_max = self.state.hmass()

        self.state.update(CP.QT_INPUTS, 0, self.T_trip + 1)
        self.s_min = self.state.smass()
        self.h_min = self.state.hmass()

        self.state.specify_phase(CP.iphase_liquid)
        p = self.p_crit
        while True:
            try:
                self.state.update(CP.PT_INPUTS, p, self.T_trip + 1)
                break
            except ValueError:
                p *= 0.999
        self.v_min = 1 / self.state.rhomass()
        self.state.unspecify_phase()

        self.quality['isolines'] = np.linspace(0, 1, 11).round(8)

        step = round(int(self.T_max - self.T_trip) / 15, -1)
        self.temperature['isolines'] = np.append(
            self.T_trip,
            np.arange(self.T_max, self.T_trip, -step)[::-1]).round(8)

        step = round(int(self.s_max - self.s_min) / 15, -1)
        self.entropy['isolines'] = np.arange(self.s_min, self.s_max, step).round(8)

        step = round(int(self.h_max - self.h_min) / 15, -1)
        self.enthalpy['isolines'] = np.arange(0, self.h_max, step).round(8)

        self.pressure['isolines'] = _isolines_log(
            self.p_trip + 1e-2, self.p_max).round(8)
        self.volume['isolines'] = _isolines_log(self.v_min, self.v_max).round(8)

    def set_isolines_from_pT(self):
        pass

    def _isolines_from_pT_boundaries(self):
        pass

    def _get_state_result_by_name(self, property_name):
        if property_name == "v":
            return 1 / self.CoolProp_results[property_name]()
        else:
            return self.CoolProp_results[property_name]()

    def _update_state(self, input_dict):
        output_dict = {}
        for property_name, value in input_dict.items():
            if property_name == "d":
                output_dict["d"] = 1 / value
            else:
                output_dict[property_name] = value

        args = [
            arg for property_name, value in output_dict.items()
            for arg in [self.CoolProp_inputs[property_name], value]
        ]

        self.state.update(*CP.CoolProp.generate_update_pair(*args))

    def set_unit_system(self, **kwargs):
        u"""Set the unit system for the fluid properties.

        Parameters
        ----------
        p : str
            Unit of pressure, units available are
            :code:`Pa, hPa, mbar, psi, kPa, bar, MPa`.

        T : str
            Unit of temperatur, units available are
            :code:`K, °C, °F`.

        s : str
            Unit of specific entropy, units available are
            :code:`J/kgK, kJ/kgK, MJ/kgK`.

        h : str
            Unit of specific enthalpy, units available are
            :code:`J/kg, kJ/kg, MJ/kg`.

        v : str
            Unit of specific volume, units available are
            :code:`m^3/kg, l/kg`.

        Q : str
            Unit of vapor mass fraction, units available are
            :code:`-, %`.
        """
        for key, value in kwargs.items():
            if value in self.converters[key].keys():
                self.units[key] = value
            else:
                msg = (
                    f'The unit {value} is not available for the  fluid '
                    f'property {self.properties[key]} . Please choose from '
                    f'{", ".join(self.converters[key].keys())}.'
                )
                raise ValueError(msg)

    def _draw_isoline_label(self, fig, ax, isoline, property, idx, x, y, x_min, x_max, y_min, y_max):
        """Draw a label for an isoline.

        Parameters
        ----------
        isoline : float
            Value of the isoline.

        property : str
            Fluid property of the isoline.

        idx : float
            Index in the array holding the isoline data, where the label
            should be plotted.

        x : ndarray
            x-values of the isoline.s

        y : ndarray
            y-values of the isoline.s
        """
        fig_size = fig.get_size_inches()
        pos = ax.get_position()
        width, height = (pos.width * fig_size[0], pos.height * fig_size[1])

        if (idx > len(x) or
                idx == 0 or
                x[idx] > x_max or x[idx] < x_min or
                y[idx] > y_max or y[idx] < y_min or
                x[idx - 1] > x_max or x[idx - 1] < x_min or
                y[idx - 1] > y_max or y[idx - 1] < y_min):
            return

        idx -= 1

        if x[idx] - x[idx - 1] == 0:
            if y[idx] > y[idx - 1]:
                alpha = 90
            else:
                alpha = -90
        elif y[idx] - y[idx - 1] == 0:
            alpha = 0
        else:
            if ax.get_xscale() == 'log':
                x_scaled = (
                    (np.log(x[idx]) - np.log(x[idx - 1]))
                    * (width / (np.log(x_max) - np.log(x_min)))
                )
            else:
                x_scaled = (x[idx] - x[idx - 1]) * (width / (x_max - x_min))

            if ax.get_yscale() == 'log':
                y_scaled = (
                    (np.log(y[idx]) - np.log(y[idx - 1]))
                    * (height / (np.log(y_max) - np.log(y_min)))
                )
            else:
                y_scaled = (y[idx] - y[idx - 1]) * (height / (y_max - y_min))

            alpha = np.arctan(y_scaled / x_scaled) / (2 * np.pi) * 360

        unit = _beautiful_unit_string(self.units[property])

        txt = f'{isoline} {unit}'
        text_bg_color = ax.get_facecolor()
        ax.text(
            x[idx], y[idx], txt, fontsize=5,
            rotation=alpha, va='center', ha='center',
            bbox=dict(facecolor=text_bg_color, edgecolor=text_bg_color, pad=0.0)
        )

    def calc_isolines(self):
        """Calculate all isolines."""
        self._isobaric()
        self._isochoric()
        self._isoquality()
        self._isenthalpic()
        self._isothermal()
        self._isentropic()

    def _isobaric(self):
        """Calculate an isoline of constant pressure."""
        isolines = self.pressure['isolines']

        iterator = 1 / np.append(
            np.geomspace(self.v_min, self.v_intermediate, 100, endpoint=False),
            np.geomspace(self.v_intermediate, self.v_max, 100))

        for p in isolines.round(8):
            if p <= self.p_crit:
                T_sat = CP.CoolProp.PropsSI("T", "P", p, "Q", 0, self.fluid)
                iterators = [
                    ("T", np.geomspace(self.T_trip, T_sat * 0.999, 120)),
                    # start in liquid and end in gas
                    ("Q", np.linspace(0, 1, 11)),
                    ("T", np.geomspace(T_sat * 1.001, self.T_max, 69))
                ]
            else:
                iterators = [
                    ("v", np.geomspace(self.v_min, self.v_intermediate, 100, endpoint=False)),
                    ("v", np.geomspace(self.v_intermediate, self.v_max, 100))
                ]
            self.pressure[p] = self._single_isoline(iterators, "p", p)

    def _isochoric(self):
        """Calculate an isoline of constant specific volume."""
        isolines = self.volume['isolines']

        iterator = np.append(
            np.geomspace(self.p_trip, self.p_crit * 0.8, 100, endpoint=False),
            np.geomspace(self.p_crit * 0.8, self.p_max, 100)
        )

        for v in isolines.round(8):
            self.volume[v] = self._single_isochoric(
                iterator, np.ones(len(iterator)) * 1 / v)

    def _isothermal(self):
        """Calculate an isoline of constant temperature."""
        isolines = self.temperature['isolines']

        for T in isolines.round(8):
            if T <= self.T_crit:
                self.state.update(CP.QT_INPUTS, 0, T)
                p_sat = self.state.p()
                iterators = [
                    ("p", np.geomspace(self.p_trip, p_sat * 0.999, 120)),
                    # start in gas and end in liquid
                    ("Q", np.linspace(1, 0, 11)),
                    ("p", np.geomspace(p_sat * 1.001, self.p_max, 69))
                ]
            elif T <= self.T_crit * 1.2:
                self.state.update(CP.PT_INPUTS, self.p_crit * 0.7, T)
                s_start = self.state.smass()
                self.state.update(CP.PT_INPUTS, self.p_crit * 1.2, T)
                s_end = self.state.smass()
                iterators = [
                    ("p", np.geomspace(self.p_trip, self.p_crit * 0.7, 80, endpoint=False)),
                    ("s", np.linspace(s_start, s_end, 40, endpoint=False)),
                    ("p", np.geomspace(self.p_crit * 1.2, self.p_max, 80)),
                ]
            else:
                iterators = [
                    ("p", np.geomspace(self.p_trip, self.p_max, 200)),
                ]

            self.temperature[T] = self._single_isoline(iterators, "T", T)

    def _isoquality(self):
        """Calculate an isoline of constant vapor mass fraction."""
        isolines = self.quality['isolines']

        iterators = [
            ("T", np.append(
                np.linspace(self.T_trip, self.T_crit * 0.97, 40, endpoint=False),
                np.linspace(self.T_crit * 0.97, self.T_crit, 40)
            ))
        ]

        for Q in isolines.round(8):
            self.quality[Q] = self._single_isoline(iterators, "Q", Q)

    def _isenthalpic(self):
        """Calculate an isoline of constant specific enthalpy."""
        isolines = self.enthalpy['isolines']

        iterators = [
            ("v", np.geomspace(self.v_min, self.v_intermediate, 100, endpoint=False)),
            ("v", np.geomspace(self.v_intermediate, self.v_max, 100))
        ]

        for h in isolines.round(8):
            self.enthalpy[h] = self._single_isoline(iterators, "h", h)

    def _isentropic(self):
        """Calculate an isoline of constant specific entropy."""
        isolines = self.entropy['isolines']

        iterators = [
            ('p', np.append(
                np.geomspace(self.p_trip, self.p_crit * 0.8, 100, endpoint=False),
                np.geomspace(self.p_crit * 0.8, self.p_max, 100)
            ))
        ]

        for s in isolines.round(8):
            self.entropy[s] = self._single_isoline(iterators, 's', s)

    def to_json(self, path):
        """Export the diagram data as json file to a path.

        Parameters
        ----------
        path : str, path-like
            Name of the file to export the data to.
        """
        data = {
            prop: {
                f"{isoline}": {
                    subprop: list(getattr(self, prop)[isoline][subprop].astype(float))
                    for subprop in getattr(self, prop)[isoline]
                } for isoline in getattr(self, prop)["isolines"]
            } for prop in self.properties.values()
        }

        data["META"] = {
            "fluid": self.fluid,
            "units": self.units,
            "CoolProp-version": CP.__version__,
            "fluprodia-version": __version__
        }

        directory = os.path.dirname(path)
        if directory != "" and not os.path.exists(directory):
            os.makedirs(directory)

        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

    def calc_individual_isoline(
            self, isoline_property=None,
            isoline_value=None,
            isoline_value_end=None,
            starting_point_property=None, ending_point_property=None,
            starting_point_value=None, ending_point_value=None):
        """Return data points of an individual isoline within.

        Pass the isoline type, its value in the diagrams unit system, as well
        as the start and the endpoint of the line. Styling can be changed using
        the line_style property.

        Parameters
        ----------
        isoline_property : str
            Type of the isoline. Choose from :code:`line_type='...'`:

            - pressure (:code:`'p'`)
            - specific volume (:code:`'v'`)
            - temperature (:code:`'T'`)
            - enthalpy (:code:`'h'`)
            - entropy (:code:`'s'`)

        isoline_value : float
            Value of the isoline specified in the respective unit.

        starting_point : dict
            Dictionary holding the starting property and its value in the
            unit used for plotting the diagram, e.g.
            :code:`starting_property{'p': 1e5}`

        ending_point : dict
            Dictionary holding the ending property and its value in the
            unit used for plotting the diagram, e.g.
            :code:`ending_property{'p': 1e4}`

        Returns
        -------
        datapoints : dict
            Dictionary holding the isoline datapoints for:

            - pressure (key=:code:`'p'`)
            - specific volume (key=:code:`'v'`)
            - temperature (key=:code:`'T'`)
            - enthalpy (key=:code:`'h'`)
            - entropy (key=:code:`'s'`)

        Example
        -------
        A full example can be found in the class documentation.
        """
        if not isinstance(isoline_property, str):
            msg = 'Parameter isoline_property must be specified as string!'
            raise ValueError(msg)
        elif (isoline_property not in self.properties.keys() or
              isoline_property == 'Q'):
            msg = 'Isoline of type ' + isoline_property + ' not available.'
            raise ValueError(msg)

        f = self.single_isoline_functions[isoline_property]

        isoline_value = self.convert_to_SI(isoline_value, isoline_property)

        if isoline_value_end is None:
            isoline_value_end = isoline_value
        else:
            isoline_value_end = self.convert_to_SI(
                isoline_value_end, isoline_property
            )

        starting_point_value = self.convert_to_SI(
            starting_point_value, starting_point_property)
        ending_point_value = self.convert_to_SI(
            ending_point_value, ending_point_property)

        if isoline_property == 'v':
            isoline_value = 1 / isoline_value
            isoline_value_end = 1 / isoline_value_end
            isoline_property = 'D'
        if starting_point_property == 'v':
            starting_point_value = 1 / starting_point_value
            starting_point_property = 'D'
        if ending_point_property == 'v':
            ending_point_value = 1 / ending_point_value
            ending_point_property = 'D'

        if isoline_property == 'D':

            if starting_point_property == 'p':
                pressure_start = starting_point_value
            else:
                pressure_start = CP.CoolProp.PropsSI(
                    'P', starting_point_property.upper(),
                    starting_point_value, 'D', isoline_value, self.fluid
                )

            if ending_point_property == 'p':
                pressure_end = ending_point_value
            else:
                pressure_end = CP.CoolProp.PropsSI(
                    'P', ending_point_property.upper(),
                    ending_point_value, 'D', isoline_value_end, self.fluid
                )

            iterator = np.geomspace(pressure_start, pressure_end, 100)

        elif isoline_property == 'T':
            if starting_point_property == 's':
                entropy_start = starting_point_value
            else:
                entropy_start = CP.CoolProp.PropsSI(
                    'S', starting_point_property.upper(),
                    starting_point_value, isoline_property.upper(),
                    isoline_value, self.fluid
                )

            if ending_point_property == 's':
                entropy_end = ending_point_value
            else:
                entropy_end = CP.CoolProp.PropsSI(
                    'S', ending_point_property.upper(),
                    ending_point_value, isoline_property.upper(),
                    isoline_value_end, self.fluid
                )

            iterator = np.linspace(entropy_start, entropy_end, 100)

        else:
            if starting_point_property == 'D':
                density_start = starting_point_value
            else:
                density_start = CP.CoolProp.PropsSI(
                    'D', starting_point_property.upper(),
                    starting_point_value, isoline_property.upper(),
                    isoline_value, self.fluid
                )

            if ending_point_property == 'D':
                density_end = ending_point_value
            else:
                density_end = CP.CoolProp.PropsSI(
                    'D', ending_point_property.upper(),
                    ending_point_value, isoline_property.upper(),
                    isoline_value_end, self.fluid
                )

            iterator = np.geomspace(density_start, density_end, 100)

        if isoline_property == 'p':
            density_change = np.append(
                0, np.diff(iterator)[::-1] / (iterator[0] - iterator[-1]))
            isoline_vector = isoline_value + density_change.cumsum() * (
                isoline_value - isoline_value_end)
        else:
            isoline_vector = np.linspace(isoline_value, isoline_value_end, 100)

        datapoints = f(iterator, isoline_vector)

        for key in datapoints.keys():
            datapoints[key] = self.convert_from_SI(datapoints[key], key)

        return datapoints

    def _single_isobaric(self, iterator, p):
        """Calculate an isoline of constant pressure."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': [], 'Q': []}

        if iterator[0] < iterator[-1]:
            rising = True
        else:
            rising = False

        for i, val in enumerate(iterator):
            try:
                self.state.update(CP.DmassP_INPUTS, val, p[i])
                datapoints['h'] += [self.state.hmass()]
                datapoints['T'] += [self.state.T()]
                datapoints['v'] += [1 / val]
                datapoints['s'] += [self.state.smass()]
                datapoints['p'] += [p[i]]
                if p[i] < self.p_crit:
                    datapoints['Q'] += [self.state.Q()]
                else:
                    datapoints['Q'] += [-1]

            except ValueError as e:
                pass

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        data = self._get_Q_crossings(datapoints, 'p', rising)
        if len(data) == 0:
            return datapoints
        else:
            return self._insert_Q_crossings(datapoints, 'p', data)

    def _single_isochoric(self, iterator, D):
        """Calculate an isoline of constant specific volume."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': [], 'Q': []}

        if iterator[0] < iterator[-1]:
            rising = True
        else:
            rising = False

        for i, val in enumerate(iterator):
            try:
                self.state.update(CP.DmassP_INPUTS, D[i], val)
                datapoints['h'] += [self.state.hmass()]
                datapoints['T'] += [self.state.T()]
                datapoints['v'] += [1 / D[i]]
                datapoints['s'] += [self.state.smass()]
                datapoints['p'] += [val]
                if val < self.p_crit:
                    datapoints['Q'] += [self.state.Q()]
                else:
                    datapoints['Q'] += [-1]
            except ValueError:
                pass

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        data = self._get_Q_crossings(datapoints, 'v', rising)
        if len(data) == 0:
            return datapoints
        else:
            return self._insert_Q_crossings(datapoints, 'v', data)

    def _single_isothermal(self, iterator, T):
        """Calculate an isoline of constant temperature."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': [], 'Q': []}

        if iterator[0] < iterator[-1]:
            rising = True
        else:
            rising = False

        for i, val in enumerate(iterator):
            try:
                self.state.update(CP.SmassT_INPUTS, val, T[i])
                datapoints['T'] += [T[i]]
                datapoints['p'] += [self.state.p()]
                datapoints['v'] += [1 / self.state.rhomass()]
                datapoints['s'] += [val]
                datapoints['h'] += [self.state.hmass()]
                if T[i] < self.T_crit:
                    datapoints['Q'] += [self.state.Q()]
                else:
                    datapoints['Q'] += [-1]
            except ValueError:
                # for some reason PropSI inputs are way more stable here
                try:
                    p = CP.CoolProp.PropsSI(
                        'P', 'T', T[i], 'S', val, self.fluid
                    )
                    self.state.update(CP.PSmass_INPUTS, val, p)
                    datapoints['T'] += [T[i]]
                    datapoints['p'] += [p]
                    datapoints['v'] += [1 / self.state.rhomass()]
                    datapoints['s'] += [val]
                    datapoints['h'] += [self.state.hmass()]
                except ValueError:
                    pass

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        data = self._get_Q_crossings(datapoints, 'T', rising)
        if len(data) == 0:
            return datapoints
        else:
            return self._insert_Q_crossings(datapoints, 'T', data)

    def _single_isoline(self, iterators, isoline_property, isoline_value):
        """Calculate an isoline of constant temperature."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': [], 'Q': []}

        num_points = sum([len(_[1]) for _ in iterators])
        datapoints[isoline_property] = np.ones(num_points) * isoline_value
        for iterator_property, iterator_values in iterators:
            self.state.unspecify_phase()
            datapoints[iterator_property] += iterator_values.tolist()

            for value in iterator_values:
                try:
                    self._update_state({isoline_property: isoline_value, iterator_property: value})
                    success = True
                except ValueError as e:
                    success = False

                result = np.nan
                for result_property in set(datapoints.keys()) - {iterator_property, isoline_property}:
                    if success:
                        result = self._get_state_result_by_name(result_property)

                    datapoints[result_property] += [result]

        for key in set(datapoints.keys()) - {isoline_property}:
            datapoints[key] = np.asarray(datapoints[key])

        return datapoints

    def _single_isenthalpic(self, iterator, h):
        """Calculate an isoline of constant specific enthalpy."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': []}

        for i, val in enumerate(iterator):
            try:
                self.state.update(CP.DmassHmass_INPUTS, val, h[i])
                datapoints['T'] += [self.state.T()]
                datapoints['p'] += [self.state.p()]
                datapoints['v'] += [1 / val]
                datapoints['s'] += [self.state.smass()]
                datapoints['h'] += [h[i]]
            except ValueError:
                pass

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        return datapoints

    def _single_isentropic(self, iterator, s):
        """Calculate an isoline of constant specific entropy."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': []}

        for i, val in enumerate(iterator):
            try:
                self.state.update(CP.DmassSmass_INPUTS, val, s[i])
                datapoints['T'] += [self.state.T()]
                datapoints['p'] += [self.state.p()]
                datapoints['v'] += [1 / val]
                datapoints['s'] += [s[i]]
                datapoints['h'] += [self.state.hmass()]
            except ValueError:
                pass

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        return datapoints

    def _get_Q_crossings(self, datapoints, property, rising):
        """Return data of Q=0 or Q=1 crossings of specified line."""
        num_points = len(datapoints['Q'])
        distance_to_gas = 1 - datapoints['Q']
        distance_to_liq = datapoints['Q']

        idx_gas = np.argwhere(
            (distance_to_gas < distance_to_liq) &
            (datapoints['Q'] > -1))[:, 0]
        idx_liq = np.argwhere(
            (distance_to_gas > distance_to_liq) &
            (datapoints['Q'] > -1))[:, 0]

        data = {}

        if idx_gas.size > 1:
            if rising is False:
                pos = -1
                next_point = -1
            else:
                pos = 0
                next_point = 1

            if idx_gas[pos] < num_points - 1 and idx_gas[pos] > 0:
                Q1 = datapoints['Q'][idx_gas[pos]]
                Q2 = datapoints['Q'][idx_gas[pos] + next_point]
                value1 = datapoints[property][idx_gas[pos]]
                value2 = datapoints[property][idx_gas[pos] + next_point]
                data[1] = value1 + (1 - Q1) / (Q2 - Q1) * (value2 - value1)

        if idx_liq.size > 1:
            if not rising or property == 'v':
                pos = 0
                next_point = 1
            else:
                pos = -1
                next_point = -1

            if idx_liq[pos] < num_points - 1 and idx_liq[pos] > 0:
                Q1 = datapoints['Q'][idx_liq[pos]]
                Q2 = datapoints['Q'][idx_liq[pos] + next_point]
                prop1 = datapoints[property][idx_liq[pos]]
                prop2 = datapoints[property][idx_liq[pos] + next_point]
                data[0] = prop1 + (0 - Q1) / (Q2 - Q1) * (prop2 - prop1)

        return data

    def _insert_Q_crossings(self, datapoints, property, data):
        """Insert data of Q=0 and Q=1 crossings into specified line."""
        rising = datapoints['s'][0] < datapoints['s'][-1]
        for Q, value in data.items():
            if property == 'v':
                value = 1 / value
            try:
                self.state.update(*CP.CoolProp.generate_update_pair(
                    self.CoolProp_inputs[property], value, CP.iQ, Q
                ))
            except ValueError:
                continue

            s = self.state.smass()

            if rising:
                position = np.searchsorted(datapoints['s'], s)
            else:
                position = np.searchsorted(-datapoints['s'], -s)

            datapoints['s'] = np.insert(datapoints['s'], position, s)
            datapoints[property] = np.insert(
                datapoints[property], position, value)
            datapoints['Q'] = np.insert(datapoints['Q'], position, Q)

            for prop in ['h', 'T', 'v', 'p']:
                if prop != property:
                    result = self.state.keyed_output(
                        self.CoolProp_inputs[prop])
                    if prop == 'v':
                        result = 1 / result
                    datapoints[prop] = np.insert(
                        datapoints[prop], position, result
                    )

        return datapoints

    def draw_isolines(self, fig, ax, diagram_type, x_min, x_max, y_min, y_max, isoline_data=None):
        """_summary_

        Parameters
        ----------
        fig : matplotlib.pyplot.figure
            Figure to draw into

        ax : matplotlib.pyplot.axes
            Axes to draw into

        diagram_type : str
            Name of the diagram

        x_min : number
            Minimum for x range

        x_max : number
            Maximum for x range

        y_min : number
            Minimum for y range

        y_max : number
            Maximum for y range

        isoline_data : dict, optional
            Dictionary holding additional data on the isolines to be drawn,
            by default None. These are

            - the isoline values with key :code:`values` and
            - the isoline style with key :code:`style`.

            The islonline style is another dictionary holding keyword arguments
            of a :code:`matplotlib.lines.Line2D` object. See
            https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html
            for more information.
        """
        if isoline_data is None:
            isoline_data = {}

        self._check_diagram_types(diagram_type)

        x_scale = self.supported_diagrams[diagram_type]['x_scale']
        y_scale = self.supported_diagrams[diagram_type]['y_scale']

        x_property = self.supported_diagrams[diagram_type]['x_property']
        y_property = self.supported_diagrams[diagram_type]['y_property']

        ax.clear()

        ax.set_ylim([y_min, y_max])
        ax.set_xlim([x_min, x_max])

        ax.set_xscale(x_scale)
        ax.set_yscale(y_scale)

        x_label = f'{x_property} in {_beautiful_unit_string(self.units[x_property])}'
        y_label = f'{y_property} in {_beautiful_unit_string(self.units[y_property])}'

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        isolines = [
            x for x in self.properties.keys()
            if x not in [x_property, y_property]
        ]

        for isoline in isolines:

            property = self.properties[isoline]
            data = getattr(self, property)

            isovalues = data['isolines']

            if isoline in isoline_data.keys():
                keys = isoline_data[isoline].keys()
                if 'style' in keys:
                    data['style'].update(isoline_data[isoline]['style'])

                if 'values' in keys:
                    isovalues = self.convert_to_SI(
                        isoline_data[isoline]['values'], isoline
                    )

                if 'label_position' in keys:
                    data['label_position'] = (
                        isoline_data[isoline]['label_position']
                    )

            for isoval in isovalues.round(8):
                if isoval not in data['isolines']:
                    msg = (
                        f'Could not find data for {property} isoline with '
                        f'value: {isoval}.'
                    )
                    continue

                x = self.convert_from_SI(data[isoval][x_property], x_property)
                y = self.convert_from_SI(data[isoval][y_property], y_property)

                indices = np.intersect1d(
                    np.where((x >= x_min) & (x <= x_max)),
                    np.where((y >= y_min) & (y <= y_max))
                )

                if len(indices) == 0:
                    continue

                gap = np.where(np.diff(indices) > 1)[0]
                if len(gap) > 0:
                    indices = np.insert(indices, gap + 1, indices[gap] + 1)
                    indices = np.insert(indices, gap + 2, indices[gap + 2] - 1)

                if indices[0] != 0:
                    indices = np.insert(indices, 0, indices[0] - 1)
                if indices[-1] < len(x) - 1:
                    indices = np.append(indices, indices[-1] + 1)

                y = y[indices]
                x = x[indices]

                ax.plot(x, y, **data['style'])

                isoval = self.convert_from_SI(isoval, isoline)

                self._draw_isoline_label(
                    fig, ax,
                    isoval.round(8), isoline,
                    int(data['label_position'] * len(x)),
                    x, y, x_min, x_max, y_min, y_max
                )

    def _check_diagram_types(self, diagram_type):
        if not isinstance(diagram_type, str):
            msg = (
                'The diagram_type must be specified as string. Available '
                f'inputs are: {", ".join(self.supported_diagrams.keys())}.'
            )
            raise TypeError(msg)
        elif diagram_type not in self.supported_diagrams.keys():
            msg = (
                'The specified diagram_type is not available. Available '
                f'types are: {", ".join(self.supported_diagrams.keys())}.'
            )
            raise ValueError(msg)

    def convert_to_SI(self, value, property):
        """Convert a value to its SI value."""
        if property == 'T':
            return (
                (
                    value + self.converters[property][self.units[property]][0]
                ) * self.converters[property][self.units[property]][1]
            )
        else:
            return value * self.converters[property][self.units[property]]

    def convert_from_SI(self, value, property):
        """Convert a SI value to value in respecive unit system."""
        if property == 'T':
            return (
                value / self.converters[property][self.units[property]][1]
                - self.converters[property][self.units[property]][0]
            )
        else:
            return value / self.converters[property][self.units[property]]
