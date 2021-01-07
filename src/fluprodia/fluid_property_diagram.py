# -*- coding: utf-8

"""Module for fluid property diagram creation.

This file is part of project fluprodia (github.com/fwitte/fluprodia). It's
copyrighted by the contributors recorded in the version control history of the
file, available from its original location
src/fluprodia/fluid_property_diagram.py

SPDX-License-Identifier: MIT
"""

import CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np


def beautiful_unit_string(unit):
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


def isolines_log(val_min, val_max):
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
    and export to Ts-, hs- and logph-diagram. The default values for width and
    height are 16/10.

    >>> from fluprodia import FluidPropertyDiagram
    >>> import numpy as np
    >>> diagram = FluidPropertyDiagram('water')
    >>> diagram.width
    16.0
    >>> diagram.height
    10.0

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

    Now we can calculate the diagram data and export after that. If you want to
    plot additional data on the diagram (e.g. measurement data from a
    thermodynamic process) you can do this on the :code:`diagram.ax` object. It
    is a :code:`matplotlib.axes._subplots.AxesSubplot` object and therefore
    you can apply standard matplotlib methods. The figure of the plot can be
    accessed via :code:`diagram.fig`

    >>> type(diagram.ax)
    <class 'matplotlib.axes._subplots.AxesSubplot'>
    >>> type(diagram.fig)
    <class 'matplotlib.figure.Figure'>
    >>> diagram.calc_isolines()

    After that it is possible to specify the view range of the plot and draw
    the isolines of a specific type of diagram. Last step is to export your
    diagram. Any file format supported by matplotlib is possible.

    >>> diagram.set_limits(x_min=0, x_max=8, y_min=0, y_max=700)
    >>> diagram.draw_isolines(diagram_type='Ts')
    >>> diagram.save('Ts_Diagramm.pdf')

    If we want to create a different diagram, e.g. hs-diagram, it is not
    necessary to recalculate the isolines. Instead, reset the limits to match
    your needs and draw the isolines for a different diagram. If limits do not
    change, there is no necessety to respecify the limits.

    >>> diagram.set_limits(y_min=0, y_max=3600)
    >>> diagram.draw_isolines(diagram_type='hs')
    >>> diagram.save('hs_Diagramm.pdf')

    >>> diagram.set_limits(x_min=0, x_max=3600, y_min=1e-2, y_max=5e2)
    >>> diagram.draw_isolines(diagram_type='logph')
    >>> diagram.save('logph_Diagramm.pdf')

    It is also possible to specify/modify the isolines to plot. For example,
    the lines of constant specific enthalpy should be plotted in red color
    and with a linewidth of 2. Also, the lines of constant specific volumen
    should not be plotted at all. For more information see the
    :py:meth:`fluprodia.fluid_property_diagram.FluidPropertyDiagram.draw_isolines`
    method.

    >>> diagram.set_limits(x_min=0, x_max=8, y_min=0, y_max=700)
    >>> diagram.draw_isolines(diagram_type='Ts',
    ... isoline_data={'h': {
    ... 'values': iso_h,
    ... 'style': {'linewidth': 2, 'color': '#ff0000'}},
    ... 'v': {'values': np.array([])}})
    >>> diagram.save('Ts_Diagramm.pdf')
    """

    def __init__(self, fluid, width=16.0, height=10.0):
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
            'bar': 1e5, 'MPa': 1e6}
        self.converters['T'] = {
            'K': [0, 1], '°C': [273.15, 1], '°F': [459.67, 5 / 9]}
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

        self.single_isoline_functions = {
            'p': self.single_isobaric,
            'v': self.single_isochoric,
            'T': self.single_isothermal,
            'h': self.single_isenthalpic,
            's': self.single_isentropic
        }

        self.CoolProp_inputs = {
            'p': CP.iP,
            'v': CP.iDmass,
            'T': CP.iT,
            'h': CP.iHmass,
            's': CP.iSmass
        }

        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None

        self.set_diagram_layout(width, height)
        self.set_unit_system()

        self.pressure = {'isolines': np.array([])}
        self.entropy = {'isolines': np.array([])}
        self.temperature = {'isolines': np.array([])}
        self.enthalpy = {'isolines': np.array([])}
        self.volume = {'isolines': np.array([])}
        self.quality = {'isolines': np.array([])}

        self.set_isoline_defaults()
        self.set_limits()
        self.default_line_layout()
        self.default_label_positioning()

    def default_line_layout(self):
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

    def default_label_positioning(self):
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
                    'The specified isoline \'' + key + '\' is not available. '
                    'Choose from: ' + str(keys) + '.')
                print(msg)

    def set_isoline_defaults(self):
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

        self.state.update(CP.PQ_INPUTS, self.p_trip + 1, 0)
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

        self.pressure['isolines'] = isolines_log(
            self.p_trip + 1e-2, self.p_max).round(8)
        self.volume['isolines'] = isolines_log(self.v_min, self.v_max).round(8)

    def set_limits(self, **kwargs):
        """Set the diagram's limits.

        Parameters
        ----------
        x_min : float
            Minimum value for x axis :code:`x_min`.

        x_max : float
            Maximum value for x axis :code:`x_max`.

        y_min : float
            Minimum value for y axis :code:`y_min`.

        y_max : float
            Maximum value for y axis :code:`y_max`.
        """
        self.x_min = kwargs.get('x_min', self.x_min)
        self.x_max = kwargs.get('x_max', self.x_max)
        self.y_min = kwargs.get('y_min', self.y_min)
        self.y_max = kwargs.get('y_max', self.y_max)

    def set_diagram_layout(self, width, height):
        """Set the diagram's width and height and create the figure.

        Parameters
        ----------
        width : float
            Width of all diagrams.

        height : float
            Height of all diagrams.
        """
        self.width = width
        self.height = height
        self.fig = plt.figure(figsize=(self.width, self.height))
        self.ax = self.fig.add_subplot()

    def set_unit_system(self, **kwargs):
        u"""Set the unit system for the fluid properties.

        Parameters
        ----------
        p : str
            Unit of pressure, units available are
            :code:`Pa, hPa, mbar, psi, bar, MPa`.

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
                    'The unit ' + str(value) + ' is not available for the '
                    'fluid property ' + self.properties[key] + '.')
                raise ValueError(msg)

    def save(self, filename='FluidPropertyDiagram.pdf', **kwargs):
        """Save the diagram.

        Parameters
        ----------
        filename : str
            Path for the exported diagram. The file extension determines the
            file format of the diagram. Available formats are the standard
            matplotlib formats.

        kwargs : misc
            Keyword arguments of the :code:`matplotlib.figure.Figure.savefig`
            method.
        """
        self.ax.set_xlim([self.x_min, self.x_max])
        self.ax.set_ylim([self.y_min, self.y_max])
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel(self.y_label)
        self.ax.grid()
        plt.tight_layout()
        self.fig.savefig(filename, **kwargs)

    def draw_isoline_label(self, isoline, property, idx, x, y):
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
        if (idx > len(x) or
                idx == 0 or
                x[idx] > self.x_max or x[idx] < self.x_min or
                y[idx] > self.y_max or y[idx] < self.y_min or
                x[idx - 1] > self.x_max or x[idx - 1] < self.x_min or
                y[idx - 1] > self.y_max or y[idx - 1] < self.y_min):
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
            if self.ax.get_xscale() == 'log':
                x_scaled = (np.log(x[idx]) - np.log(x[idx - 1])) * (
                    self.width / (np.log(self.x_max) - np.log(self.x_min)))
            else:
                x_scaled = (x[idx] - x[idx - 1]) * (
                    self.width / (self.x_max - self.x_min))

            if self.ax.get_yscale() == 'log':
                y_scaled = (np.log(y[idx]) - np.log(y[idx - 1])) * (
                    self.height / (np.log(self.y_max) - np.log(self.y_min)))
            else:
                y_scaled = (y[idx] - y[idx - 1]) * (
                    self.height / (self.y_max - self.y_min))

            alpha = np.arctan(y_scaled / x_scaled) / (2 * np.pi) * 360

        unit = beautiful_unit_string(self.units[property])

        txt = str(isoline) + ' ' + unit
        self.ax.text(
            x[idx], y[idx], txt, fontsize=5,
            rotation=alpha, va='center', ha='center',
            bbox=dict(facecolor='white', edgecolor='white', pad=0.0))

    def calc_isolines(self):
        """Calculate all isolines."""
        self.isobaric()
        self.isochoric()
        self.isoquality()
        self.isenthalpic()
        self.isothermal()
        self.isentropic()

    def isobaric(self):
        """Calculate an isoline of constant pressure."""
        isolines = self.pressure['isolines']

        iterator = 1 / np.append(
            np.geomspace(self.v_min, self.v_intermediate, 100, endpoint=False),
            np.geomspace(self.v_intermediate, self.v_max, 100))

        for p in isolines.round(8):
            self.pressure[p] = self.single_isobaric(iterator, np.ones(len(iterator)) * p)

    def isochoric(self):
        """Calculate an isoline of constant specific volume."""
        isolines = self.volume['isolines']

        iterator = np.append(
            np.geomspace(self.p_trip, self.p_crit * 0.8, 100, endpoint=False),
            np.geomspace(self.p_crit * 0.8, self.p_max, 100))

        for v in isolines.round(8):
            self.volume[v] = self.single_isochoric(
                iterator, np.ones(len(iterator)) * 1 / v)

    def isothermal(self):
        """Calculate an isoline of constant temperature."""
        isolines = self.temperature['isolines']

        iterator = np.linspace(self.s_min, self.s_max, 200)

        for T in isolines.round(8):
            self.temperature[T] = self.single_isothermal(
                iterator, np.ones(len(iterator)) * T)

    def isoquality(self):
        """Calculate an isoline of constant vapor mass fraction."""
        isolines = self.quality['isolines']

        iterator = np.append(
            np.linspace(self.T_trip, self.T_crit * 0.97, 40, endpoint=False),
            np.linspace(self.T_crit * 0.97, self.T_crit, 40))

        for Q in isolines.round(8):
            self.quality[Q] = {
                'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in iterator:
                try:
                    self.state.update(CP.QT_INPUTS, Q, val)
                    self.quality[Q]['T'] += [val]
                    self.quality[Q]['h'] += [self.state.hmass()]
                    self.quality[Q]['p'] += [self.state.p()]
                    self.quality[Q]['v'] += [1 / self.state.rhomass()]
                    self.quality[Q]['s'] += [self.state.smass()]
                except ValueError:
                    continue

            self.quality[Q]['h'] = np.asarray(self.quality[Q]['h'])
            self.quality[Q]['p'] = np.asarray(self.quality[Q]['p'])
            self.quality[Q]['v'] = np.asarray(self.quality[Q]['v'])
            self.quality[Q]['s'] = np.asarray(self.quality[Q]['s'])
            self.quality[Q]['T'] = np.asarray(self.quality[Q]['T'])

    def isenthalpic(self):
        """Calculate an isoline of constant specific enthalpy."""
        isolines = self.enthalpy['isolines']

        iterator = 1 / np.append(
            np.geomspace(self.v_min, self.v_intermediate, 100, endpoint=False),
            np.geomspace(self.v_intermediate, self.v_max, 100))

        for h in isolines.round(8):
            self.enthalpy[h] = self.single_isenthalpic(
                    iterator, np.ones(len(iterator)) * h)

    def isentropic(self):
        """Calculate an isoline of constant specific entropy."""
        isolines = self.entropy['isolines']

        iterator = 1 / np.append(
            np.geomspace(self.v_min, self.v_intermediate, 100, endpoint=False),
            np.geomspace(self.v_intermediate, self.v_max, 100))

        for s in isolines.round(8):
            self.entropy[s] = self.single_isentropic(
                        iterator, np.ones(len(iterator)) * s)

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

    def single_isobaric(self, iterator, p):
        """Calculate an isoline of constant pressure."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': [], 'Q': []}

        if iterator[0] < iterator[-1]:
            rising = True
        else:
            rising = False

        i = 0
        for val in iterator:
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
            i += 1

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        data = self.get_Q_crossings(datapoints, 'p', rising)
        return self.insert_Q_crossings(datapoints, 'p', data)

    def single_isochoric(self, iterator, D):
        """Calculate an isoline of constant specific volume."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': [], 'Q': []}

        if iterator[0] < iterator[-1]:
            rising = True
        else:
            rising = False

        i = 0
        for val in iterator:
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
            i += 1

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        data = self.get_Q_crossings(datapoints, 'v', rising)
        return self.insert_Q_crossings(datapoints, 'v', data)

    def single_isothermal(self, iterator, T):
        """Calculate an isoline of constant temperature."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': [], 'Q': []}

        if iterator[0] < iterator[-1]:
            rising = True
        else:
            rising = False

        i = 0
        for val in iterator:
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
                        'P', 'T', T[i], 'S', val, self.fluid)
                    self.state.update(CP.PSmass_INPUTS, val, p)
                    datapoints['T'] += [T[i]]
                    datapoints['p'] += [p]
                    datapoints['v'] += [1 / self.state.rhomass()]
                    datapoints['s'] += [val]
                    datapoints['h'] += [self.state.hmass()]
                except ValueError:
                    pass
            i += 1

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        data = self.get_Q_crossings(datapoints, 'T', rising)
        return self.insert_Q_crossings(datapoints, 'T', data)

    def single_isenthalpic(self, iterator, h):
        """Calculate an isoline of constant specific enthalpy."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': []}

        i = 0
        for val in iterator:
            try:
                self.state.update(CP.DmassHmass_INPUTS, val, h[i])
                datapoints['T'] += [self.state.T()]
                datapoints['p'] += [self.state.p()]
                datapoints['v'] += [1 / val]
                datapoints['s'] += [self.state.smass()]
                datapoints['h'] += [h[i]]
            except ValueError:
                pass

            i += 1

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        return datapoints

    def single_isentropic(self, iterator, s):
        """Calculate an isoline of constant specific entropy."""
        datapoints = {'h': [], 'T': [], 'v': [], 's': [], 'p': []}

        i = 0
        for val in iterator:
            try:
                self.state.update(CP.DmassSmass_INPUTS, val, s[i])
                datapoints['T'] += [self.state.T()]
                datapoints['p'] += [self.state.p()]
                datapoints['v'] += [1 / val]
                datapoints['s'] += [s[i]]
                datapoints['h'] += [self.state.hmass()]
            except ValueError:
                pass

            i += 1

        for key in datapoints.keys():
            datapoints[key] = np.asarray(datapoints[key])

        return datapoints

    def get_Q_crossings(self, datapoints, property, rising):
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

    def insert_Q_crossings(self, datapoints, property, data):
        """Insert data of Q=0 and Q=1 crossings into specified line."""
        rising = datapoints['s'][0] < datapoints['s'][-1]
        for Q, value in data.items():
            if property == 'v':
                value = 1 / value
            try:
                self.state.update(*CP.CoolProp.generate_update_pair(
                    self.CoolProp_inputs[property], value, CP.iQ, Q))
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

    def draw_isolines(self, diagram_type, isoline_data={}):
        """Draw the isolines of a specific diagram type.

        Parameters
        ----------
        diagram_type : str
            Which type of diagram should be drawn.

        isoline_data : dict
            Dictionary holding additional data on the isolines to be drawn.
            These are

            - the isoline values with key :code:`values` and
            - the isoline style with key :code:`style`.

            The islonline style is another dictionary holding keyword arguments
            of a :code:`matplotlib.lines.Line2D` object. See
            https://matplotlib.org/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
            for more information.
        """
        if not isinstance(diagram_type, str):
            msg = (
                'The diagram_type must be specified as string. Available '
                'inputs are: ' + str(self.supported_diagrams.keys()) + '.')
            raise TypeError(msg)
        elif diagram_type not in self.supported_diagrams.keys():
            msg = (
                'The specified diagram_type is not available. Available '
                'inputs are: ' + str(self.supported_diagrams.keys()) + '.')
            raise ValueError(msg)

        self.ax.clear()

        x_property = self.supported_diagrams[diagram_type]['x_property']
        y_property = self.supported_diagrams[diagram_type]['y_property']
        x_scale = self.supported_diagrams[diagram_type]['x_scale']
        y_scale = self.supported_diagrams[diagram_type]['y_scale']
        self.ax.set_xscale(x_scale)
        self.ax.set_yscale(y_scale)

        self.x_label = (
            x_property + ' in ' +
            beautiful_unit_string(self.units[x_property]))
        self.y_label = (
            y_property + ' in ' +
            beautiful_unit_string(self.units[y_property]))

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
                        isoline_data[isoline]['values'], isoline)

                if 'label_position' in keys:
                    data['label_position'] = (
                        isoline_data[isoline]['label_position'])

            for isoval in isovalues.round(8):
                if isoval not in data['isolines']:
                    msg = (
                        'Could not find data for ' + property + ' isoline '
                        'with value: ' + str(isoval) + '.')
                    continue

                x = self.convert_from_SI(data[isoval][x_property], x_property)
                y = self.convert_from_SI(data[isoval][y_property], y_property)

                indices = np.intersect1d(
                    np.where((x >= self.x_min) & (x <= self.x_max)),
                    np.where((y >= self.y_min) & (y <= self.y_max))
                )

                if len(indices) == 0:
                    continue

                gap = np.where(np.diff(indices) > 1)[0]
                if len(gap) > 0:
                    indices = np.insert(
                        indices, gap + 1, indices[gap] + 1)
                    indices = np.insert(
                        indices, gap + 2, indices[gap + 2] - 1)

                if indices[0] != 0:
                    indices = np.insert(indices, 0, indices[0] - 1)
                if indices[-1] < len(x) - 1:
                    indices = np.append(indices, indices[-1] + 1)

                y = y[indices]
                x = x[indices]

                self.ax.plot(x, y, **data['style'])

                isoval = self.convert_from_SI(isoval, isoline)

                self.draw_isoline_label(
                    isoval.round(8), isoline,
                    int(data['label_position'] * len(x)), x, y)

    def convert_to_SI(self, value, property):
        """Convert a value to its SI value."""
        if property == 'T':
            return (
                (value +
                 self.converters[property][self.units[property]][0]) *
                self.converters[property][self.units[property]][1])
        else:
            return value * self.converters[property][self.units[property]]

    def convert_from_SI(self, value, property):
        """Convert a SI value to value in respecive unit system."""
        if property == 'T':
            return (
                value / self.converters[property][self.units[property]][1] -
                self.converters[property][self.units[property]][0])
        else:
            return value / self.converters[property][self.units[property]]
