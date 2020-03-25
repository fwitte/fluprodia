import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy.interpolate
from matplotlib.patches import Rectangle


def beautiful_unit_string(unit):
    if '/' in unit:
        numerator = unit.split('/')[0]
        denominator = unit.split('/')[1]
        unit = '$\\frac{' + numerator + '}{' + denominator + '}$'

    return unit


def IF97_Dmass(
    rhotarget, state, inputfix,
    inputfix_value, inputvar, inputvar_value, val_min, val_max):
    iter = 0
    inputvar_value = max(min(inputvar_value, val_max * 0.999), val_min * 1.001)
    while True:
        inputpair = CP.CoolProp.generate_update_pair(
            inputfix, inputfix_value, inputvar, inputvar_value)
        state.update(*inputpair)
        residual = rhotarget - state.rhomass()

        d = 0.01
        inputpair = CP.CoolProp.generate_update_pair(
            inputfix, inputfix_value, inputvar, inputvar_value + d)
        state.update(*inputpair)
        rho1 = state.rhomass()
        inputpair = CP.CoolProp.generate_update_pair(
            inputfix, inputfix_value, inputvar, inputvar_value - d)
        state.update(*inputpair)
        rho2 = state.rhomass()

        derivative = -(rho1 - rho2) / (2 * d)

        inputvar_value -= residual / derivative
        inputvar_value = max(min(inputvar_value, val_max * 0.999), val_min * 1.001)

        if abs(residual) < 1e-10:
            return inputvar_value
        elif iter > 10:
            return np.nan
        else:
            iter += 1


class StatesDiagram:

    def __init__(self, fluid, backend='HEOS', width=16, height=10):
        self.state = CP.AbstractState(backend, fluid)
        self.backend = backend
        self.fluid = fluid
        self.converters = {}
        self.converters['p'] = {'Pa': 1, 'mbar': 1e2, 'bar': 1e5, 'MPa': 1e6}
        self.converters['T'] = {'K': 0, '°C': 273.15}
        self.converters['s'] = {'J/kgK': 1, 'kJ/kgK': 1e3, 'MJ/kgK': 1e6}
        self.converters['h'] = {'J/kg': 1, 'kJ/kg': 1e3, 'MJ/kg': 1e6}
        self.converters['v'] = {'m^3/kg': 1, 'l/kg': 1e-3}
        self.converters['Q'] = {'%': 0.01}
        self.properties = {
            'p': 'pressure',
            'v': 'volume',
            'T': 'temperature',
            'h': 'enthalpy',
            's': 'entropy',
            'Q': 'quality'}
        self.units = {}
        self.supported_diagrams = ['Ts', 'hs', 'ph', 'pv']

        self.set_diagram_layout(width, height)
        self.set_unit_system()

        self.fig, self.ax = plt.subplots(figsize=(self.width, self.height))
        self.get_ax_size()

        self.pressure = {'isolines': np.array([])}
        self.entropy = {'isolines': np.array([])}
        self.temperature = {'isolines': np.array([])}
        self.enthalpy = {'isolines': np.array([])}
        self.volume = {'isolines': np.array([])}
        self.quality = {'isolines': np.array([])}

        p_min = self.state.trivial_keyed_output(CP.iP_min)
        T_max = self.state.trivial_keyed_output(CP.iT_max)
        self.state.update(CP.PT_INPUTS, p_min, T_max)
        s_max = self.state.smass()
        self.iterator = np.linspace(0, s_max, 100)

        self.set_isoline_defaults()
        self.set_limits()

    def set_isolines(self, **kwargs):
        keys = ['p', 'T', 'Q', 's', 'h', 'v']
        for key in kwargs:
            if key in keys:
                obj = getattr(self, self.properties[key])
                if key == 'T':
                    obj['isolines'] = kwargs[key] + self.converters[key][self.units[key]]
                else:
                    obj['isolines'] = kwargs[key] * self.converters[key][self.units[key]]
            else:
                msg = (
                    'The specified isoline \'' + key + '\' is not available. '
                    'Choose from: ' + str(keys) + '.')
                print(msg)

    def set_isoline_defaults(self):
        self.p_trip = self.state.trivial_keyed_output(CP.iP_triple)
        self.p_max = self.state.trivial_keyed_output(CP.iP_max)
        self.T_trip = self.state.trivial_keyed_output(CP.iT_triple)
        self.T_max = self.state.trivial_keyed_output(CP.iT_max)
        self.state.update(CP.PT_INPUTS, self.p_trip, self.T_max)
        self.v_max = 1 / self.state.rhomass()
        self.s_max = self.state.smass()
        self.state.update(CP.PT_INPUTS, self.p_max, self.T_trip)
        self.v_min = 1 / self.state.rhomass()

        self.p_crit = self.state.trivial_keyed_output(CP.iP_critical)
        self.T_crit = self.state.trivial_keyed_output(CP.iT_critical)
        self.v_crit = 1 / self.state.trivial_keyed_output(CP.irhomass_critical)

        self.pressure['isolines'] = np.geomspace(self.p_trip, self.p_max, 11)
        self.temperature['isolines'] = np.linspace(self.T_trip, self.T_max, 11)
        self.quality['isolines'] = np.linspace(0, 1, 11).round(1)
        self.entropy['isolines'] = np.linspace(0, self.s_max, 11)
        self.enthalpy = {'isolines': np.array([])}
        self.volume = {'isolines': np.geomspace(self.v_min, self.v_max, 11)}

    def set_limits(self, x_min=None, x_max=None, y_min=None, y_max=None):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

    def set_diagram_layout(self, width, height):
        self.width = width
        self.height = height

    def set_unit_system(self, p_unit='Pa', T_unit='K', s_unit='J/kgK', h_unit='J/kg', v_unit='m^3/kg', Q_unit='%'):
        self.units['p'] = p_unit
        self.units['T'] = T_unit
        self.units['s'] = s_unit
        self.units['h'] = h_unit
        self.units['v'] = v_unit
        self.units['Q'] = Q_unit

    def save(self, filename='StatesDiagram.pdf', **kwargs):
        self.ax.set_xlim([self.x_min, self.x_max])
        self.ax.set_ylim([self.y_min, self.y_max])
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel(self.y_label)
        self.ax.grid()
        plt.tight_layout()
        self.fig.savefig(filename, **kwargs)

    def draw_isoline_label(self, isoline, property, idx, x, y):

        if idx > len(x):
            return

        alpha = np.arctan(
            (y[idx] - y[idx - 1]) * (self.height / (self.y_max - self.y_min)) /
            ((x[idx] - x[idx - 1]) * (self.width / (self.x_max - self.x_min)))
        ) / (2 * np.pi) * 360

        unit = beautiful_unit_string(self.units[property])

        txt = str(isoline) + ' ' + unit
        self.ax.text(
            x[idx], y[idx], txt, fontsize=4,
            rotation=alpha, va='center', ha='center',
            bbox=dict(facecolor='white', edgecolor='white', pad=0.0))

    def get_ax_size(self):
        bbox = self.ax.get_window_extent().transformed(
            self.fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        self.width *= self.fig.dpi
        self.height *= self.fig.dpi

    def isobar(self, isolines=None):
        if isolines is None:
            isolines = self.pressure['isolines']
        else:
            self.pressure['isolines'] = isolines

        for p in isolines:
            self.pressure[p] = {'h': [], 'T': [], 'v': [], 's': [], 'p': []}
            for val in self.iterator:
                try:
                    self.state.update(CP.PSmass_INPUTS, p, val)
                    self.pressure[p]['h'] += [self.state.hmass()]
                    self.pressure[p]['T'] += [self.state.T()]
                    self.pressure[p]['v'] += [1 / self.state.rhomass()]
                    self.pressure[p]['s'] += [val]
                    self.pressure[p]['p'] += [p]
                except ValueError:
                    continue

            self.pressure[p]['h'] = np.array(self.pressure[p]['h'])
            self.pressure[p]['T'] = np.array(self.pressure[p]['T'])
            self.pressure[p]['v'] = np.array(self.pressure[p]['v'])
            self.pressure[p]['s'] = np.array(self.pressure[p]['s'])
            self.pressure[p]['p'] = np.array(self.pressure[p]['p'])

            if p <= self.p_crit:
                for Q in [0, 1]:
                    self.state.update(CP.PQ_INPUTS, p, Q)
                    s = self.state.smass()

                    idx = np.searchsorted(self.pressure[p]['s'], s)
                    self.pressure[p]['h'] = np.insert(
                        self.pressure[p]['h'], idx, self.state.hmass())
                    self.pressure[p]['T'] = np.insert(
                        self.pressure[p]['T'], idx, self.state.T())
                    self.pressure[p]['v'] = np.insert(
                        self.pressure[p]['v'], idx, 1 / self.state.rhomass())
                    self.pressure[p]['s'] = np.insert(
                        self.pressure[p]['s'], idx, s)
                    self.pressure[p]['p'] = np.insert(
                        self.pressure[p]['p'], idx, s)

    def isochor(self, isolines=None):
        if isolines is None:
            isolines = self.volume['isolines']
        else:
            self.volume['isolines'] = isolines

        for v in isolines:
            self.volume[v] = {'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in self.iterator:
                try:
                    self.state.update(CP.DmassSmass_INPUTS, 1 / v, val)
                    self.volume[v]['h'] += [self.state.hmass()]
                    self.volume[v]['T'] += [self.state.T()]
                    self.volume[v]['p'] += [self.state.p()]
                    self.volume[v]['s'] += [val]
                    self.volume[v]['v'] += [v]
                except ValueError:
                    continue

            self.volume[v]['h'] = np.array(self.volume[v]['h'])
            self.volume[v]['T'] = np.array(self.volume[v]['T'])
            self.volume[v]['p'] = np.array(self.volume[v]['p'])
            self.volume[v]['s'] = np.array(self.volume[v]['s'])
            self.volume[v]['v'] = np.array(self.volume[v]['v'])

            if v >= self.v_crit:
                try:
                    self.state.update(CP.DmassQ_INPUTS, 1 / v, 1)
                    s = self.state.smass()

                    idx = np.searchsorted(self.volume[v]['s'], s)
                    self.volume[v]['h'] = np.insert(
                        self.volume[v]['h'], idx, self.state.hmass())
                    self.volume[v]['T'] = np.insert(
                        self.volume[v]['T'], idx, self.state.T())
                    self.volume[v]['p'] = np.insert(
                        self.volume[v]['p'], idx, self.state.p())
                    self.volume[v]['s'] = np.insert(
                        self.volume[v]['s'], idx, s)
                    self.volume[v]['v'] = np.insert(
                        self.volume[v]['v'], idx, s)
                except ValueError:
                    continue

    def isoquality(self, isolines=None):
        if isolines is None:
            isolines = self.quality['isolines']
        else:
            self.quality['isolines'] = isolines

        iterator = np.append(
            np.linspace(self.T_trip, self.T_crit * 0.97, 40, endpoint=False),
            np.linspace(self.T_crit * 0.97, self.T_crit, 40))

        for Q in isolines:
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

            self.quality[Q]['h'] = np.array(self.quality[Q]['h'])
            self.quality[Q]['p'] = np.array(self.quality[Q]['p'])
            self.quality[Q]['v'] = np.array(self.quality[Q]['v'])
            self.quality[Q]['s'] = np.array(self.quality[Q]['s'])
            self.quality[Q]['T'] = np.array(self.quality[Q]['T'])

    def isoenthalpy(self, isolines=None):
        if isolines is None:
            isolines = self.enthalpy['isolines']
        else:
            self.enthalpy['isolines'] = isolines

        iterator = np.geomspace(
            self.pressure['isolines'].min(),
            self.pressure['isolines'].max(), 100)

        for h in isolines:
            self.enthalpy[h] = {
                'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in iterator:
                try:
                    self.state.update(CP.HmassP_INPUTS, h, val)
                    self.enthalpy[h]['T'] += [self.state.T()]
                    self.enthalpy[h]['p'] += [self.state.p()]
                    self.enthalpy[h]['v'] += [1 / self.state.rhomass()]
                    self.enthalpy[h]['s'] += [self.state.smass()]
                    self.enthalpy[h]['h'] += [h]
                except ValueError:
                    continue

            self.enthalpy[h]['T'] = np.array(self.enthalpy[h]['T'])
            self.enthalpy[h]['p'] = np.array(self.enthalpy[h]['p'])
            self.enthalpy[h]['v'] = np.array(self.enthalpy[h]['v'])
            self.enthalpy[h]['s'] = np.array(self.enthalpy[h]['s'])
            self.enthalpy[h]['h'] = np.array(self.enthalpy[h]['h'])

            for Q in [0, 1]:
                try:
                    self.state.update(CP.HmassQ_INPUTS, h, Q)
                    s = self.state.smass()

                    idx = np.searchsorted(self.enthalpy[h]['s'], s)
                    self.enthalpy[h]['h'] = np.insert(
                        self.enthalpy[h]['h'], idx, h)
                    self.enthalpy[h]['T'] = np.insert(
                        self.enthalpy[h]['T'], idx, self.state.T())
                    self.enthalpy[h]['p'] = np.insert(
                        self.enthalpy[h]['p'], idx, self.state.p())
                    self.enthalpy[h]['s'] = np.insert(
                        self.enthalpy[h]['s'], idx, s)
                    self.enthalpy[h]['v'] = np.insert(
                        self.enthalpy[h]['v'], idx, 1 / self.state.rhomass())
                except ValueError:
                    continue

    def isotherm(self, isolines=None):
        if isolines is None:
            isolines = self.temperature['isolines']
        else:
            self.temperature['isolines'] = isolines

        iterator = np.geomspace(
            self.volume['isolines'].min(),
            self.volume['isolines'].max(), 150)

        for T in isolines:
            self.temperature[T] = {
                'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in iterator:
                try:
                    self.state.update(CP.DmassT_INPUTS, 1 / val, T)
                    self.temperature[T]['T'] += [T]
                    self.temperature[T]['p'] += [self.state.p()]
                    self.temperature[T]['v'] += [1 / self.state.rhomass()]
                    self.temperature[T]['s'] += [self.state.smass()]
                    self.temperature[T]['h'] += [self.state.hmass()]
                except ValueError:
                    continue

            self.temperature[T]['T'] = np.array(self.temperature[T]['T'])
            self.temperature[T]['p'] = np.array(self.temperature[T]['p'])
            self.temperature[T]['v'] = np.array(self.temperature[T]['v'])
            self.temperature[T]['s'] = np.array(self.temperature[T]['s'])
            self.temperature[T]['h'] = np.array(self.temperature[T]['h'])

            if T <= self.T_crit:
                for Q in [0, 1]:
                    try:
                        self.state.update(CP.QT_INPUTS, Q, T)
                        s = self.state.smass()

                        idx = np.searchsorted(self.temperature[T]['s'], s)
                        self.temperature[T]['h'] = np.insert(
                            self.temperature[T]['h'], idx, self.state.hmass())
                        self.temperature[T]['T'] = np.insert(
                            self.temperature[T]['T'], idx, T)
                        self.temperature[T]['p'] = np.insert(
                            self.temperature[T]['p'], idx, self.state.p())
                        self.temperature[T]['s'] = np.insert(
                            self.temperature[T]['s'], idx, s)
                        self.temperature[T]['v'] = np.insert(
                            self.temperature[T]['v'], idx, 1 / self.state.rhomass())
                    except ValueError:
                        continue

    def isoentropy(self, isolines=None):
        if isolines is None:
            isolines = self.entropy['isolines']
        else:
            self.entropy['isolines'] = isolines

        iterator = np.linspace(
            self.temperature['isolines'].min(),
            self.temperature['isolines'].max(), 100)

        for s in isolines:
            self.entropy[s] = {
                'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in iterator:
                try:
                    self.state.update(CP.SmassT_INPUTS, s, val)
                    self.entropy[s]['T'] += [val]
                    self.entropy[s]['p'] += [self.state.p()]
                    self.entropy[s]['v'] += [1 / self.state.rhomass()]
                    self.entropy[s]['s'] += [s]
                    self.entropy[s]['h'] += [self.state.hmass()]
                except ValueError:
                    continue

            self.entropy[s]['T'] = np.array(self.entropy[s]['T'])
            self.entropy[s]['p'] = np.array(self.entropy[s]['p'])
            self.entropy[s]['v'] = np.array(self.entropy[s]['v'])
            self.entropy[s]['s'] = np.array(self.entropy[s]['s'])
            self.entropy[s]['h'] = np.array(self.entropy[s]['h'])

    def draw_isolines(self, isolines, diagram_type=None):
        if diagram_type not in self.supported_diagrams:
            msg = 'This diagram is not supported.'
            raise ValueError(msg)

        x_property = diagram_type[1]
        y_property = diagram_type[0]

        self.x_label = x_property + ' in ' + beautiful_unit_string(self.units[x_property])
        self.y_label = y_property + ' in ' + beautiful_unit_string(self.units[y_property])

        for isoline in isolines:

            property = self.properties[isoline]
            data = getattr(self, property)

            isoline_conv = self.converters[isoline][self.units[isoline]]
            x_conv = self.converters[x_property][self.units[x_property]]
            y_conv = self.converters[y_property][self.units[y_property]]

            for isoval in data['isolines']:
                if x_property == 'T':
                    x = data[isoval][x_property] - x_conv
                else:
                    x = data[isoval][x_property] / x_conv

                if y_property == 'T':
                    y = data[isoval][y_property] - y_conv
                else:
                    y = data[isoval][y_property] / y_conv

                indices = np.intersect1d(
                    np.where((x >= self.x_min) & (x <= self.x_max)),
                    np.where((y >= self.y_min) & (y <= self.y_max))
                )

                if len(indices) > 0:
                    if indices[0] != 0:
                        indices = np.insert(indices, 0, indices[0] - 1)
                    if indices[-1] < len(x) - 1:
                        indices = np.append(indices, indices[-1] + 1)

                    y = y[indices]
                    x = x[indices]

                    self.ax.plot(
                        x, y,
                        linestyle='-.', color='#363636', linewidth=0.5)
                    idx = 20
                    if isoline == 'T':
                        isoval -= isoline_conv
                    else:
                        isoval /= isoline_conv
                    self.draw_isoline_label(isoval, isoline, idx, x, y)
