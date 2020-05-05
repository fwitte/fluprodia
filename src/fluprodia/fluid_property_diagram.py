import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt


def beautiful_unit_string(unit):
    if '/' in unit:
        numerator = unit.split('/')[0]
        denominator = unit.split('/')[1]
        unit = '$\\frac{' + numerator + '}{' + denominator + '}$'

    return unit


def isolines_log(val_min, val_max):

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

    def __init__(self, fluid, width=16, height=10):
        self.state = CP.AbstractState('HEOS', fluid)
        self.fluid = fluid
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
        self.properties = {
            'p': 'pressure',
            'v': 'volume',
            'T': 'temperature',
            'h': 'enthalpy',
            's': 'entropy',
            'Q': 'quality'}
        self.units = {}
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
        self.pressure['label_position'] = 0.85
        self.volume['label_position'] = 0.7
        self.quality['label_position'] = 0.225
        self.entropy['label_position'] = 0.95
        self.enthalpy['label_position'] = 0.85
        self.temperature['label_position'] = 0.95

    def set_isolines(self, **kwargs):
        keys = ['p', 'T', 'Q', 's', 'h', 'v']
        for key in kwargs:
            if key in keys:
                obj = getattr(self, self.properties[key])
                if key == 'T':
                    obj['isolines'] = (
                        (kwargs[key] +
                         self.converters[key][self.units[key]][0]) *
                        self.converters[key][self.units[key]][1])
                else:
                    obj['isolines'] = (
                        kwargs[key] * self.converters[key][self.units[key]])
                obj['isolines'] = obj['isolines'].round(8)
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
        self.h_max = self.state.hmass()
        self.iterator = np.linspace(0, self.s_max, 100)
        self.state.update(CP.PT_INPUTS, self.p_max, self.T_trip)
        self.v_min = 1 / self.state.rhomass()

        self.p_crit = self.state.trivial_keyed_output(CP.iP_critical)
        self.T_crit = self.state.trivial_keyed_output(CP.iT_critical)
        self.v_crit = 1 / self.state.trivial_keyed_output(CP.irhomass_critical)

        self.quality['isolines'] = np.linspace(0, 1, 11).round(8)

        step = round(int(self.T_max - self.T_trip) / 15, -1)
        self.temperature['isolines'] = np.append(
            self.T_trip,
            np.arange(self.T_max, self.T_trip, -step)[::-1]).round(8)

        step = round(int(self.s_max) / 15, -1)
        self.entropy['isolines'] = np.arange(0, self.s_max, step).round(8)

        step = round(int(self.h_max) / 15, -1)
        self.enthalpy['isolines'] = np.arange(0, self.h_max, step).round(8)

        self.pressure['isolines'] = isolines_log(
            self.p_trip, self.p_max).round(8)
        self.volume['isolines'] = isolines_log(self.v_min, self.v_max).round(8)

    def set_limits(self, x_min=None, x_max=None, y_min=None, y_max=None):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

    def set_diagram_layout(self, width, height):
        self.width = width
        self.height = height
        self.fig = plt.figure(figsize=(self.width, self.height))
        self.ax = self.fig.add_subplot()

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

        if (x[idx] > self.x_max or x[idx] < self.x_min or
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
        self.isobar()
        self.isochor()
        self.isoquality()
        self.isoenthalpy()
        self.isotherm()
        self.isoentropy()

    def isobar(self):
        isolines = self.pressure['isolines']

        for p in isolines.round(8):
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

            self.pressure[p]['h'] = np.asarray(self.pressure[p]['h'])
            self.pressure[p]['T'] = np.asarray(self.pressure[p]['T'])
            self.pressure[p]['v'] = np.asarray(self.pressure[p]['v'])
            self.pressure[p]['s'] = np.asarray(self.pressure[p]['s'])
            self.pressure[p]['p'] = np.asarray(self.pressure[p]['p'])

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

    def isochor(self):
        isolines = self.volume['isolines']

        iterator = np.geomspace(self.p_trip, self.p_max, 100)

        for v in isolines.round(8):
            self.volume[v] = {'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in iterator:
                try:
                    self.state.update(CP.DmassP_INPUTS, 1 / v, val)
                    self.volume[v]['h'] += [self.state.hmass()]
                    self.volume[v]['T'] += [self.state.T()]
                    self.volume[v]['p'] += [val]
                    self.volume[v]['s'] += [self.state.smass()]
                    self.volume[v]['v'] += [v]
                except ValueError:
                    continue

            self.volume[v]['h'] = np.asarray(self.volume[v]['h'])
            self.volume[v]['T'] = np.asarray(self.volume[v]['T'])
            self.volume[v]['p'] = np.asarray(self.volume[v]['p'])
            self.volume[v]['s'] = np.asarray(self.volume[v]['s'])
            self.volume[v]['v'] = np.asarray(self.volume[v]['v'])

            for Q in [0, 1]:
                try:
                    self.state.update(CP.DmassQ_INPUTS, 1 / v, Q)
                    val = self.state.p()

                    idx = np.searchsorted(self.volume[v]['p'], val)
                    self.volume[v]['h'] = np.insert(
                        self.volume[v]['h'], idx, self.state.hmass())
                    self.volume[v]['T'] = np.insert(
                        self.volume[v]['T'], idx, self.state.T())
                    self.volume[v]['p'] = np.insert(
                        self.volume[v]['p'], idx, val)
                    self.volume[v]['s'] = np.insert(
                        self.volume[v]['s'], idx, self.state.smass())
                    self.volume[v]['v'] = np.insert(
                        self.volume[v]['v'], idx, v)
                except ValueError:
                    continue

    def isoquality(self):
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

    def isoenthalpy(self):
        isolines = self.enthalpy['isolines']

        iterator = np.geomspace(self.p_trip, self.p_max, 100)

        for h in isolines.round(8):
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

            self.enthalpy[h]['T'] = np.asarray(self.enthalpy[h]['T'])
            self.enthalpy[h]['p'] = np.asarray(self.enthalpy[h]['p'])
            self.enthalpy[h]['v'] = np.asarray(self.enthalpy[h]['v'])
            self.enthalpy[h]['s'] = np.asarray(self.enthalpy[h]['s'])
            self.enthalpy[h]['h'] = np.asarray(self.enthalpy[h]['h'])

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

    def isotherm(self):
        isolines = self.temperature['isolines']

        iterator = np.geomspace(self.p_trip, self.p_max, 300)

        for T in isolines.round(8):
            self.temperature[T] = {
                'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in iterator:
                try:
                    self.state.update(CP.PT_INPUTS, val, T)
                    self.temperature[T]['T'] += [T]
                    self.temperature[T]['p'] += [val]
                    self.temperature[T]['v'] += [1 / self.state.rhomass()]
                    self.temperature[T]['s'] += [self.state.smass()]
                    self.temperature[T]['h'] += [self.state.hmass()]
                except ValueError:
                    continue

            self.temperature[T]['T'] = np.asarray(self.temperature[T]['T'])
            self.temperature[T]['p'] = np.asarray(self.temperature[T]['p'])
            self.temperature[T]['v'] = np.asarray(self.temperature[T]['v'])
            self.temperature[T]['s'] = np.asarray(self.temperature[T]['s'])
            self.temperature[T]['h'] = np.asarray(self.temperature[T]['h'])

            if T <= self.T_crit:
                for Q in [0, 1]:
                    try:
                        self.state.update(CP.QT_INPUTS, Q, T)
                        p = self.state.p()
                        idx = np.searchsorted(self.temperature[T]['p'], p)

                        self.temperature[T]['h'] = np.insert(
                            self.temperature[T]['h'], idx, self.state.hmass())
                        self.temperature[T]['T'] = np.insert(
                            self.temperature[T]['T'], idx, T)
                        self.temperature[T]['p'] = np.insert(
                            self.temperature[T]['p'], idx, p)
                        self.temperature[T]['s'] = np.insert(
                            self.temperature[T]['s'], idx, self.state.smass())
                        self.temperature[T]['v'] = np.insert(
                            self.temperature[T]['v'], idx,
                            1 / self.state.rhomass())
                    except ValueError:
                        continue

    def isoentropy(self):
        isolines = self.entropy['isolines']

        iterator = np.geomspace(self.p_trip, self.p_max, 200)

        for s in isolines.round(8):
            self.entropy[s] = {
                'h': [], 'T': [], 'p': [], 's': [], 'v': []}
            for val in iterator:
                try:
                    self.state.update(CP.PSmass_INPUTS, val, s)
                    self.entropy[s]['T'] += [self.state.T()]
                    self.entropy[s]['p'] += [val]
                    self.entropy[s]['v'] += [1 / self.state.rhomass()]
                    self.entropy[s]['s'] += [s]
                    self.entropy[s]['h'] += [self.state.hmass()]
                except ValueError:
                    continue

            self.entropy[s]['T'] = np.asarray(self.entropy[s]['T'])
            self.entropy[s]['p'] = np.asarray(self.entropy[s]['p'])
            self.entropy[s]['v'] = np.asarray(self.entropy[s]['v'])
            self.entropy[s]['s'] = np.asarray(self.entropy[s]['s'])
            self.entropy[s]['h'] = np.asarray(self.entropy[s]['h'])

    def draw_isolines(self, diagram_type, isoline_data=None):

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

        if isoline_data is None:
            isolines = [
                x for x in self.properties.keys()
                if x not in [x_property, y_property]
            ]
        elif not isinstance(isoline_data, dict):
            msg = ('The keyword isolines must either be a dictionary or None.')
            raise TypeError(msg)
        else:
            isolines = isoline_data.keys()

        for isoline in isolines:

            property = self.properties[isoline]
            data = getattr(self, property)

            isoline_conv = self.converters[isoline][self.units[isoline]]
            x_conv = self.converters[x_property][self.units[x_property]]
            y_conv = self.converters[y_property][self.units[y_property]]

            isovalues = data['isolines']

            if isoline_data is not None:
                if 'style' in isoline_data[isoline].keys():
                    data['style'].update(isoline_data[isoline]['style'])

                if 'values' in isoline_data[isoline].keys():
                    if isoline == 'T':
                        isovalues = (
                            isoline_data[isoline]['values'] / isoline_conv[1] -
                            isoline_conv[0])
                    else:
                        isovalues = (
                            isoline_data[isoline]['values'] * isoline_conv)

            for isoval in isovalues.round(8):
                if isoval not in data.keys():
                    msg = (
                        'Could not find data for ' + property + ' isoline '
                        'with value: ' + str(isoval) + '.')
                    print(msg)
                    continue

                if x_property == 'T':
                    x = data[isoval][x_property] / x_conv[1] - x_conv[0]
                else:
                    x = data[isoval][x_property] / x_conv

                if y_property == 'T':
                    y = data[isoval][y_property] / y_conv[1] - y_conv[0]
                else:
                    y = data[isoval][y_property] / y_conv

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

                num_points = len(indices)

                if num_points > 1:
                    if isoline == 'T':
                        isoval = isoval / isoline_conv[1] - isoline_conv[0]
                    else:
                        isoval /= isoline_conv

                    if num_points < 5:
                        self.draw_isoline_label(
                            isoval.round(8), isoline, num_points - 2, x, y)
                    else:
                        self.draw_isoline_label(
                            isoval.round(8), isoline,
                            int(data['label_position'] * len(x)), x, y)
