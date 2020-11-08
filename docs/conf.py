# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os
import traceback
import fluprodia
from fluprodia import FluidPropertyDiagram
import numpy as np
from tespy.components import (
    compressor, cycle_closer, heat_exchanger_simple, valve)
from tespy.connections import connection
from tespy.networks import network


def run_simple_heat_pump_model():
    nw = network(['NH3'], T_unit='C', p_unit='bar', h_unit='kJ / kg')
    nw.set_attr(iterinfo=False)
    cp = compressor('compressor')
    cc = cycle_closer('cycle_closer')
    cd = heat_exchanger_simple('condenser')
    va = valve('expansion valve')
    ev = heat_exchanger_simple('evaporator')

    cp.char_warnings = False
    cd.char_warnings = False
    va.char_warnings = False
    ev.char_warnings = False

    cc_cd = connection(cc, 'out1', cd, 'in1')
    cd_va = connection(cd, 'out1', va, 'in1')
    va_ev = connection(va, 'out1', ev, 'in1')
    ev_cp = connection(ev, 'out1', cp, 'in1')
    cp_cc = connection(cp, 'out1', cc, 'in1')

    nw.add_conns(cc_cd, cd_va, va_ev, ev_cp, cp_cc)

    cd.set_attr(pr=0.95, Q=-1e6)
    ev.set_attr(pr=0.9)
    cp.set_attr(eta_s=0.9)

    cc_cd.set_attr(fluid={'NH3': 1})
    cd_va.set_attr(Td_bp=-5, T=85)
    ev_cp.set_attr(Td_bp=5, T=15)
    nw.solve('design')

    result_dict = {
        'compressor': {
            'isoline_property': 's',
            'isoline_value': ev_cp.s.val,
            'isoline_value_end': cp_cc.s.val,
            'starting_point_property': 'p',
            'starting_point_value': ev_cp.p.val,
            'ending_point_property': 'p',
            'ending_point_value': cp_cc.p.val
        },
        'condenser': {
            'isoline_property': 'p',
            'isoline_value': cc_cd.p.val,
            'isoline_value_end': cd_va.p.val,
            'starting_point_property': 's',
            'starting_point_value': cc_cd.s.val,
            'ending_point_property': 's',
            'ending_point_value': cd_va.s.val
        },
        'valve': {
            'isoline_property': 'h',
            'isoline_value': cd_va.h.val,
            'isoline_value_end': va_ev.h.val,
            'starting_point_property': 'p',
            'starting_point_value': cd_va.p.val,
            'ending_point_property': 'p',
            'ending_point_value': va_ev.p.val
        },
        'evaporator': {
            'isoline_property': 'p',
            'isoline_value': va_ev.p.val,
            'isoline_value_end': ev_cp.p.val,
            'starting_point_property': 's',
            'starting_point_value': va_ev.s.val,
            'ending_point_property': 's',
            'ending_point_value': ev_cp.s.val
        }
    }
    return result_dict


diagram = FluidPropertyDiagram('NH3')
diagram.set_unit_system(T='°C', p='bar', h='kJ/kg')
iso_T = np.arange(-50, 201, 25)
diagram.set_isolines(T=iso_T)
diagram.calc_isolines()
diagram.set_limits(x_min=0, x_max=8000, y_min=-50, y_max=200)
diagram.draw_isolines('Ts')
diagram.save('reference/_images/Ts_diagram.svg')

diagram.set_limits(x_min=0, x_max=2100, y_min=1e-1, y_max=2e2)
diagram.draw_isolines('logph')
diagram.save('reference/_images/logph_diagram.svg')

print(diagram.supported_diagrams.keys())

T = np.arange(-50, 201, 5)
Q = np.linspace(0, 1, 41)
diagram.set_isolines(T=T, Q=Q)
diagram.calc_isolines()

diagram.set_limits(x_min=0, x_max=2100, y_min=1e-1, y_max=2e2)
mydata = {
    'Q': {'values': np.linspace(0, 1, 11)},
    'T': {'values': np.arange(-50, 201, 25)}}
diagram.draw_isolines('logph', isoline_data=mydata)
diagram.save('reference/_images/logph_NH3_full.svg')

diagram.set_limits(x_min=1000, x_max=1500, y_min=1, y_max=2e2)
mydata = {
    'T': {
        'style': {'color': '#ff0000'},
        'values': np.arange(-50, 201, 5)},
    'v': {'values': np.array([])}
    }
diagram.draw_isolines('logph', isoline_data=mydata)
diagram.save('reference/_images/logph_NH3_zoomed.svg')

diagram.set_limits(x_min=1000, x_max=1500, y_min=1, y_max=2e2)
mydata = {
    'T': {
        'style': {'color': '#ff0000'},
        'values': np.arange(-50, 201, 5),
        'label_position': 0.8},
    'v': {'values': np.array([])}
    }
diagram.draw_isolines('logph', isoline_data=mydata)
diagram.save('reference/_images/logph_NH3_zoomed_temperature_labels.svg')

diagram.set_limits(x_min=0, x_max=2100, y_min=1e0, y_max=2e2)
mydata = {
    'Q': {'values': np.linspace(0, 1, 11)},
    'T': {
        'values': np.arange(-25, 226, 25),
        'style': {'color': '#000000'}
    }
}
diagram.draw_isolines('logph', isoline_data=mydata)

tespy_results = run_simple_heat_pump_model()
for key, data in tespy_results.items():
    tespy_results[key]['datapoints'] = diagram.calc_individual_isoline(**data)

for key in tespy_results.keys():
    datapoints = tespy_results[key]['datapoints']
    diagram.ax.plot(datapoints['h'], datapoints['p'], color='#ff0000')
    diagram.ax.scatter(datapoints['h'][0], datapoints['p'][0], color='#ff0000')
diagram.save('reference/_images/logph_diagram_states.svg')

diagram.set_limits(x_min=2000, x_max=7000, y_min=-50, y_max=225)
diagram.draw_isolines('Ts')

for key in tespy_results.keys():
    datapoints = tespy_results[key]['datapoints']
    diagram.ax.plot(datapoints['s'], datapoints['T'], color='#ff0000')
    diagram.ax.scatter(datapoints['s'][0], datapoints['T'][0], color='#ff0000')
diagram.save('reference/_images/Ts_diagram_states.svg')


data = {
    'isobaric': {
        'isoline_property': 'p',
        'isoline_value': 10,
        'starting_point_property': 'T',
        'starting_point_value': -25,
        'ending_point_property': 'T',
        'ending_point_value': 150
    },
    'isochoric': {
        'isoline_property': 'v',
        'isoline_value': 0.035,
        'starting_point_property': 'h',
        'starting_point_value': 750,
        'ending_point_property': 'T',
        'ending_point_value': 150
    },
    'isothermal': {
        'isoline_property': 'T',
        'isoline_value': 65,
        'starting_point_property': 'v',
        'starting_point_value': 0.01,
        'ending_point_property': 'v',
        'ending_point_value': 0.5
    },
    'isenthalpic': {
        'isoline_property': 'h',
        'isoline_value': 850,
        'starting_point_property': 'p',
        'starting_point_value': 200,
        'ending_point_property': 'v',
        'ending_point_value': 0.5
    },
    'isentropic': {
        'isoline_property': 's',
        'isoline_value': 4700,
        'starting_point_property': 'T',
        'starting_point_value': -20,
        'ending_point_property': 'T',
        'ending_point_value': 150
    }
}

for name, specs in data.items():
    data[name]['datapoints'] = diagram.calc_individual_isoline(**specs)

diagram.set_limits(x_min=0, x_max=2100, y_min=1e-1, y_max=2e2)
mydata = {
    'Q': {'values': np.linspace(0, 1, 11)},
    'T': {'values': np.arange(-50, 201, 25)}}
diagram.draw_isolines('logph', isoline_data=mydata)
for key, specs in data.items():
    datapoints = specs['datapoints']
    diagram.ax.plot(specs['datapoints']['h'], specs['datapoints']['p'], label=key)
diagram.ax.legend(loc='lower right')
diagram.save('reference/_images/logph_NH3_isolines.svg')

diagram.set_limits(x_min=0, x_max=7000, y_min=-50, y_max=201)
diagram.draw_isolines('Ts')
for key, specs in data.items():
    datapoints = specs['datapoints']
    diagram.ax.plot(specs['datapoints']['s'], specs['datapoints']['T'], label=key)
diagram.ax.legend(loc='lower right')
diagram.save('reference/_images/Ts_NH3_isolines.svg')

data = {
    'isoline_property': 'p',
    'isoline_value': 10,
    'isoline_value_end': 9,
    'starting_point_property': 'h',
    'starting_point_value': 350,
    'ending_point_property': 'h',
    'ending_point_value': 1750
}
datapoints = diagram.calc_individual_isoline(**data)
diagram.draw_isolines('Ts')
for specs in data.values():
    diagram.ax.plot(datapoints['s'], datapoints['T'])
diagram.save('reference/_images/Ts_NH3_pressure_loss.svg')


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
]
source_suffix = '.rst'
master_doc = 'index'
project = 'fluprodia'
year = '2020'
author = 'Francesco Witte'
copyright = '{0}, {1}'.format(year, author)

version = release = fluprodia.__version__

pygments_style = 'trac'
templates_path = ['.']
extlinks = {
    'issue': ('https://github.com/fwitte/fluprodia/issues/%s', '#'),
    'pr': ('https://github.com/fwitte/fluprodia/pull/%s', 'PR #'),
}
# on_rtd is whether we are on readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only set the theme if we're building docs locally
    html_theme = 'sphinx_rtd_theme'

html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
   '**': ['searchbox.html', 'globaltoc.html', 'sourcelink.html'],
}
html_short_title = '%s-%s' % (project, version)

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False
