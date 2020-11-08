=====
Usage
=====

Introduction
^^^^^^^^^^^^

After installation of fluprodia you can easily create fluid property diagrams
for all pure and pseudo-pure fluids from the CoolProp fluid property database.
For a list of available fluids please refer to the online documentation of
`CoolProp <http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids>`_.

In order to start, import the package and create an object of the class
:py:class:`fluprodia.fluid_property_diagram.FluidPropertyDiagram` by passing
the alias of the fluid. After that, it is possible to specify a unit system
for all fluid properties available with the
:py:meth:`fluprodia.fluid_property_diagram.FluidPropertyDiagram.set_unit_system`
method. The fluid properties available are:

- pressure :code:`p`
- specific enthalpy :code:`h`
- specific entropy :code:`s`
- specific volume :code:`v`
- temperature :code:`T`
- vapor mass fraction :code:`Q`

.. code-block:: python

    from fluprodia import FluidPropertyDiagram
    import numpy as np
    diagram = FluidPropertyDiagram('NH3')
    diagram.set_unit_system(T='°C', p='bar', h='kJ/kg')

After that, you can use the default isolines or specify your own lines by
using the
:py:meth:`fluprodia.fluid_property_diagram.FluidPropertyDiagram.set_isolines`
method. If you do not specify custom isolines, generic isolines will be used
instead. Next step is to calculate the isolines, drawing them and exporting the
diagram in your favorite format. The formats available are the matplotlib file
formats for figures. You will also need to specify the limits in order to
determine the view. Also, different diagrams will have different value ranges
for their x- and y-axes.

.. code-block:: python

    iso_T = np.arange(-50, 201, 25)
    diagram.set_isolines(T=iso_T)
    diagram.calc_isolines()
    diagram.set_limits(x_min=0, x_max=8000, y_min=-50, y_max=200)
    diagram.draw_isolines('Ts')
    diagram.save('Ts_diagram.svg')

.. figure:: reference/_images/Ts_diagram.svg
    :align: center

As all fluid properties will be stored in the object referenced by
:code:`diagram`, it is possible to change the diagram type and export a new
diagram without recalculating the isolines. Only if you wish to draw a
different set of isolines unlike specified in the :code:`set_isolines()` method
call, you need to recalculate the isolines.

.. code-block:: python

    diagram.set_limits(x_min=0, x_max=2100, y_min=1e-1, y_max=2e2)
    diagram.draw_isolines('logph')
    diagram.save('logph_diagram.svg')

.. figure:: reference/_images/logph_diagram.svg
    :align: center

All available diagram types can be displayed by printing the following line.

.. code-block:: python

    print(diagram.supported_diagrams.keys())

Customizing the Display
^^^^^^^^^^^^^^^^^^^^^^^

Customization is possible regarding

- the isovalues of the isolines,
- the isolines displayed,
- the linestyle of the isolines and
- the position of the isolines' labels.

Isoline values available
************************

As already mentioned, you can set the isolines for your diagram like this. All
isolines you specify are available for drawing the diagram later. Therefore,
the more values you specify, the more lines can be displayed. Also, the
computation time will rise.

Still, it might be useful to specify a lot of values. E.g., if we want to
create a full view of a logph diagram for NH3 and a zoomed view in the two
phase region with lines of constant vapor mass fraction for every 2.5 % and
lines of constant temperature every 5 K.

.. code-block:: python

    T = np.arange(-50, 201, 5)
    Q = np.linspace(0, 1, 41)
    diagram.set_isolines(T=T, Q=Q)
    diagram.calc_isolines()

The following sections shows how to select from all isolines available.

Lines displayed and Linestyle
*****************************

As we do not want to display all values for temperature and vapor mass fraction
for the full view diagram, we specify the values to be displayed for these
properties. This is done by using the isoline_data property, which must be
a dictionary holding the required information.

.. code-block:: python

    diagram.set_limits(x_min=0, x_max=2100, y_min=1e-1, y_max=2e2)
    mydata = {
        'Q': {'values': np.linspace(0, 1, 11)},
        'T': {'values': np.arange(-50, 201, 25)}}
    diagram.draw_isolines('logph', isoline_data=mydata)
    diagram.save('logph_NH3_full.svg')

.. figure:: reference/_images/logph_NH3_full.svg
    :align: center

Now, for the zoomed diagram we want the full temperature and vapor mass
fraction data. At the same time, you might want to change the color or the
linestyle of an isoline. For this example, we will color the lines of constant
temperature in red. Additionally, the lines of constant specific volume should
not be displayed at all. This can be done by passing an empty list or an empty
numpy array.

.. code-block:: python

    diagram.set_limits(x_min=1000, x_max=1500, y_min=1, y_max=2e2)
    mydata = {
        'T': {
            'style': {'color': '#ff0000'},
            'values': np.arange(-50, 201, 5)},
        'v': {'values': np.array([])}
        }
    diagram.draw_isolines('logph', isoline_data=mydata)
    diagram.save('logph_NH3_zoomed.svg')

.. figure:: reference/_images/logph_NH3_zoomed.svg
    :align: center

.. note::

    For changing the style of a specific isoline pass the respective keyword
    and value pairs in a dictionary. The keywords available are the keywords
    of a :code:`matplotlib.lines.Line2D` object. See
    https://matplotlib.org/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
    for more information.

Positioning of the isoline lables
*********************************

In the last section we briefly describe, how to change the placing of the
labels for the isolines. Looking at the zoomed diagram, you see that some of
the temperature labels are missing.

You can specify a positioning value between 0 and 1. Every label of an
isoline type (e.g. constant temerature) will be placed at the relative position
of each isoline within the limits of the view.

.. code-block:: python

    diagram.set_limits(x_min=1000, x_max=1500, y_min=1, y_max=2e2)
    mydata = {
        'T': {
            'style': {'color': '#ff0000'},
            'values': np.arange(-50, 201, 5),
            'label_position': 0.8},
        'v': {'values': np.array([])}
        }
    diagram.draw_isolines('logph', isoline_data=mydata)
    diagram.save('logph_NH3_zoomed_temperature_labels.svg')

.. figure:: reference/_images/logph_NH3_zoomed_temperature_labels.svg
    :align: center

.. note::

    The placing method of the labels is not fully satisfactory at the moment.
    If you have ideas, how to place the labels in an improved way, we are
    looking forward for you suggestions.

Plotting individual isolines (and isolike lines)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FluProDia offers a method to generate data for individual isolines with a
specified starting and a specified ending point. Use the method
:py:meth:`fluprodia.fluid_property_diagram.FluidPropertyDiagram.calc_individual_isoline`
to create datapoints for the isoline. The method returns a dictionary
containing the datapoints in numpy arrays using the property name as
respective key. Therefore, independent of the diagram you want to draw, you
will have all data available. Following, we will draw all available isolines
into a Ts and a logph diagram. Each property value must be passed in the
diagram's respective unit system.

.. code-block:: python

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

With these data, it is possible to plot to your diagram simply by plotting on
the :code:`diagram.ax` object, which is a
:code:`matplotlib.axes._subplots.AxesSubplot` object. Therefore all matplolib
plotting functionalities are available. Simply pass the data of the x and y
property of your diagram, e.g. to the :code:`plot()` method.

.. code-block:: python

    diagram.set_limits(x_min=0, x_max=2100, y_min=1e-1, y_max=2e2)
    mydata = {
        'Q': {'values': np.linspace(0, 1, 11)},
        'T': {'values': np.arange(-50, 201, 25)}}
    diagram.draw_isolines('logph', isoline_data=mydata)
    for key, specs in data.items():
        datapoints = specs['datapoints']
        diagram.ax.plot(specs['datapoints']['h'], specs['datapoints']['p'], label=key)
    diagram.ax.legend(loc='lower right')
    diagram.save('logph_NH3_isolines.svg')

    diagram.set_limits(x_min=0, x_max=7000, y_min=-50, y_max=201)
    diagram.draw_isolines('Ts')
    for key, specs in data.items():
        datapoints = specs['datapoints']
        diagram.ax.plot(specs['datapoints']['s'], specs['datapoints']['T'], label=key)
    diagram.ax.legend(loc='lower right')
    diagram.save('Ts_NH3_isolines.svg')

.. figure:: reference/_images/logph_NH3_isolines.svg
    :align: center

.. figure:: reference/_images/Ts_NH3_isolines.svg
    :align: center

.. note::

    Note that the :code:`starting_point_property` and the
    :code:`ending_point_property` do not need to be identical! E.g., you can
    draw an isobaric line starting at a specific entropy and ending at a
    specific temperature.

On top of that, e.g. in order to display a pressure loss in a heat exchanger,
you can have different values for the (iso)line at the starting and the ending
points. The (then former) isoline property will be changed linearly to either
change in entropy (for isobars and isothermals) or change in pressure (for all
other lines). This functionality is only supposed to display the change in a
beautiful way, it does not represent the actual process connecting your
starting point with your ending point as this would require perfect knowledge
of the process. In order to generate these data, you need to pass the
:code:`'isoline_value_end'` keyword to the
:py:meth:`fluprodia.fluid_property_diagram.FluidPropertyDiagram.calc_individual_isoline`
method.

.. code-block:: python

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
    diagram.save('Ts_NH3_pressure_loss.svg')

.. figure:: reference/_images/Ts_NH3_pressure_loss.svg
    :align: center

Plotting States into the Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For instance, if you want to plot two different states of :code:`NH3` into your
diagram, you could use the :code:`scatter()` method. If you want to have
connected states, you will need the :code:`plot()` method. In this example, we
will plot from a simple heat pump simulation in TESPy [1]_ (for more
information on TESPy see the
`online documentation <https://tespy.readthedocs.io/en/master>`_) into a logph
and a Ts diagram.

.. code-block:: python

    tespy_results = run_simple_heat_pump_model()
    for key, data in tespy_results.items():
        tespy_results[key]['datapoints'] = diagram.calc_individual_isoline(**data)

    diagram.set_limits(x_min=0, x_max=2100, y_min=1e0, y_max=2e2)
    mydata = {
        'Q': {'values': np.linspace(0, 1, 11)},
        'T': {
            'values': np.arange(-25, 226, 25),
            'style': {'color': '#000000'}
        }
    }
    diagram.draw_isolines('logph', isoline_data=mydata)

    for key in tespy_results.keys():
        datapoints = tespy_results[key]['datapoints']
        diagram.ax.plot(datapoints['h'], datapoints['p'], color='#ff0000')
        diagram.ax.scatter(datapoints['h'][0], datapoints['p'][0], color='#ff0000')
    diagram.save('logph_diagram_states.svg')

    diagram.set_limits(x_min=2000, x_max=7000, y_min=-50, y_max=225)
    diagram.draw_isolines('Ts')

    for key in tespy_results.keys():
        datapoints = tespy_results[key]['datapoints']
        diagram.ax.plot(datapoints['s'], datapoints['T'], color='#ff0000')
        diagram.ax.scatter(datapoints['s'][0], datapoints['T'][0], color='#ff0000')
    diagram.save('Ts_diagram_states.svg')

.. figure:: reference/_images/logph_diagram_states.svg
    :align: center

.. figure:: reference/_images/Ts_diagram_states.svg
    :align: center

The script to generate the results is the following code snippet. Just add it
into your plotting code and it will create the results shown. An interface
automatically generating a dictionary for every component of the network is
planned in future versions of TESPy.

.. code-block:: python

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

.. note::

    The values for plotting must be passed in the diagrams unit system.

.. [1] Witte, F., 2020. Thermal Engineering Systems in Python (Version v0.3.4). Zenodo. http://doi.org/10.5281/zenodo.3699275
