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

Plotting States into the Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to plot to your diagram simply by plotting on the
:code:`diagram.ax` object, which is a
:code:`matplotlib.axes._subplots.AxesSubplot` object. Therefore all matplolib
plotting functionalities are available.

For instance, if you want to plot two different states of :code:`NH3` into your
diagram, you could use the :code:`scatter()` method. If you want to have
connected states, you will need the :code:`plot()` method. In this example, we
will plot from a simple heat pump simulation in TESPy [1]_ (for more
information on TESPy see the
`online documentation <https://tespy.readthedocs.io/en/master>`_) into a logph
diagram.

.. code-block:: python

	diagram.set_limits(x_min=0, x_max=2100, y_min=1e-1, y_max=2e2)
	mydata = {
		'Q': {'values': np.linspace(0, 1, 11)},
		'T': {'values': np.arange(-50, 201, 25)}}
	diagram.draw_isolines('logph', isoline_data=mydata)

	tespy_results = run_simple_heat_pump_model()

	diagram.ax.scatter(tespy_results['h'], tespy_results['p'])
	diagram.ax.plot(tespy_results['h'], tespy_results['p'])
	diagram.save('logph_diagram_states.svg')

.. figure:: reference/_images/logph_diagram_states.svg
	:align: center

The script to generate the results is the following code snippet. Just add it
into your plotting code and it will create the results shown.

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

		cd.set_attr(pr=1, Q=-1e6)
		ev.set_attr(pr=1)
		cp.set_attr(eta_s=0.9)

		cc_cd.set_attr(fluid={'NH3': 1})
		cd_va.set_attr(Td_bp=-5, T=85)
		ev_cp.set_attr(Td_bp=5, T=15)
		nw.solve('design')

		result_dict = {
			prop: [conn.get_attr(prop).val for conn in nw.conns.index]
			for prop in ['p', 'h']
		}
		return result_dict

.. note::

    The values for plotting must be passed in the diagrams unit system.

.. [1] Witte, F., 2020. Thermal Engineering Systems in Python (Version v0.2.2). Zenodo. http://doi.org/10.5281/zenodo.3699275
