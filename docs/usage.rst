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

    iso_T = np.arange(-25, 151, 25)
    diagram.set_isolines(T=iso_T)
    diagram.calc_isolines()
    diagram.set_limits(x_min=0, x_max=8000, y_min=0, y_max=700)
    diagram.draw_isolines('Ts')
    diagram.save('Ts_diagram.pdf')

.. figure:: reference/_images/Ts_diagram.svg
    :align: center

As all fluid properties will be stored in the object referenced by
:code:`diagram`, it is possible to change the diagram type and export a new
diagram without recalculating the isolines. Only if you wish to draw a
different set of isolines unlike specified in the :code:`set_isolines()` method
call, you need to recalculate the isolines.

.. code-block:: python

    diagram.set_limits(x_min=0, x_max=2000, y_min=1e-1, y_max=1e2)
    diagram.draw_isolines('logph')
    diagram.save('logph_diagram.pdf')

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

    T = np.arange(-25, 151, 5)
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

    diagram.set_limits(x_min=0, x_max=2000, y_min=1e-1, y_max=1e2)
    mydata = {
        'Q': {'values': np.linspace(0, 1, 11)},
        'T': {'values': np.arange(-25, 151, 25)}}
    diagram.drawing('logph', isoline_data=mydata)
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

    diagram.set_limits(x_min=800, x_max=1600, y_min=1, y_max=75)
    mydata = {
        'T': {'style': {'color': '#ff0000'}},
        'v': {'values': []}
	}
    diagram.drawing('logph', isoline_data=mydata)
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
labels for the isolines. This placing of the labels is not fully satisfactory
at the moment. If you have ideas, how to place the labels in an improved way,
we are looking forward for you suggestions.

You can specify a positioning value between 0 and 1. Every label of an
isoline type (e.g. constant enthalpy) will be placed at the relative position
of each isoline within the limits of the view.

.. code-block:: python
	mydata = {
		'T': {'label_position': 0.5}
	}
	diagram.drawing('logph', isoline_data=mydata)
	diagram.save('logph_NH3_zoomed_temperature_labels.svg')

.. figure:: reference/_images/logph_NH3_zoomed_temperature_labels.svg
    :align: center

Plotting States into the Diagram
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to plot to your diagram simply by plotting on the
:code:`diagram.ax` object, which is a
:code:`matplotlib.axes._subplots.AxesSubplot` object. Therefore all matplolib
plotting functionalities are available.

For instance, if you want to plot two different states of :code:`NH3` into your
diagram, you could use the :code:`scatter()` method. If you want to have
connected states, you will need the :code:`plot()` method. In this example, we
will plot into a logph diagram.

.. code-block:: python

    diagram.set_limits(x_min=0, x_max=2000, y_min=1e-1, y_max=1e2)
    diagram.draw_isolines('logph')
    p = [1e1, 1e1]
    h = [700, 1600]

    diagram.ax.scatter(h, p)
    diagram.ax.plot(h, p)
    diagram.save('logph_diagram.pdf')

.. note::

    You need to make sure to pass your values in the diagrams units.

.. figure:: reference/_images/logph_diagram_NH3_plot_states.svg
    :align: center
