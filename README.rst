=======================
Fluid Property Diagrams
=======================

Create custom and beautiful Fluid Property Diagrams with fluprodia. The package
implements fluid property data from CoolProp [1]_. Plotting is handled by
matplotlib [2]_, all calculations are performed with numpy [3]_.
The list of fluids available can be found at
`CoolProp <http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids>`_.

fluprodia is licensed under the MIT software license.

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - package
      - | |version| |wheel|
        | |supported-versions|
        | |zenodo|

.. |docs| image:: https://readthedocs.org/projects/fluprodia/badge/?style=flat
    :target: https://readthedocs.org/projects/fluprodia
    :alt: Documentation Status

.. |version| image:: https://img.shields.io/pypi/v/fluprodia.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/fluprodia

.. |wheel| image:: https://img.shields.io/pypi/wheel/fluprodia.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/fluprodia

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/fluprodia.svg
    :alt: Supported versions
    :target: https://pypi.org/project/fluprodia

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3795771.svg
    :alt: Zenodo
    :target: https://doi.org/10.5281/zenodo.3795771

.. end-badges

Installation
============

.. code-block:: bash

    pip install fluprodia

Usage
=====

To create a diagram import the library, specify unit system and isolines to be
calculated and run the calculation:

.. code-block:: python

    >>> from fluprodia import FluidPropertyDiagram
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np

    >>> diagram = FluidPropertyDiagram(fluid='H2O')
    >>> diagram.set_unit_system(T='°C', h='kJ/kg', p='bar')
    >>> Q = np.linspace(0, 1, 11)
    >>> T = np.arange(25, 501, 25)
    >>> p = np.geomspace(0.01, 1000, 6) * 1e5
    >>> v = np.geomspace(0.001, 10, 5)
    >>> s = np.linspace(1000, 10000, 10)
    >>> h = np.linspace(0, 3600, 19)
    >>> diagram.set_isolines(Q=Q, T=T, p=p, v=v, s=s, h=h)
    >>> diagram.calc_isolines()

Then you can plot the data to different types of plots, e.g. logph diagram:

.. code-block:: python

    >>> fig, ax = plt.subplots(1, figsize=(8, 5))
    >>> diagram.draw_isolines(diagram_type='logph', fig=fig, ax=ax, x_min=0, x_max=3000, y_min=0.01, y_max=1000)
    >>> plt.tight_layout()
    >>> fig.savefig('logph_diagram_H2O.svg')
    >>> fig.savefig('logph_diagram_H2O.png', dpi=300)

.. figure:: https://raw.githubusercontent.com/fwitte/fluprodia/master/docs/reference/_images/logph_diagram_H2O.svg
    :align: center

Or, a Ts-diagram:

.. code-block:: python

    >>> fig, ax = plt.subplots(1, figsize=(8, 5))
    >>> diagram.draw_isolines(diagram_type='Ts', fig=fig, ax=ax, x_min=0, x_max=8000, y_min=0, y_max=700)
    >>> plt.tight_layout()
    >>> fig.savefig('Ts_diagram_H2O.svg')
    >>> fig.savefig('Ts_diagram_H2O.png', dpi=300)

.. figure:: https://raw.githubusercontent.com/fwitte/fluprodia/master/docs/reference/_images/Ts_diagram_H2O.svg
    :align: center

The fluids are available through CoolProp. To generate a diagram for a new fluid
simply change the name. Isolines come with defaults as well.

.. code-block:: python

    >>> diagram = FluidPropertyDiagram(fluid='R290')
    >>> diagram.set_unit_system(T='°C', h='kJ/kg', p='bar')
    >>> diagram.calc_isolines()
    >>> fig, ax = plt.subplots(1, figsize=(8, 5))
    >>> diagram.draw_isolines(diagram_type='logph', fig=fig, ax=ax, x_min=0, x_max=800, y_min=1e-1, y_max=1e2)
    >>> plt.tight_layout()
    >>> fig.savefig('logph_diagram_R290.png', dpi=300)
    >>> fig.savefig('logph_diagram_R290.svg')

.. figure:: https://raw.githubusercontent.com/fwitte/fluprodia/master/docs/reference/_images/logph_diagram_R290.svg
    :align: center

Documentation
=============

For further examples and usage please refer to the online documentation at
https://fluprodia.readthedocs.io/.

Citation
========

Every version of fluprodia is archived at zenodo. You can cite the latest or
a specific version. For citation info and more details please go to the
`zenodo entry <https://zenodo.org/record/3795771>`_ of fluprodia.

References
==========

This software depends on the packages CoolProp, matplolib and numpy.

.. [1] Bell, I., Wronski, J., Quoilin, S. and Lemort, V., 2014. Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp. *Industrial & Engineering Chemistry Research*, 53(6), pp. 2498-2508.
.. [2] Hunter, J., 2007. Matplotlib: A 2D Graphics Environment. *Computing in Science & Engineering*, 9(3), pp. 90-95.
.. [3] van der Walt, S., Colbert, S. and Varoquaux, G., 2011. The NumPy Array: A Structure for Efficient Numerical Computation. *Computing in Science & Engineering*, 13(2), pp. 22-30.
