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
    * - tests
      - | |travis|
        | |coveralls|
    * - package
      - | |version| |wheel|
        | |supported-versions|
        | |zenodo|

.. |docs| image:: https://readthedocs.org/projects/fluprodia/badge/?style=flat
    :target: https://readthedocs.org/projects/fluprodia
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.org/fwitte/fluprodia.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/fwitte/fluprodia

.. |coveralls| image:: https://coveralls.io/repos/fwitte/fluprodia/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/github/fwitte/fluprodia

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

Installation and Usage
======================

.. code-block:: bash

    pip install fluprodia

.. code-block:: python

    from fluprodia import FluidPropertyDiagram

    diagram = FluidPropertyDiagram(fluid='H2O')
    diagram.set_unit_system(T='°C', h='kJ/kg', p='bar')
    diagram.calc_isolines()
    diagram.set_limits(x_min=0, x_max=8000, y_min=0, y_max=700)
    diagram.draw_isolines(diagram_type='Ts')
    diagram.save('Ts_diagram_H2O.svg')
    diagram.save('Ts_diagram_H2O.png', dpi=300)

.. figure:: docs/reference/_images/Ts_diagram_H2O.png
    :align: center

.. code-block:: python

    diagram = FluidPropertyDiagram(fluid='NH3')
    diagram.set_unit_system(T='°C', h='kJ/kg', p='bar')
    diagram.calc_isolines()
    diagram.set_limits(x_min=0, x_max=2000, y_min=1e-1, y_max=2e2)
    diagram.draw_isolines(diagram_type='logph')
    diagram.save('logph_diagram_NH3.png', dpi=300)
    diagram.save('logph_diagram_NH3.svg')

.. figure:: docs/reference/_images/logph_diagram_NH3.png
    :align: center

Documentation
=============

For further examples and usage please refer to the online documentation at
https://fluprodia.readthedocs.io/en/latest/.

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
