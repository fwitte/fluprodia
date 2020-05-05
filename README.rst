=======================
Fluid Property Diagrams
=======================

Create custom Fluid Proprety Diagrams with fluprodia using CoolProp [1]_,
matplotlib [2]_ and numpy [3]_.

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

.. |docs| image:: https://readthedocs.org/projects/fluprodia/badge/?style=flat
    :target: https://readthedocs.org/projects/fluprodia
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.org/fwitte/fluprodia.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/fwitte/fluprodia

.. |coveralls| image:: https://coveralls.io/repos/fwitte/fluprodia/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/fwitte/fluprodia

.. |version| image:: https://img.shields.io/pypi/v/fluprodia.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/fluprodia

.. |wheel| image:: https://img.shields.io/pypi/wheel/fluprodia.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/fluprodia

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/fluprodia.svg
    :alt: Supported versions
    :target: https://pypi.org/project/fluprodia



.. end-badges

* Free software: MIT license

Installation and Quickstart
===========================

.. code-block:: bash

    pip install fluprodia

.. code-block:: python

    from fluprodia import FluidPropertyDiagram

    diagram = FluidPropertyDiagram(fluid='NH3')
    diagram.calc_isolines()
    diagram.set_limits(x_min=0, x_max=7000, y_min=200, y_max=450)
    diagram.draw_isolines(diagram_type='Ts')
    diagram.save('Ts_Diagramm.pdf')

Documentation
=============

See the online documentation at https://fluprodia.readthedocs.io/.

References
==========

.. [1] Bell, I., Wronski, J., Quoilin, S. and Lemort, V., 2014. Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp. *Industrial & Engineering Chemistry Research*, 53(6), pp.2498-2508.
.. [2] Hunter, J., 2007. Matplotlib: A 2D Graphics Environment. *Computing in Science & Engineering*, 9(3), pp.90-95.
.. [3] van der Walt, S., Colbert, S. and Varoquaux, G., 2011. The NumPy Array: A Structure for Efficient Numerical Computation. *Computing in Science & Engineering*, 13(2), pp.22-30.
