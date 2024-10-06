
Changelog
=========

v3.4 (October, 6, 2024)
-----------------------

* Refactor the generators for the isolines to make them more resilient and
  distribute points better within the given range of each isoline.
* Add a method that allows the creation of isolines only within a specific
  range of the fluid properties.

v3.3 (July, 8, 2024)
--------------------

* Make fluprodia compatible with numpy version 2.0.

v3.2 (June, 30, 2024)
---------------------

* A class method to construct a `FluidPropertyDiagram` from the `to_json` data
  export is available.

v3.1 (May, 23, 2024)
--------------------

* The `draw_isolines` method is now compatible with matplotlib darkmode style.

v3.0 (April, 26, 2024)
----------------------

* You can export the underlying data of your diagram using the `to_json` method.

v2.0 (November, 22, 2023)
-------------------------

* The API changed in a way, that you have to create the figure and axes
  externally using matplotlib. These are then passed to the `draw_isolines`
  method. See the README or in the examples for the new API version.

v1.6 (December, 02, 2022)
-------------------------

* Remove upper Python version limit.

v1.5 (July, 28, 2021)
---------------------

* Update documentation on pressure units.
* Improve error message for not available units.

v1.4 (July, 28, 2021)
---------------------

* Add kPa to pressure unit system.
* Fix TESPy API calls in the example usage.

v1.3 (January, 7, 2021)
-----------------------

* Reduce the number of datapoints for isolines to 200 for faster performance.

v1.2 (December, 8, 2020)
------------------------

* Fix minimum volume value for iterators.

v1.1 (November, 10, 2020)
-------------------------

* Change the iterator for isobaric, isenthalpic and isentropic to specific volume.
* Adjust individual isoline plotting iterators and isolines accordingly.

v1.0 (November, 8, 2020)
------------------------

* Add method to calculate datapoints of individual isolines and isolike lines.

v0.1.2 (October, 2, 2020)
-------------------------

* Minor bug fixes for isochoric drawing.
* Change in default values for isobarics.

v0.1.1 (May, 13, 2020)
----------------------

* Catch exceptions in calculation of minimum specific volume for default
  isoline generation.
* Allow Python 3.8 usage.

v0.1.0 (May, 6, 2020)
---------------------

* First release on PyPI.
