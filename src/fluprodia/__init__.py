# -*- coding: utf-8
"""Create beautiful fluid property diagrams using CoolProp"""

import sys
import warnings

__version__ = '4.2'

if sys.version_info < (3, 11):
    warnings.warn(
        f"Python {sys.version_info.major}.{sys.version_info.minor} is no "
        f"longer supported as of the next major release of fluprodia. "
        f"Please upgrade to Python 3.11 or newer.",
        FutureWarning
    )

from .fluid_property_diagram import FluidPropertyDiagram  # noqa: F401
