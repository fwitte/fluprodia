# -*- coding: utf-8

"""Module for unit handling.

This file is part of project fluprodia (github.com/fwitte/fluprodia). It's
copyrighted by the contributors recorded in the version control history of the
file, available from its original location src/fluprodia/_units.py

SPDX-License-Identifier: MIT
"""
import shutil
import sys
import warnings

import pint
import platformdirs


# Map legacy/non-standard unit strings used in fluprodia to pint-parseable
# equivalents.  Strings already understood by pint (e.g. "Pa", "bar",
# "J/kgK" after defining kgK) are NOT listed here - they pass through as-is.
_UNIT_ALIASES = {
    "°C": "degC",
    "°F": "degF",
    "-": "1",
    "%": "percent",
    "m^3/kg": "m**3/kg",
    "l/kg": "L/kg",
}

# Map long-form / tespy-style property names to fluprodia's short internal
# keys.  Extra tespy quantities that have no fluprodia equivalent (mass_flow,
# power, ...) are simply absent - they will be silently skipped.
_LONG_KEY_MAP = {
    "pressure": "p",
    "temperature": "T",
    "volume": "vol",
    "specific_volume": "vol",    # tespy key
    "enthalpy": "h",
    "entropy": "s",
    "quality": "Q",
    "vapor_mass_fraction": "Q",  # tespy backwards-compat key
}


class Units:
    """Unit registry for fluid property diagrams.

    Wraps a :class:`pint.UnitRegistry` and stores the active unit for each
    fluid property.  The interface mirrors
    :class:`tespy.tools.units.Units` so that both libraries behave
    consistently.

    Parameters are set via :meth:`set_defaults`.

    Long-form property names (``pressure``, ``temperature``, ``volume``,
    ``enthalpy``, ``entropy``, ``quality``) are accepted alongside the short
    forms (``p``, ``T``, ``vol``, ``h``, ``s``, ``Q``).  tespy-style names
    such as ``specific_volume`` are also recognised.  Unknown keys are
    silently skipped, so a tespy :class:`tespy.tools.units.Units` object's
    ``default`` dictionary can be passed directly without pre-filtering.

    Example
    -------
    >>> from fluprodia._units import Units
    >>> u = Units()
    >>> u.default['T']
    'K'
    >>> u.set_defaults(T='°C', p='bar')
    >>> u.default['T']
    '°C'
    >>> round(u.to_SI(100.0, 'T'), 2)
    373.15
    >>> round(u.from_SI(373.15, 'T'), 2)
    100.0
    >>> u.set_defaults(temperature='K', pressure='MPa')
    >>> u.default['p']
    'MPa'
    """

    @classmethod
    def from_json(cls, data):
        """Restore a :class:`Units` instance from a serialised dictionary.

        Parameters
        ----------
        data : dict
            Dictionary as produced by :meth:`_serialize`.
        """
        instance = cls()
        instance.set_defaults(**data)
        return instance

    def __init__(self):
        # Default units - identical to the SI base units so that no
        # conversion is applied unless the user explicitly changes them.
        self.default = {
            "p": "Pa",
            "T": "K",
            "s": "J/kgK",
            "h": "J/kg",
            "vol": "m^3/kg",
            "Q": "-",
        }
        major = sys.version_info.major
        minor = sys.version_info.minor
        cache_path = platformdirs.user_cache_dir(
            "fluprodia", False, f"py{major}{minor}pint{pint.__version__}"
        )
        try:
            self._ureg = pint.UnitRegistry(cache_folder=cache_path)
        except FileNotFoundError:
            # The cache folder may point into a deleted venv; recreate it.
            shutil.rmtree(cache_path)
            self._ureg = pint.UnitRegistry(cache_folder=cache_path)
        # kgK is used in entropy units (J/kgK, kJ/kgK, ...) - same definition
        # as in tespy.tools.units
        self._ureg.define("kgK = kg * K")
        # m3 mirrors the tespy registry so unit strings like "m3/kg" work
        self._ureg.define("m3 = m ** 3")

        # Reference quantities (at SI units) used for compatibility checks
        self._quantities = {
            k: self._ureg.Quantity(1, _UNIT_ALIASES.get(v, v))
            for k, v in self.default.items()
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_defaults(self, **kwargs):
        """Set the active unit for one or more fluid properties.

        Both short-form keys (``p``, ``T``, ``vol``, ``h``, ``s``, ``Q``) and
        long-form keys (``pressure``, ``temperature``, ``volume``,
        ``enthalpy``, ``entropy``, ``quality``) are accepted.  tespy-style
        names (``specific_volume``, ``vapor_mass_fraction``) are recognised as
        well.  Keys that do not correspond to any fluprodia property are
        silently skipped, so a tespy ``Units.default`` dictionary can be
        unpacked here without pre-filtering.

        Parameters
        ----------
        p / pressure : str
            Unit of pressure, e.g. ``'Pa'``, ``'bar'``, ``'MPa'``.
        T / temperature : str
            Unit of temperature, e.g. ``'K'``, ``'°C'``, ``'°F'``.
        s / entropy : str
            Unit of specific entropy, e.g. ``'J/kgK'``, ``'kJ/kgK'``.
        h / enthalpy : str
            Unit of specific enthalpy, e.g. ``'J/kg'``, ``'kJ/kg'``.
        vol / volume / specific_volume : str
            Unit of specific volume, e.g. ``'m^3/kg'``, ``'l/kg'``.
        Q / quality : str
            Unit of vapour quality, e.g. ``'-'``, ``'%'``.

        Any unit string understood by pint and dimensionally compatible with
        the property is accepted, not just the examples listed above.
        """
        if 'v' in kwargs:
            warnings.warn(
                "The key 'v' for specific volume is deprecated and will be "
                "removed in a future release. Use 'vol' instead.",
                FutureWarning
            )
            kwargs['vol'] = kwargs.pop('v')

        # Normalise long-form / tespy-style keys to short internal keys.
        # Keys that map to nothing in fluprodia are silently skipped.
        normalised = {}
        for key, value in kwargs.items():
            short = _LONG_KEY_MAP.get(key, key)
            if short in self.default:
                normalised[short] = value
            # else: unknown key - silently skip (e.g. tespy's mass_flow, ...)

        for key, value in normalised.items():
            pint_unit = _UNIT_ALIASES.get(value, value)
            try:
                q = self._ureg.Quantity(1, pint_unit)
            except pint.UndefinedUnitError:
                raise ValueError(
                    f"The unit '{value}' is not recognised by pint. Please "
                    f"provide a valid unit string."
                )
            if not q.is_compatible_with(self._quantities[key]):
                raise ValueError(
                    f"The unit '{value}' is not compatible with the fluid "
                    f"property '{key}'."
                )
            self.default[key] = value

    def to_SI(self, value, property):
        """Convert *value* from the currently active unit to SI.

        Works element-wise on scalars and NumPy arrays alike.

        Parameters
        ----------
        value : float or numpy.ndarray
        property : str
            Short property key, e.g. ``'T'``, ``'p'``.
        """
        unit = _UNIT_ALIASES.get(self.default[property], self.default[property])
        si_unit = _UNIT_ALIASES.get(SI_UNITS[property], SI_UNITS[property])
        return self._ureg.Quantity(value, unit).to(si_unit).magnitude

    def from_SI(self, value, property):
        """Convert *value* from SI to the currently active unit.

        Works element-wise on scalars and NumPy arrays alike.

        Parameters
        ----------
        value : float or numpy.ndarray
        property : str
            Short property key, e.g. ``'T'``, ``'p'``.
        """
        unit = _UNIT_ALIASES.get(self.default[property], self.default[property])
        si_unit = _UNIT_ALIASES.get(SI_UNITS[property], SI_UNITS[property])
        return self._ureg.Quantity(value, si_unit).to(unit).magnitude

    def _serialize(self):
        """Return a plain dict suitable for JSON serialisation."""
        return dict(self.default)

    # ------------------------------------------------------------------
    # Dict-like access for backwards compatibility
    # (diagram.units['T'] still works)
    # ------------------------------------------------------------------

    def __getitem__(self, key):
        return self.default[key]

    def __setitem__(self, key, value):
        self.default[key] = value

    @property
    def ureg(self):
        """The underlying :class:`pint.UnitRegistry`."""
        return self._ureg


# Module-level singleton - same pattern as tespy's _UNITS / SI_UNITS.
# SI_UNITS is the reference for all to_SI / from_SI conversions.
_UNITS = Units()
SI_UNITS = _UNITS.default.copy()
