"""Unit tests for fluprodia._units.Units."""
import numpy as np
import pytest

from fluprodia._units import SI_UNITS, Units


class TestDefaults:
    def test_default_keys(self):
        u = Units()
        assert set(u.default) == {"p", "T", "s", "h", "vol", "Q"}

    def test_default_values_are_si(self):
        u = Units()
        assert u.default == SI_UNITS

    def test_dict_style_read(self):
        u = Units()
        assert u["T"] == "K"
        assert u["p"] == "Pa"

    def test_dict_style_write(self):
        u = Units()
        u["T"] = "°C"
        assert u.default["T"] == "°C"

    def test_si_units_module_constant(self):
        assert SI_UNITS["p"] == "Pa"
        assert SI_UNITS["T"] == "K"
        assert SI_UNITS["s"] == "J/kgK"
        assert SI_UNITS["h"] == "J/kg"
        assert SI_UNITS["vol"] == "m^3/kg"
        assert SI_UNITS["Q"] == "-"


class TestSetDefaults:
    def test_pressure_units(self):
        for unit in ("Pa", "hPa", "mbar", "kPa", "bar", "MPa", "psi"):
            u = Units()
            u.set_defaults(p=unit)
            assert u.default["p"] == unit

    def test_temperature_units(self):
        for unit in ("K", "°C", "°F"):
            u = Units()
            u.set_defaults(T=unit)
            assert u.default["T"] == unit

    def test_entropy_units(self):
        for unit in ("J/kgK", "kJ/kgK", "MJ/kgK"):
            u = Units()
            u.set_defaults(s=unit)
            assert u.default["s"] == unit

    def test_enthalpy_units(self):
        for unit in ("J/kg", "kJ/kg", "MJ/kg"):
            u = Units()
            u.set_defaults(h=unit)
            assert u.default["h"] == unit

    def test_specific_volume_units(self):
        for unit in ("m^3/kg", "l/kg"):
            u = Units()
            u.set_defaults(vol=unit)
            assert u.default["vol"] == unit

    def test_quality_units(self):
        for unit in ("-", "%"):
            u = Units()
            u.set_defaults(Q=unit)
            assert u.default["Q"] == unit

    def test_unknown_property_silently_skipped(self):
        # Unknown keys (e.g. from tespy's extra properties) are silently ignored.
        u = Units()
        u.set_defaults(mass_flow="kg/s", x_unknown_key="Pa")
        assert u.default == {"p": "Pa", "T": "K", "s": "J/kgK", "h": "J/kg", "vol": "m^3/kg", "Q": "-"}

    def test_incompatible_unit_raises(self):
        u = Units()
        with pytest.raises(ValueError, match="not compatible"):
            u.set_defaults(p="K")  # temperature unit for pressure property

    def test_incompatible_unit_raises_entropy(self):
        u = Units()
        with pytest.raises(ValueError, match="not compatible"):
            u.set_defaults(s="bar")

    def test_unknown_unit_string_raises(self):
        u = Units()
        with pytest.raises(ValueError, match="not recognised"):
            u.set_defaults(p="nonsense_unit_xyz")

    def test_v_key_deprecated(self):
        u = Units()
        with pytest.warns(FutureWarning, match="'v'.*deprecated"):
            u.set_defaults(v="l/kg")
        assert u.default["vol"] == "l/kg"

    def test_arbitrary_pint_unit_accepted(self):
        """Any dimensionally-consistent pint unit should be accepted."""
        u = Units()
        u.set_defaults(p="N/m**2")  # equivalent to Pa
        assert u.default["p"] == "N/m**2"


class TestLongFormKeys:
    """Long-form and tespy-style property names should be accepted."""

    def test_long_form_pressure(self):
        u = Units()
        u.set_defaults(pressure="bar")
        assert u.default["p"] == "bar"

    def test_long_form_temperature(self):
        u = Units()
        u.set_defaults(temperature="°C")
        assert u.default["T"] == "°C"

    def test_long_form_enthalpy(self):
        u = Units()
        u.set_defaults(enthalpy="kJ/kg")
        assert u.default["h"] == "kJ/kg"

    def test_long_form_entropy(self):
        u = Units()
        u.set_defaults(entropy="kJ/kgK")
        assert u.default["s"] == "kJ/kgK"

    def test_long_form_volume(self):
        u = Units()
        u.set_defaults(volume="l/kg")
        assert u.default["vol"] == "l/kg"

    def test_long_form_quality(self):
        u = Units()
        u.set_defaults(quality="%")
        assert u.default["Q"] == "%"

    def test_tespy_specific_volume_key(self):
        u = Units()
        u.set_defaults(specific_volume="l/kg")
        assert u.default["vol"] == "l/kg"

    def test_tespy_vapor_mass_fraction_key(self):
        u = Units()
        u.set_defaults(vapor_mass_fraction="%")
        assert u.default["Q"] == "%"

    def test_tespy_full_default_dict(self):
        """Unpacking a tespy-style defaults dict must not raise."""
        tespy_defaults = {
            "temperature": "kelvin",
            "pressure": "Pa",
            "specific_volume": "m3/kg",
            "enthalpy": "J/kg",
            "entropy": "J/kg/K",
            "quality": "1",
            # extra tespy keys that fluprodia doesn't know about:
            "mass_flow": "kg/s",
            "power": "W",
            "efficiency": "1",
        }
        u = Units()
        u.set_defaults(**tespy_defaults)
        assert u.default["p"] == "Pa"
        assert u.default["T"] == "kelvin"
        assert u.default["vol"] == "m3/kg"
        assert u.default["h"] == "J/kg"
        assert u.default["s"] == "J/kg/K"
        assert u.default["Q"] == "1"

    def test_long_form_conversion_roundtrip(self):
        u = Units()
        u.set_defaults(pressure="bar", temperature="°C", enthalpy="kJ/kg")
        assert u.to_SI(1.0, "p") == pytest.approx(1e5)
        assert u.to_SI(0.0, "T") == pytest.approx(273.15)
        assert u.to_SI(1.0, "h") == pytest.approx(1e3)


class TestToSI:
    def test_pressure_bar(self):
        u = Units()
        u.set_defaults(p="bar")
        assert u.to_SI(1.0, "p") == pytest.approx(1e5)

    def test_pressure_mpa(self):
        u = Units()
        u.set_defaults(p="MPa")
        assert u.to_SI(1.0, "p") == pytest.approx(1e6)

    def test_temperature_celsius(self):
        u = Units()
        u.set_defaults(T="°C")
        assert u.to_SI(0.0, "T") == pytest.approx(273.15)
        assert u.to_SI(100.0, "T") == pytest.approx(373.15)

    def test_temperature_fahrenheit(self):
        u = Units()
        u.set_defaults(T="°F")
        assert u.to_SI(32.0, "T") == pytest.approx(273.15, rel=1e-4)
        assert u.to_SI(212.0, "T") == pytest.approx(373.15, rel=1e-4)

    def test_entropy_kj(self):
        u = Units()
        u.set_defaults(s="kJ/kgK")
        assert u.to_SI(1.0, "s") == pytest.approx(1e3)

    def test_enthalpy_kj(self):
        u = Units()
        u.set_defaults(h="kJ/kg")
        assert u.to_SI(1.0, "h") == pytest.approx(1e3)

    def test_specific_volume_liters(self):
        u = Units()
        u.set_defaults(vol="l/kg")
        assert u.to_SI(1.0, "vol") == pytest.approx(1e-3)

    def test_quality_percent(self):
        u = Units()
        u.set_defaults(Q="%")
        assert u.to_SI(50.0, "Q") == pytest.approx(0.5)
        assert u.to_SI(100.0, "Q") == pytest.approx(1.0)

    def test_si_default_is_identity(self):
        u = Units()
        assert u.to_SI(1e5, "p") == pytest.approx(1e5)
        assert u.to_SI(300.0, "T") == pytest.approx(300.0)
        assert u.to_SI(1000.0, "s") == pytest.approx(1000.0)

    def test_numpy_array_pressure(self):
        u = Units()
        u.set_defaults(p="bar")
        result = u.to_SI(np.array([1.0, 2.0, 10.0]), "p")
        np.testing.assert_allclose(result, [1e5, 2e5, 1e6])

    def test_numpy_array_temperature(self):
        """Temperature offset conversion must work element-wise on arrays."""
        u = Units()
        u.set_defaults(T="°C")
        result = u.to_SI(np.array([0.0, 100.0]), "T")
        np.testing.assert_allclose(result, [273.15, 373.15])

    def test_numpy_array_quality(self):
        u = Units()
        u.set_defaults(Q="%")
        result = u.to_SI(np.array([0.0, 50.0, 100.0]), "Q")
        np.testing.assert_allclose(result, [0.0, 0.5, 1.0])


class TestFromSI:
    def test_pressure_bar(self):
        u = Units()
        u.set_defaults(p="bar")
        assert u.from_SI(1e5, "p") == pytest.approx(1.0)

    def test_temperature_celsius(self):
        u = Units()
        u.set_defaults(T="°C")
        assert u.from_SI(273.15, "T") == pytest.approx(0.0, abs=1e-9)
        assert u.from_SI(373.15, "T") == pytest.approx(100.0)

    def test_specific_volume_liters(self):
        u = Units()
        u.set_defaults(vol="l/kg")
        assert u.from_SI(1e-3, "vol") == pytest.approx(1.0)

    def test_quality_percent(self):
        u = Units()
        u.set_defaults(Q="%")
        assert u.from_SI(0.5, "Q") == pytest.approx(50.0)

    def test_numpy_array_enthalpy(self):
        u = Units()
        u.set_defaults(h="kJ/kg")
        result = u.from_SI(np.array([0.0, 1e3, 2e3]), "h")
        np.testing.assert_allclose(result, [0.0, 1.0, 2.0])

    def test_numpy_array_temperature(self):
        u = Units()
        u.set_defaults(T="°C")
        result = u.from_SI(np.array([273.15, 373.15]), "T")
        np.testing.assert_allclose(result, [0.0, 100.0])

    def test_roundtrip_scalar(self):
        u = Units()
        u.set_defaults(T="°C", p="bar", s="kJ/kgK", h="kJ/kg", vol="l/kg", Q="%")
        cases = [
            ("T", 25.0),
            ("p", 3.5),
            ("s", 7.2),
            ("h", 1200.0),
            ("vol", 0.5),
            ("Q", 80.0),
        ]
        for prop, val in cases:
            assert u.from_SI(u.to_SI(val, prop), prop) == pytest.approx(val, rel=1e-9)

    def test_roundtrip_array(self):
        u = Units()
        u.set_defaults(T="°C", p="bar")
        arr_T = np.array([0.0, 25.0, 100.0])
        np.testing.assert_allclose(u.from_SI(u.to_SI(arr_T, "T"), "T"), arr_T)
        arr_p = np.array([1.0, 5.0, 10.0])
        np.testing.assert_allclose(u.from_SI(u.to_SI(arr_p, "p"), "p"), arr_p)


class TestSerializeFromJson:
    def test_serialize_returns_dict_copy(self):
        u = Units()
        d = u._serialize()
        assert isinstance(d, dict)
        d["p"] = "bar"
        assert u.default["p"] == "Pa"  # original not mutated

    def test_serialize_default(self):
        u = Units()
        assert u._serialize() == SI_UNITS

    def test_from_json_default(self):
        u = Units.from_json(SI_UNITS)
        assert u.default == SI_UNITS

    def test_from_json_custom_units(self):
        u = Units()
        u.set_defaults(T="°C", p="bar", h="kJ/kg")
        restored = Units.from_json(u._serialize())
        assert restored.default == u.default

    def test_from_json_converts_correctly(self):
        u = Units.from_json({"p": "bar", "T": "°C", "s": "J/kgK",
                             "h": "J/kg", "vol": "m^3/kg", "Q": "-"})
        assert u.to_SI(1.0, "p") == pytest.approx(1e5)
        assert u.to_SI(0.0, "T") == pytest.approx(273.15)

    def test_from_json_v_key_backwards_compat(self):
        """JSON files written before the vol rename should still load."""
        data = {"p": "Pa", "T": "K", "s": "J/kgK",
                "h": "J/kg", "v": "m^3/kg", "Q": "-"}
        with pytest.warns(FutureWarning):
            u = Units.from_json(data)
        assert u.default["vol"] == "m^3/kg"
