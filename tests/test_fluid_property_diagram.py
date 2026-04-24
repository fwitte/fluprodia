"""Unit tests for FluidPropertyDiagram."""
import json
import warnings

import numpy as np
import pytest

from fluprodia import FluidPropertyDiagram
from fluprodia._units import Units

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def diagram():
    """Diagram for water without computed isolines (fast to construct)."""
    return FluidPropertyDiagram("water")


@pytest.fixture(scope="module")
def diagram_r290():
    """R290 diagram with computed isolines (shared across tests in module)."""
    d = FluidPropertyDiagram("R290")
    d.calc_isolines()
    return d


# ---------------------------------------------------------------------------
# Unit system
# ---------------------------------------------------------------------------

class TestUnitSystemIntegration:
    def test_units_is_units_instance(self, diagram):
        assert isinstance(diagram.units, Units)

    def test_set_unit_system_updates_units(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(T="°C", p="bar", s="kJ/kgK", h="kJ/kg")
        assert d.units["T"] == "°C"
        assert d.units["p"] == "bar"
        assert d.units["s"] == "kJ/kgK"
        assert d.units["h"] == "kJ/kg"

    def test_set_unit_system_v_key_deprecated(self):
        d = FluidPropertyDiagram("water")
        with pytest.warns(FutureWarning, match="'v'.*deprecated"):
            d.set_unit_system(v="l/kg")
        assert d.units["vol"] == "l/kg"

    def test_set_unit_system_incompatible_raises(self):
        d = FluidPropertyDiagram("water")
        with pytest.raises(ValueError, match="not compatible"):
            d.set_unit_system(p="K")

    def test_set_unit_system_long_form_keys(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(pressure="bar", temperature="°C", enthalpy="kJ/kg", entropy="kJ/kgK")
        assert d.units["p"] == "bar"
        assert d.units["T"] == "°C"
        assert d.units["h"] == "kJ/kg"
        assert d.units["s"] == "kJ/kgK"

    def test_set_unit_system_tespy_units_object(self):
        """A tespy-style Units object (with .default dict) can be passed directly."""
        class FakeTespyUnits:
            default = {
                "temperature": "degC",
                "pressure": "bar",
                "specific_volume": "m3/kg",
                "enthalpy": "kJ/kg",
                "entropy": "J/kg/K",
                "quality": "1",
                "mass_flow": "kg/s",
                "power": "W",
            }

        d = FluidPropertyDiagram("water")
        d.set_unit_system(FakeTespyUnits())
        assert d.units["p"] == "bar"
        assert d.units["T"] == "degC"
        assert d.units["h"] == "kJ/kg"

    def test_set_unit_system_tespy_units_kwargs_override(self):
        """Keyword arguments override the tespy Units values."""
        class FakeTespyUnits:
            default = {"pressure": "bar", "temperature": "kelvin"}

        d = FluidPropertyDiagram("water")
        d.set_unit_system(FakeTespyUnits(), pressure="MPa")
        assert d.units["p"] == "MPa"
        assert d.units["T"] == "kelvin"

    def test_set_unit_system_invalid_units_arg_raises(self):
        d = FluidPropertyDiagram("water")
        with pytest.raises(TypeError, match="default"):
            d.set_unit_system("not_a_units_object")

    def test_units_dict_access(self, diagram):
        assert diagram.units["T"] == "K"
        assert diagram.units["p"] == "Pa"


class TestConvertToFromSI:
    def test_convert_to_si_pressure(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(p="bar")
        assert d.convert_to_SI(1.0, "p") == pytest.approx(1e5)

    def test_convert_from_si_pressure(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(p="bar")
        assert d.convert_from_SI(1e5, "p") == pytest.approx(1.0)

    def test_convert_to_si_temperature_celsius(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(T="°C")
        assert d.convert_to_SI(0.0, "T") == pytest.approx(273.15)
        assert d.convert_to_SI(100.0, "T") == pytest.approx(373.15)

    def test_convert_from_si_temperature_celsius(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(T="°C")
        assert d.convert_from_SI(373.15, "T") == pytest.approx(100.0)

    def test_convert_to_si_array(self):
        """Vectorised conversion must work for numpy arrays."""
        d = FluidPropertyDiagram("water")
        d.set_unit_system(T="°C")
        result = d.convert_to_SI(np.array([0.0, 100.0]), "T")
        np.testing.assert_allclose(result, [273.15, 373.15])

    def test_convert_from_si_array(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(p="bar")
        result = d.convert_from_SI(np.array([1e5, 2e5]), "p")
        np.testing.assert_allclose(result, [1.0, 2.0])

    def test_roundtrip_scalar(self):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(T="°C", p="bar", s="kJ/kgK", h="kJ/kg")
        cases = [("T", 25.0), ("p", 5.0), ("s", 6.5), ("h", 800.0)]
        for prop, val in cases:
            assert d.convert_from_SI(d.convert_to_SI(val, prop), prop) == pytest.approx(val, rel=1e-9)


# ---------------------------------------------------------------------------
# vol / v backwards compatibility
# ---------------------------------------------------------------------------

class TestVolKey:
    def test_set_isolines_vol(self):
        d = FluidPropertyDiagram("water")
        import numpy as np
        d.set_isolines(vol=np.array([1e-3, 1e-2, 1e-1]))
        assert len(d.volume["isolines"]) == 3

    def test_set_isolines_v_deprecated(self):
        d = FluidPropertyDiagram("water")
        import numpy as np
        with pytest.warns(FutureWarning, match="'v'.*deprecated"):
            d.set_isolines(v=np.array([1e-3, 1e-2]))
        assert len(d.volume["isolines"]) == 2

    def test_calc_individual_isoline_vol(self, diagram_r290):
        """calc_individual_isoline accepts 'vol' as isoline/endpoint property."""
        result = diagram_r290.calc_individual_isoline(
            isoline_property="s",
            isoline_value=diagram_r290.convert_from_SI(diagram_r290.state.smass(), "s"),
            starting_point_property="vol",
            starting_point_value=diagram_r290.convert_from_SI(diagram_r290.v_min * 2, "vol"),
            ending_point_property="vol",
            ending_point_value=diagram_r290.convert_from_SI(diagram_r290.v_max / 2, "vol"),
        )
        assert "vol" in result
        assert "p" in result
        assert "T" in result

    def test_calc_individual_isoline_v_deprecated(self, diagram_r290):
        """calc_individual_isoline still accepts 'v' with a FutureWarning."""
        with pytest.warns(FutureWarning, match="'v'.*deprecated"):
            result = diagram_r290.calc_individual_isoline(
                isoline_property="s",
                isoline_value=diagram_r290.convert_from_SI(diagram_r290.state.smass(), "s"),
                starting_point_property="v",
                starting_point_value=diagram_r290.convert_from_SI(diagram_r290.v_min * 2, "vol"),
                ending_point_property="v",
                ending_point_value=diagram_r290.convert_from_SI(diagram_r290.v_max / 2, "vol"),
            )
        assert "vol" in result


# ---------------------------------------------------------------------------
# JSON roundtrip
# ---------------------------------------------------------------------------

class TestJsonRoundtrip:
    def test_to_json_units_serialized(self, tmp_path):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(T="°C", p="bar")
        path = tmp_path / "test.json"
        d.to_json(str(path))
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        assert data["META"]["units"]["T"] == "°C"
        assert data["META"]["units"]["p"] == "bar"

    def test_from_json_restores_units(self, tmp_path):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(T="°C", p="bar", h="kJ/kg")
        d.calc_isolines()
        path = tmp_path / "water.json"
        d.to_json(str(path))

        d2 = FluidPropertyDiagram.from_json(str(path))
        assert d2.units["T"] == "°C"
        assert d2.units["p"] == "bar"
        assert d2.units["h"] == "kJ/kg"
        assert d2.fluid == "water"

    def test_from_json_converts_correctly_after_restore(self, tmp_path):
        d = FluidPropertyDiagram("water")
        d.set_unit_system(p="bar")
        d.calc_isolines()
        path = tmp_path / "water_bar.json"
        d.to_json(str(path))

        d2 = FluidPropertyDiagram.from_json(str(path))
        assert d2.convert_to_SI(1.0, "p") == pytest.approx(1e5)

    def test_from_json_isolines_preserved(self, tmp_path):
        d = FluidPropertyDiagram("water")
        d.calc_isolines()
        n_pressure = len(d.pressure["isolines"])
        path = tmp_path / "water_iso.json"
        d.to_json(str(path))

        d2 = FluidPropertyDiagram.from_json(str(path))
        assert len(d2.pressure["isolines"]) == n_pressure


# ---------------------------------------------------------------------------
# Legacy JSON import (fluprodia 3.x flat format)
# ---------------------------------------------------------------------------

def _make_legacy_json(path, fluid="water", version="3.3"):
    """Write a minimal v3.x-style JSON file to *path*."""
    data = {
        "pressure": {
            "100000.0": {"h": [1e5, 2e5], "T": [300.0, 400.0], "v": [0.001, 0.01],
                         "s": [100.0, 200.0], "p": [1e5, 1e5], "Q": [-1.0, -1.0]},
            "200000.0": {"h": [1e5, 2e5], "T": [300.0, 420.0], "v": [0.0008, 0.008],
                         "s": [90.0, 190.0], "p": [2e5, 2e5], "Q": [-1.0, -1.0]},
        },
        "volume": {},
        "temperature": {},
        "enthalpy": {},
        "entropy": {},
        "quality": {},
        "META": {
            "fluid": fluid,
            "units": {"p": "Pa", "T": "K", "s": "J/kgK", "h": "J/kg", "v": "m^3/kg", "Q": "-"},
            "CoolProp-version": "6.4.1",
            "fluprodia-version": version,
        },
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)


class TestLegacyJsonImport:
    def test_legacy_raises_user_warning(self, tmp_path):
        path = tmp_path / "legacy.json"
        _make_legacy_json(path)
        with pytest.warns(UserWarning, match="3.3"):
            FluidPropertyDiagram.from_json(str(path))

    def test_legacy_isolines_migrated(self, tmp_path):
        path = tmp_path / "legacy.json"
        _make_legacy_json(path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            d = FluidPropertyDiagram.from_json(str(path))
        assert len(d.pressure["isolines"]) == 2
        assert len(d.pressure["isoline_data"]) == 2

    def test_current_format_no_warning(self, tmp_path):
        """Loading a current-format JSON must not emit a legacy UserWarning."""
        d = FluidPropertyDiagram("water")
        d.calc_isolines()
        path = tmp_path / "current.json"
        d.to_json(str(path))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            FluidPropertyDiagram.from_json(str(path))
        legacy_warnings = [
            x for x in w
            if issubclass(x.category, UserWarning)
            and "migrated" in str(x.message)
        ]
        assert len(legacy_warnings) == 0


# ---------------------------------------------------------------------------
# Supported diagram types
# ---------------------------------------------------------------------------

class TestSupportedDiagrams:
    def test_plogv_uses_vol_property(self):
        d = FluidPropertyDiagram("water")
        assert d.supported_diagrams["plogv"]["x_property"] == "vol"

    def test_all_diagram_types_present(self):
        d = FluidPropertyDiagram("water")
        assert set(d.supported_diagrams) == {"Ts", "hs", "logph", "Th", "plogv"}

    def test_invalid_diagram_type_raises(self, diagram_r290):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1)
        with pytest.raises(ValueError, match="not available"):
            diagram_r290.draw_isolines(
                diagram_type="invalid", fig=fig, ax=ax,
                x_min=0, x_max=1, y_min=0, y_max=1
            )
        plt.close(fig)
