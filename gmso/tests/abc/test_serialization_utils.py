import pytest
import unyt as u

from gmso.abc.serialization_utils import dict_to_unyt, unyt_to_dict
from gmso.tests.base_test import BaseTest


class TestSerializationUtils(BaseTest):
    @pytest.fixture
    def custom_json_serializable_class(self):
        pass

    def test_unyt_to_dict_quantity(self):
        quantity = 5.0 * u.nm
        quantity_dict = unyt_to_dict(quantity)
        assert quantity_dict["array"] == 5.0
        assert quantity_dict["unit"] == "nm"

    def test_unyt_to_dict_array(self):
        array = [2.0, 3.5, 5.6] * u.K
        array_dict = unyt_to_dict(array)
        assert array_dict["array"] == [2.0, 3.5, 5.6]
        assert array_dict["unit"] == "K"

    def test_unyt_to_dict_non_unyt(self):
        with pytest.raises(TypeError):
            unyt_to_dict(5.0)

    def test_dict_to_unyt_shallow(self):
        unyt_dict = {"my_quantity": {"array": [7.0, 2.5, 6.3], "unit": "K"}}
        dict_to_unyt(unyt_dict)
        assert u.allclose_units(
            unyt_dict["my_quantity"], u.unyt_array([7.0, 2.5, 6.3], u.K)
        )

    def test_dict_to_unyt_nested(self):
        unyt_dict = {
            "level1": {
                "my_quantity": {"array": [7.0, 2.5, 6.3], "unit": "K"},
                "level2": {
                    "my_quantity": {"array": [-10, 2.5, 100], "unit": "C"},
                },
            }
        }
        dict_to_unyt(unyt_dict)

        assert u.allclose_units(
            unyt_dict["level1"]["my_quantity"],
            u.unyt_array([7.0, 2.5, 6.3], u.K),
        )
        assert u.allclose_units(
            unyt_dict["level1"]["level2"]["my_quantity"],
            u.unyt_array([-10, 2.5, 100], u.C),
        )
