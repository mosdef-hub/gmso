import pytest
import unyt as u
from sympy import symbols, sympify
from unyt.testing import assert_allclose_units
import numpy as np

from gmso.core.atom import Atom
from gmso.core.virtual_site import VirtualSite
from gmso.core.virtual_type import (
    VirtualPositionType,
    VirtualPotentialType,
    VirtualType,
)
from gmso.tests.base_test import BaseTest


class TestVirturalSite(BaseTest):
    @pytest.fixture(scope="session")
    def virtual_site(self):
        site = Atom(position=[1,1,1])
        return VirtualSite(parent_sites=[site])

    @pytest.fixture(scope="session")
    def virtual_type(self):
        v_pot = VirtualPotentialType()
        v_pos = VirtualPositionType()
        v_type = VirtualType(
            virtual_potential=v_pot, virtual_position=v_pos, member_types=("C",)
        )
        return v_type

    def test_new_site(self, water_system):
        v_site = VirtualSite(parent_sites=water_system.sites)
        assert len(v_site.parent_sites) == 3
        for site in v_site.parent_sites:
            assert site in water_system.sites

    def test_virtual_position(self, virtual_site):
        # Check position as a function of virtual_position_type
        # TODO: check for arrays for b, sin norm
        # TODO: checkk all gromacs potential forms

        v_pot = VirtualPotentialType(
            expression="5*a*b",
            independent_variables={"a"},
            parameters={"b": 1 * u.kJ},
        )
        v_pos = VirtualPositionType(
            expression="ri*cos(b)",
            independent_variables=["ri"],
            parameters={"b": np.pi*u.radian},
        )
        assert v_pos
        v_type = VirtualType(
            virtual_potential=v_pot, virtual_position=v_pos
        )
        virtual_site.virtual_type = v_type # assign virtual type
        assert_allclose_units(virtual_site.position(), -1*([1,1,1]*u.nm))


    def test_tip4p_water(self):
        pass

    def test_tip5p_water(self):
        pass


    def test_virtual_type(self):
        v_pot = VirtualPotentialType(
            expression="5*a*b",
            independent_variables={"a"},
            parameters={"b": 1 * u.kJ},
        )
        assert v_pot
        v_pos = VirtualPositionType(
            expression="ri+rj+b",
            independent_variables=["ri", "rj"],
            parameters={"b": [1, 0, 0] * u.nm},
        )
        assert v_pos
        v_type = VirtualType(
            virtual_potential=v_pot, virtual_position=v_pos, member_types=("C", "C")
        )
        assert v_type

        v_pos2 = VirtualPositionType(
            expression="ri+rj+b",
            independent_variables=["ri", "rj"],
            parameters={"b": [1, 0, 0] * u.nm},
        )
        assert v_pos2 == v_pos

    def test_expression(self, virtual_type):
        assert virtual_type.virtual_position.expression == sympify(
            "ri + b*(rj-ri+a*(rk-rj))"
        )
        assert virtual_type.virtual_potential.expression == sympify(
            "4*epsilon*((sigma/r)**12 - (sigma/r)**6)"
        )

    def test_members(self, virtual_type):
        assert virtual_type.member_types == ("C",)
        assert virtual_type.member_classes is None

    def test_charge(self, virtual_type):
        assert virtual_type.charge == 0.0 * u.elementary_charge

    def test_doi(self, virtual_type):
        assert virtual_type.doi == ""
        v_type = virtual_type.clone()
        v_type.doi = "10.2.example"
        assert v_type.doi == "10.2.example"

    def test_equality(self, virtual_type):
        v_pot = VirtualPotentialType()
        v_pos = VirtualPositionType()
        v_type = VirtualType(
            virtual_potential=v_pot, virtual_position=v_pos, member_types=("C",)
        )
        assert v_type == virtual_type

    def test_clone(self, virtual_type):
        v_type = virtual_type.clone()
        assert v_type == virtual_type

    def test_setters(self):
        new_type = VirtualType()
        new_type.name = "SettingName"
        new_type.charge = -1.0 * u.Coulomb

        new_type.virtual_potential.independent_variables = "r"
        new_type.virtual_potential.parameters = {
            "sigma": 1 * u.nm,
            "epsilon": 10 * u.Unit("kcal / mol"),
        }
        new_type.virtual_potential.expression = "r * sigma * epsilon"
        assert new_type.name == "SettingName"
        assert_allclose_units(new_type.charge, -1.0 * u.Coulomb, rtol=1e-5, atol=1e-8)
        assert new_type.virtual_potential.independent_variables == {symbols("r")}
        assert new_type.virtual_potential.parameters == {
            "sigma": 1 * u.nm,
            "epsilon": 10 * u.Unit("kcal / mol"),
        }
        assert new_type.virtual_potential.expression == sympify("r * sigma * epsilon")
