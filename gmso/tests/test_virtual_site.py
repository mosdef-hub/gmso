import numpy as np
import pytest
import unyt as u
from sympy import symbols, sympify
from unyt.testing import assert_allclose_units

from gmso import ForceField, Topology
from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.element import element_by_symbol
from gmso.core.virtual_site import VirtualSite
from gmso.core.virtual_type import (
    VirtualPositionType,
    VirtualPotentialType,
    VirtualType,
)
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn


class TestVirturalSite(BaseTest):
    @pytest.fixture(scope="session")
    def virtual_site(self):
        site = Atom(position=[1, 1, 1])
        return VirtualSite(parent_sites=[site])

    @pytest.fixture(scope="session")
    def virtual_type(self):
        v_pot = VirtualPotentialType()
        v_pos = VirtualPositionType()
        v_type = VirtualType(
            virtual_potential=v_pot, virtual_position=v_pos, member_types=("C",)
        )
        return v_type

    def test_virtual_molecule(self, virtual_site):
        assert virtual_site.molecule == virtual_site.parent_sites[0].molecule

    def test_new_site(self, water_system):
        v_site = VirtualSite(parent_sites=water_system.sites)
        assert len(v_site.parent_sites) == 3
        for site in v_site.parent_sites:
            assert site in water_system.sites

    def test_virtual_position(self, virtual_site):
        # Check position as a function of virtual_position_type

        v_pot = VirtualPotentialType(
            expression="5*a*b",
            independent_variables={"a"},
            parameters={"b": 1 * u.kJ},
        )
        v_pos = VirtualPositionType(
            expression="ri*cos(b)",
            independent_variables=["ri"],
            parameters={"b": np.pi * u.radian},
        )
        assert v_pos
        v_type = VirtualType(virtual_potential=v_pot, virtual_position=v_pos)
        virtual_site.virtual_type = v_type  # assign virtual type
        assert_allclose_units(virtual_site.position(), -1 * ([1, 1, 1] * u.nm))

    def test_tip4p_water(self):
        water_ff = ForceField(get_fn("gmso_xmls/test_molecules/tip4p_ice.xml"))
        assert (
            water_ff.virtual_types["HW~OW~HW"].charge == -1.1794 * u.elementary_charge
        )
        doh = 0.09572  # nm
        dhh = 0.15139  # nm
        dom = 0.01577  # nm
        ahoh = 104.52  # degrees
        top = Topology(name="water")
        h1posx = doh * np.sin(ahoh / 360 * np.pi)
        h1posy = doh * np.cos(ahoh / 360 * np.pi)
        poso1 = np.array([0.00000, 0, 0.00000])
        posh1 = np.array([h1posx, h1posy, 0.00000])
        posh2 = np.array([-1 * h1posx, h1posy, 0.00000])
        o1 = Atom(name="O", element=element_by_symbol("O"), position=poso1)
        h1 = Atom(element=element_by_symbol("H"), position=posh1)
        h2 = Atom(element=element_by_symbol("H"), position=posh2)
        np.testing.assert_allclose(
            np.linalg.norm(o1.position - h1.position), doh, rtol=1e-4
        )
        np.testing.assert_allclose(
            np.linalg.norm(h1.position - h2.position), dhh, rtol=1e-4
        )
        for part in [o1, h1, h2]:
            top.add_site(part)
        top.add_connection(Bond(connection_members=(o1, h1)))
        top.add_connection(Bond(connection_members=(o1, h2)))

        ptop = apply(top, water_ff, ignore_params=["bond"], remove_untyped=False)
        np.testing.assert_allclose(
            np.linalg.norm(o1.position - h1.position).value, doh, rtol=1e-4
        )
        np.testing.assert_allclose(
            np.linalg.norm(h1.position - h2.position).value, dhh, rtol=1e-4
        )
        np.testing.assert_allclose(
            np.asin(dhh / 2 / doh), ahoh / 360 * np.pi, rtol=1e-4
        )
        # manually calculate virtual position
        # use position function to get virtual_position
        np.testing.assert_allclose(
            ptop.virtual_sites[0].position().value,
            (
                o1.position
                + 0.13458 * (h1.position - o1.position)
                + 0.13458 * (h2.position - o1.position)
            ).value,
            rtol=1e-4,
        )
        np.testing.assert_allclose(
            np.linalg.norm(ptop.virtual_sites[0].position() - o1.position).value,
            dom,
            rtol=1e-4,
        )
        assert ptop.virtual_sites[0].charge == -1.1794 * u.elementary_charge

        ptop.virtual_sites[0].charge = 1.0 * u.elementary_charge
        assert ptop.virtual_sites[0].charge == 1 * u.elementary_charge

    def test_tip5p_water(self):
        water_ff = ForceField(get_fn("gmso_xmls/test_molecules/tip5p.xml"))
        assert_allclose_units(
            water_ff.virtual_types["HW~OW~HW"].charge, -0.241 * u.elementary_charge
        )
        doh = 0.09572  # nm
        dhh = 0.1513900654  # nm
        dom = 0.070  # nm
        ahoh = 104.52  # degrees
        # amom = 109.47  # degrees
        top = Topology(name="water")
        h1posx = doh * np.sin(ahoh / 360 * np.pi)
        h1posy = doh * np.cos(ahoh / 360 * np.pi)
        poso1 = np.array([0.00000, 0, 0.00000])
        posh1 = np.array([h1posx, h1posy, 0.00000])
        posh2 = np.array([-1 * h1posx, h1posy, 0.00000])
        o1 = Atom(name="O", element=element_by_symbol("O"), position=poso1)
        h1 = Atom(element=element_by_symbol("H"), position=posh1)
        h2 = Atom(element=element_by_symbol("H"), position=posh2)
        np.testing.assert_allclose(
            np.linalg.norm(o1.position - h1.position), doh, rtol=1e-4
        )
        np.testing.assert_allclose(
            np.linalg.norm(h1.position - h2.position), dhh, rtol=1e-4
        )
        for part in [o1, h1, h2]:
            top.add_site(part)
        top.add_connection(Bond(connection_members=(o1, h1)))
        top.add_connection(Bond(connection_members=(o1, h2)))

        ptop = apply(top, water_ff, ignore_params=["bond"], remove_untyped=False)
        np.testing.assert_allclose(
            np.linalg.norm(o1.position - h1.position).value, doh, rtol=1e-4
        )
        np.testing.assert_allclose(
            np.linalg.norm(h1.position - h2.position).value, dhh, rtol=1e-4
        )
        np.testing.assert_allclose(
            np.asin(dhh / 2 / doh), ahoh / 360 * np.pi, rtol=1e-4
        )
        # manually calculate virtual position
        # use position function to get virtual_position

        assert ptop.virtual_sites[0].position().shape == (3,)
        np.testing.assert_allclose(
            ptop.virtual_sites[0].position().value,
            [0.0, -0.0404151276, -0.0571543301],
            atol=1e-4,
        )
        np.testing.assert_allclose(
            ptop.virtual_sites[0].position().value,
            (
                o1.position.value
                - 0.344908268 * (h1.position - o1.position).value
                - 0.344908268 * (h2.position - o1.position).value
                - 6.4437903
                * np.cross(
                    (h1.position - o1.position).value, (h2.position - o1.position).value
                )
            ),
            atol=1e-4,
        )
        np.testing.assert_allclose(
            np.linalg.norm(ptop.virtual_sites[0].position() - o1.position).value,
            dom,
            rtol=1e-4,
        )
        assert_allclose_units(
            ptop.virtual_sites[0].charge, -0.241 * u.elementary_charge
        )
        # assert len(ptop.virtual_sites) == 2 # TODO: This fails because the virtual sites have the same member_classes, so they get the same key

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
