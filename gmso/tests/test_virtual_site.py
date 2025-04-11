import unyt as u

from gmso.core.virtual_site import VirtualSite
from gmso.core.virtual_type import (
    VirtualPositionType,
    VirtualPotentialType,
    VirtualType,
)
from gmso.tests.base_test import BaseTest


class TestVirturalSite(BaseTest):
    def test_new_site(self, water_system):
        v_site = VirtualSite(parent_atoms=water_system.sites)
        assert len(v_site.parent_atoms) == 3
        for site in v_site.parent_atoms:
            assert site in water_system.sites

    def test_virtual_position(self):
        # TODO: Check position as a function of virtual_position_type
        pass

    def test_virtual_type(self):
        v_pot = VirtualPotentialType(
            charge=0.1,
            expression="5*a*b",
            independent_variables={"a"},
            parameters={"b": 1 * u.kJ},
        )
        assert v_pot
        v_pos = VirtualPositionType(
            expression="ri+b",
            independent_variables="ri",
            parameters={"b": [1, 0, 0] * u.nm},
        )
        assert v_pos
        v_type = VirtualType(
            virtual_potential=v_pot, virtual_position=v_pos, member_types=("C", "C")
        )
        assert v_type

        v_pos2 = VirtualPositionType(
            expression="ri+b",
            independent_variables="ri",
            parameters={"b": [1, 0, 0] * u.nm},
        )
        assert v_pos2 == v_pos
