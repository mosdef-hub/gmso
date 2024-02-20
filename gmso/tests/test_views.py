import pytest
import unyt as u

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.topology import Topology
from gmso.core.views import PotentialFilters
from gmso.tests.base_test import BaseTest
from gmso.utils.misc import unyt_to_hashable


def str_expr_dict_param(potential):
    pot_id = (
        str(potential.expression),
        frozenset(potential.parameters),
        tuple(unyt_to_hashable(list(potential.parameters.values()))),
    )
    return pot_id


class TestViews(BaseTest):
    @pytest.fixture(scope="session")
    def custom_top(self):
        custom_top = Topology()
        atom_type1 = AtomType(
            name="dummy_1",
            expression="x+y*z",
            parameters={"x": 2.0 * u.nm, "z": 33000 * u.dimensionless},
            independent_variables={"y"},
        )

        atom_type2 = AtomType(
            name="dummy_2",
            expression="x+y*z",
            parameters={"x": 5.0 * u.nm, "z": 33000 * u.dimensionless},
            independent_variables={"y"},
        )

        for j in range(200):
            custom_top.add_site(
                Atom(
                    atom_type=(
                        atom_type1.clone() if j % 2 == 0 else atom_type2.clone()
                    )
                )
            )

        for j in range(0, 200, 2):
            bond = Bond(
                connection_members=[
                    custom_top.sites[j],
                    custom_top.sites[j + 1],
                ],
                bond_type=BondType(member_classes=(str(j), str(j + 1))),
            )
            custom_top.add_connection(bond)

        custom_top.update_topology()
        return custom_top

    def test_view_atom_types_typed_ar_system(self, n_typed_ar_system):
        atom_types = n_typed_ar_system().atom_types()
        assert len(atom_types) == 1
        atom_types_unique = n_typed_ar_system().atom_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(atom_types_unique) == 1

    def test_ethane_views(self, typed_ethane):
        # test filters
        atom_types = typed_ethane.atom_types
        unique_atomtypes = atom_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(atom_types) == 8
        assert len(unique_atomtypes) == 2

        bond_types = typed_ethane.bond_types
        unique_bondtypes = typed_ethane.bond_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(bond_types) == 7
        assert len(unique_bondtypes) == 2
        assert typed_ethane._potentials_count["bond_types"] == len(bond_types)

        angle_types = typed_ethane.angle_types
        unique_angletypes = typed_ethane.angle_types(
            filter_by=PotentialFilters.UNIQUE_SORTED_NAMES
        )
        unique_angletypes_no_symmetries = typed_ethane.angle_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(angle_types) == 12
        assert len(unique_angletypes) == 2
        assert len(unique_angletypes_no_symmetries) == 3
        assert typed_ethane._potentials_count["angle_types"] == len(angle_types)

        dihedral_types = typed_ethane.dihedral_types
        unique_dihedraltypes = typed_ethane.dihedral_types(
            filter_by=PotentialFilters.UNIQUE_SORTED_NAMES
        )
        assert len(dihedral_types) == 9
        assert len(unique_dihedraltypes) == 1
        assert typed_ethane._potentials_count["dihedral_types"] == len(
            dihedral_types
        )

    def test_custom_filter(self, custom_top):
        assert len(custom_top.atom_types) == 200
        assert len(custom_top.atom_types(filter_by=str_expr_dict_param)) == 2

    def test_get_index(self, custom_top):
        assert custom_top.get_index(custom_top.sites[50].atom_type) == 50

    def test_bondtype_views(self, custom_top):
        assert len(custom_top.bond_types) == 100
        assert len(custom_top.bond_types(filter_by=str_expr_dict_param)) == 1

    def test_call(self, custom_top):
        atom_types = custom_top.atom_types
        atom_types_new = atom_types()
        atom_types_new_different_filter = atom_types(
            filter_by=str_expr_dict_param
        )

        assert id(atom_types) == id(atom_types_new)
        assert id(atom_types) != atom_types_new_different_filter

    def test_default_filters(self, custom_top):
        bond_types = custom_top.bond_types(
            filter_by=PotentialFilters.UNIQUE_EXPRESSION
        )
        bond_types_params = bond_types(
            filter_by=PotentialFilters.UNIQUE_PARAMETERS
        )
        bond_types_name_class = bond_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )

        bond_types_repeat = bond_types(
            filter_by=PotentialFilters.REPEAT_DUPLICATES
        )

        assert len(bond_types) == 1
        assert len(bond_types_params) == 1
        assert len(bond_types_name_class) == 100
        assert len(bond_types_repeat) == 100
