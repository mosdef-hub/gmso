import unyt as u

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.topology import Topology
from gmso.core.views import PotentialFilters
from gmso.tests.base_test import BaseTest
from gmso.utils.misc import unyt_to_hashable


class TestViews(BaseTest):
    def test_view_atom_types_typed_ar_system(self, n_typed_ar_system):
        atom_types = n_typed_ar_system().atom_types()
        assert len(atom_types) == 1
        atom_types_unique = n_typed_ar_system().atom_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(atom_types_unique) == 1

    def test_ethane_views(self, typed_ethane):
        atom_types = typed_ethane.atom_types
        unique_atomtypes = atom_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(atom_types) == len(unique_atomtypes)

        bond_types = typed_ethane.bond_types
        unique_bondtypes = typed_ethane.bond_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(bond_types) == len(unique_bondtypes)
        assert typed_ethane._potentials_count["bond_types"] == len(bond_types)

        angle_types = typed_ethane.angle_types
        unique_angletypes = typed_ethane.angle_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(angle_types) == len(unique_angletypes)
        assert typed_ethane._potentials_count["angle_types"] == len(bond_types)

        dihedral_types = typed_ethane.dihedral_types
        unique_dihedraltypes = typed_ethane.dihedral_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(unique_dihedraltypes) == len(dihedral_types)
        assert typed_ethane._potentials_count["dihedral_types"] == len(
            dihedral_types
        )

    def test_custom_filter(self):
        def str_expr_dict_param(potentials):
            visited = set()
            for potential in potentials:
                print(potential.parameters.values())
                pot_id = (
                    str(potential.expression),
                    frozenset(potential.parameters),
                    tuple(
                        unyt_to_hashable(list(potential.parameters.values()))
                    ),
                )
                if pot_id not in visited:
                    visited.add(pot_id)
                    yield potential

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
                    atom_type=atom_type1.clone()
                    if j % 2 == 0
                    else atom_type2.clone()
                )
            )

        custom_top.update_topology()
        assert len(custom_top.atom_types) == 200
        assert len(custom_top.atom_types(filter_by=str_expr_dict_param)) == 2
