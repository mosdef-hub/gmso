from gmso.core.views import PotentialFilters
from gmso.tests.base_test import BaseTest


class TestViews(BaseTest):
    def test_view_atom_types_typed_ar_system(self, n_typed_ar_system):
        atom_types = n_typed_ar_system().atom_types(repeat=True)
        assert len(atom_types) == 100
        atom_types_unique = n_typed_ar_system().atom_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        )
        assert len(atom_types_unique) == 1

    def test_ethane_views(self, typed_ethane):
        atom_types = typed_ethane.atom_types
        unique_atomtypes = atom_types(repeat=True)
        assert len(atom_types) == len(unique_atomtypes)
