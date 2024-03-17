import pytest
from pydantic import ValidationError

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.topology import Topology
from gmso.tests.base_test import BaseTest


class TestImproper(BaseTest):
    def test_improper_nonparametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        connect = Improper(connection_members=[atom1, atom2, atom3, atom4])

        assert connect.connection_type is None

    def test_improper_parametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        improper_type = ImproperType()

        connect = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
            improper_type=improper_type,
            name="improper_name",
        )

        assert len(connect.connection_members) == 4
        assert connect.connection_type is not None
        assert connect.improper_type is not None
        assert connect.name == "improper_name"

    def test_improper_fake(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(TypeError):
            Improper(connection_members=["fakeatom1", "fakeatom2", 4.2])

    def test_improper_fake_impropertype(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(ValidationError):
            Improper(
                connection_members=[atom1, atom2, atom3, atom4],
                connection_type="Fake impropertype",
            )

    def test_improper_constituent_types(self):
        atom1 = Atom(
            name="atom1", position=[0, 0, 0], atom_type=AtomType(name="A")
        )
        atom2 = Atom(
            name="atom2", position=[1, 0, 0], atom_type=AtomType(name="B")
        )
        atom3 = Atom(
            name="atom3", position=[1, 1, 0], atom_type=AtomType(name="C")
        )
        atom4 = Atom(
            name="atom4", position=[1, 1, 4], atom_type=AtomType(name="D")
        )
        imptype = ImproperType(
            member_types=[
                atom1.atom_type.name,
                atom2.atom_type.name,
                atom3.atom_type.name,
                atom4.atom_type.name,
            ]
        )
        imp = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )
        imp.improper_type = imp.connection_type = imptype
        assert "A" in imp.connection_type.member_types
        assert "B" in imp.connection_type.member_types
        assert "C" in imp.connection_type.member_types
        assert "D" in imp.connection_type.member_types

    def test_improper_eq(self):
        atom1 = Atom(name="atom1", position=[0, 0, 0])
        atom2 = Atom(name="atom2", position=[1, 0, 0])
        atom3 = Atom(name="atom3", position=[1, 1, 0])
        atom4 = Atom(name="atom4", position=[1, 1, 1])

        ref_improper = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        same_improper = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        diff_improper = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        assert ref_improper != same_improper
        assert ref_improper != diff_improper

    def test_add_equivalent_connections(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        improper = Improper(connection_members=[atom1, atom2, atom3, atom4])
        improper_eq = Improper(connection_members=[atom1, atom3, atom2, atom4])
        improper_not_eq = Improper(
            connection_members=[atom2, atom3, atom1, atom4]
        )

        top = Topology()
        top.add_connection(improper)
        top.add_connection(improper_eq)
        assert top.n_impropers == 1
        top.add_connection(improper_not_eq)
        assert top.n_impropers == 2

    def test_equivalent_members_set(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        improper = Improper(connection_members=[atom1, atom2, atom3, atom4])
        improper_eq = Improper(connection_members=[atom1, atom3, atom2, atom4])
        improper_not_eq = Improper(
            connection_members=[atom2, atom3, atom1, atom4]
        )

        assert (
            tuple(improper_eq.connection_members)
            in improper.equivalent_members()
        )
        assert (
            tuple(improper.connection_members)
            in improper_eq.equivalent_members()
        )
        assert not (
            tuple(improper.connection_members)
            in improper_not_eq.equivalent_members()
        )

    def test_sort_improper_types(self):
        from gmso.utils.sorting import sort_by_classes, sort_by_types

        atom1 = Atom(
            name="atom1", position=[0, 0, 0], atom_type=AtomType(name="A")
        )
        atom2 = Atom(
            name="atom2", position=[1, 0, 0], atom_type=AtomType(name="B")
        )
        atom3 = Atom(
            name="atom3", position=[1, 1, 0], atom_type=AtomType(name="C")
        )
        atom4 = Atom(
            name="atom4", position=[1, 1, 4], atom_type=AtomType(name="D")
        )

        consituentList = [
            atom2.atom_type.name,
            atom4.atom_type.name,
            atom3.atom_type.name,
            atom1.atom_type.name,
        ]
        imptype = ImproperType(
            member_types=consituentList, member_classes=consituentList
        )

        expected_sortingList = tuple(
            [atom.atom_type.name for atom in [atom2, atom3, atom4, atom1]]
        )
        assert sort_by_classes(imptype) == expected_sortingList
        assert sort_by_types(imptype) == expected_sortingList

    def test_sorting_improper_based_on_impropertype(self):
        from gmso.exceptions import MissingParameterError
        from gmso.utils.sorting import sort_by_classes, sort_by_types

        def sort_improper_connection_members(improper):
            if improper.improper_type is None:
                return improper.connection_members
            improper_classes = improper.improper_type.member_classes
            improper_types = improper.improper_type.member_types
            if improper_classes:
                order_improperList = improper_classes
                orderStr = "atomclass"  # String to access site attribute
            elif improper_types:
                order_improperList = improper_types
                orderStr = "name"  # String to access site attribute
            else:
                missing_types = [site.atom_type.atomclass for site in improper]
                raise MissingParameterError(
                    improper.improper_type, missing_types
                )

            # get the site atomtypes and make a dictionary map to match to the order_improperList
            cmemList = improper.connection_members
            assert order_improperList[0] == getattr(
                cmemList[0].atom_type, orderStr
            )  # first atoms should be the same
            first_site = cmemList[0]
            middle_sitesList = []
            for site in cmemList[1:]:
                if getattr(site.atom_type, orderStr) == order_improperList[-1]:
                    last_site = site
                else:
                    middle_sitesList.append(site)
            assert (
                len(middle_sitesList) == 2
            ), f"The improper_type {improper.improper_type} could not find 2 middle sites from {middle_sitesList}"
            middle_sitesList = sorted(
                middle_sitesList,
                key=lambda site: getattr(site.atom_type, orderStr),
            )
            return [first_site] + middle_sitesList + [last_site]

        atom1 = Atom(
            name="atom1",
            position=[0, 0, 0],
            atom_type=AtomType(name="A", atomclass="A"),
        )
        atom2 = Atom(
            name="atom2",
            position=[1, 0, 0],
            atom_type=AtomType(name="B", atomclass="B"),
        )
        atom3 = Atom(
            name="atom3",
            position=[1, 1, 0],
            atom_type=AtomType(name="C", atomclass="C"),
        )
        atom4 = Atom(
            name="atom4",
            position=[1, 1, 4],
            atom_type=AtomType(name="D", atomclass="D"),
        )

        connect = Improper(connection_members=[atom2, atom1, atom4, atom3])

        consituentList = [
            atom2.atom_type.name,
            atom4.atom_type.name,
            atom3.atom_type.name,
            atom1.atom_type.name,
        ]

        imptype = ImproperType(
            member_types=consituentList, member_classes=consituentList
        )
        connect.improper_type = imptype
        expected_membersList = [atom2, atom3, atom4, atom1]
        assert sort_by_types(connect.improper_type) == tuple(
            [site.atom_type.name for site in expected_membersList]
        )
        assert (
            tuple([site.atom_type for site in connect.connection_members])
            != imptype.member_types
        )
        assert sort_improper_connection_members(connect) == expected_membersList

    def test_applied_improper_updates_connection_members(self, benzeneTopology):
        improper = benzeneTopology.impropers[0]
        classes_connectionList = tuple(
            [site.atom_type.atomclass for site in improper.connection_members]
        )
        classes_typeList = improper.improper_type.member_classes
        assert classes_connectionList == classes_typeList
