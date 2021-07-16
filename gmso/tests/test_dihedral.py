import pytest
import unyt as u
from pydantic import ValidationError

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.dihedral import Dihedral, LayeredDihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.topology import Topology
from gmso.tests.base_test import BaseTest


class TestDihedral(BaseTest):
    def test_dihedral_nonparametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom1")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        connect = Dihedral(connection_members=[atom1, atom2, atom3, atom4])

        assert connect.connection_type is None

    def test_dihedral_parametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        dihedral_type = DihedralType()

        connect = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
            dihedral_type=dihedral_type,
            name="dihedral_name",
        )

        assert len(connect.connection_members) == 4
        assert connect.connection_type is not None
        assert connect.dihedral_type is not None
        assert connect.name == "dihedral_name"

    def test_dihedral_fake(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(ValidationError):
            Dihedral(connection_members=["fakeatom1", "fakeatom2", 4.2])

    def test_dihedral_fake_dihedraltype(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(ValidationError):
            Dihedral(
                connection_members=[atom1, atom2, atom3, atom4],
                dihedral_type="Fake dihedraltype",
            )

    def test_dihedral_constituent_types(self):
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
        dihtype = DihedralType(
            member_types=[
                atom1.atom_type.name,
                atom2.atom_type.name,
                atom3.atom_type.name,
                atom4.atom_type.name,
            ]
        )
        dih = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )
        dih.dihedral_type = dihtype
        assert "A" in dih.connection_type.member_types
        assert "B" in dih.connection_type.member_types
        assert "C" in dih.connection_type.member_types
        assert "D" in dih.connection_type.member_types

    def test_dihedral_eq(self):
        atom1 = Atom(name="atom1", position=[0, 0, 0])
        atom2 = Atom(name="atom2", position=[1, 0, 0])
        atom3 = Atom(name="atom3", position=[1, 1, 0])
        atom4 = Atom(name="atom4", position=[1, 1, 1])

        ref_dihedral = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        same_dihedral = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        diff_dihedral = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        assert ref_dihedral != same_dihedral
        assert ref_dihedral != diff_dihedral

    def test_add_equivalent_connections(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        dihedral = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
        dihedral_eq = Dihedral(connection_members=[atom4, atom3, atom2, atom1])
        dihedral_not_eq = Dihedral(
            connection_members=[atom4, atom2, atom3, atom1]
        )

        top = Topology()
        top.add_connection(dihedral)
        top.add_connection(dihedral_eq)
        assert top.n_dihedrals == 1

        top.add_connection(dihedral_not_eq)
        assert top.n_dihedrals == 2

    def test_equivalent_members_set(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        dihedral = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
        dihedral_eq = Dihedral(connection_members=[atom4, atom3, atom2, atom1])
        dihedral_not_eq = Dihedral(
            connection_members=[atom4, atom2, atom3, atom1]
        )

        assert (
            tuple(dihedral_eq.connection_members)
            in dihedral.equivalent_members()
        )
        assert (
            tuple(dihedral.connection_members)
            in dihedral_eq.equivalent_members()
        )
        assert not (
            tuple(dihedral.connection_members)
            in dihedral_not_eq.equivalent_members()
        )

    def test_layered_dihedrals(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        dihedral_type1 = DihedralType(
            name=f"layer1",
            expression="kn * (1 + cos(n * a - a0))",
            independent_variables="a",
            parameters={
                "kn": 1.0 * u.K * u.kb,
                "n": 1 * u.dimensionless,
                "a0": 30.0 * u.degree,

            }
        )
        dihedral_type2 = DihedralType(
            name=f"layer2",
            expression="kn * (1 + cos(n * a - a0))",
            independent_variables="a",
            parameters={
                "kn": 1.0 * u.K * u.kb,
                "n": 2 * u.dimensionless,
                "a0": 30.0 * u.degree,

            }
        )
        dihedral_type3 = DihedralType(
            name=f"layer3",
            expression="kn * (1 + cos(n * a - a0))",
            independent_variables="a",
            parameters={
                "kn": 1.0 * u.K * u.kb,
                "n": 3 * u.dimensionless,
                "a0": 30.0 * u.degree,

            }
        )

        connect = LayeredDihedral(
            connection_members=[atom1, atom2, atom3, atom4],
            dihedral_types=[dihedral_type1, dihedral_type2, dihedral_type3],
            name="dihedral_name",
        )

        assert dihedral_type1 in connect.dihedral_types
        assert dihedral_type2 in connect.dihedral_types
        assert dihedral_type3 in connect.dihedral_types

        assert connect.dihedral_types[0].parameters["n"] == 1
        assert connect.dihedral_types[1].parameters["n"] == 2
        assert connect.dihedral_types[2].parameters["n"] == 3
