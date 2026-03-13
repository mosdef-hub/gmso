import mbuild as mb
import numpy as np
import pytest
import unyt as u

from gmso.core.atom import Atom
from gmso.core.dihedral import Dihedral
from gmso.utils.geometry import moment_of_inertia
from gmso.utils.io import run_from_ipython
from gmso.utils.misc import unyt_to_hashable
from gmso.utils.slicing import (
    slice_topology_by_molecule,
    slice_topology_by_residue,
)
from gmso.utils.sorting import sort_connection_members, sort_connection_strings


def test_unyt_to_hashable():
    hash(unyt_to_hashable(None))
    hash(unyt_to_hashable(1 * u.nm))
    hash(unyt_to_hashable([4, 4] * u.nm))

    assert hash(unyt_to_hashable(1 * u.nm)) == hash(unyt_to_hashable(10 * u.angstrom))
    assert hash(unyt_to_hashable(1 * u.kg)) == hash(unyt_to_hashable(1000 * u.g))

    assert hash(unyt_to_hashable(1 * u.nm)) != hash(unyt_to_hashable(1.01 * u.nm))
    assert hash(unyt_to_hashable(1 * u.nm)) != hash(unyt_to_hashable(1.01 * u.second))
    assert hash(unyt_to_hashable(1 * u.nm)) != hash(unyt_to_hashable([1, 1] * u.nm))


def test_has_ipython():
    assert run_from_ipython() is False


def test_sorting():
    with pytest.raises(TypeError):
        sort_connection_members([1])
    atom1 = Atom(name="atom1")
    atom2 = Atom(name="atom2")
    atom3 = Atom(name="atom3")
    atom4 = Atom(name="atom4")

    connect = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
    with pytest.raises(ValueError):
        sort_connection_members(connect, sort_by="error")

    bondList = [atom3.name, atom2.name]
    angleList = [atom3.name, atom2.name, atom1.name]
    dihList = [atom3.name, atom2.name, atom1.name, atom4.name]

    assert sort_connection_strings(bondList) == ("atom2", "atom3")
    assert sort_connection_strings(angleList) == ("atom1", "atom2", "atom3")
    assert sort_connection_strings(dihList) == (
        "atom4",
        "atom1",
        "atom2",
        "atom3",
    )
    assert sort_connection_strings(dihList, improperBool=True) == (
        "atom3",
        "atom1",
        "atom2",
        "atom4",
    )


def test_moment_of_inertia():
    xyz = [
        np.array([-0.5, -0.5, 0]),
        np.array([0.5, -0.5, 0]),
        np.array([0.5, 0.5, 0]),
        np.array([-0.5, 0.5, 0]),
    ]
    tensor = moment_of_inertia(xyz=xyz, masses=np.array([1.0 for i in xyz]))
    assert np.array_equal(tensor, np.array([1, 1, 2]))

    xyz = [
        np.array([0, 0, 0]),
        np.array([1, 0, 0]),
        np.array([1, 1, 0]),
        np.array([0, 1, 0]),
    ]
    tensor = moment_of_inertia(
        xyz=xyz,
        center=np.array((0.5, 0.5, 0)),
        masses=np.array([1.0 for i in xyz]),
    )
    assert np.array_equal(tensor, np.array([1, 1, 2]))


def test_slice_by_molecule():
    benzene = mb.load("c1ccccc1", smiles=True)
    benzene.name = "Benzene"
    ethane = mb.load("CC", smiles=True)
    ethane.name = "Ethane"

    system = mb.fill_box(compound=[benzene, ethane], n_compounds=[2, 2], box=[2, 2, 2])
    topology = system.to_gmso()
    topology.identify_connections()

    single_benzene_top = slice_topology_by_molecule(topology, "Benzene", 0)
    assert single_benzene_top.n_sites == 12

    all_benzene_top = slice_topology_by_molecule(topology, "Benzene")
    assert all_benzene_top.n_sites == 24


def test_slice_by_residue():
    benzene = mb.load("c1ccccc1", smiles=True)
    benzene.name = "Benzene"
    ethane = mb.load("CC", smiles=True)
    ethane.name = "Ethane"

    system = mb.fill_box(compound=[benzene, ethane], n_compounds=[2, 2], box=[2, 2, 2])
    topology = system.to_gmso()
    topology.identify_connections()

    single_ethane_top = slice_topology_by_residue(topology, "Ethane", 0)
    assert single_ethane_top.n_sites == 8

    all_ethane_top = slice_topology_by_residue(topology, "Ethane")
    assert all_ethane_top.n_sites == 16
