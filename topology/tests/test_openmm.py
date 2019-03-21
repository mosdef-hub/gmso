import numpy as np
import pytest 
import unyt as u
import simtk.unit

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.box import Box
from topology.testing.utils import allclose
from topology.core.atom_type import AtomType
from topology.core.element import Element
from topology.external.convert_openmm import to_openmm


def test_openmm_modeller():
    top = Topology()
    top.box = Box(lengths=[1,1,1])
    H = Element(name='H', symbol='H', mass=1)
    site1 = Site(name='site1',
            element=H,
            atom_type=AtomType(name="at1",
                               mass=H.mass)
            )
    top.add_site(site1)
    to_openmm(top, openmm_object='modeller')


def test_openmm_topology():
    top = Topology()
    top.box = Box(lengths=[1,1,1])
    H = Element(name='H', symbol='H', mass=1)
    site1 = Site(name='site1',
            element=H,
            atom_type=AtomType(name="at1",
                               mass=H.mass)
            )
    top.add_site(site1)
    to_openmm(top, openmm_object='topology')

def test_n_atoms():
    top = Topology()
    top.box = Box(lengths=[1,1,1])
    H = Element(name='H', symbol='H', mass=1)
    site1 = Site(name='site1',
            element=H,
            atom_type=AtomType(name="at1",
                               mass=H.mass)
            )
    for i in range(10):
        top.add_site(site1)

    n_topology_sites = len(top.site_list)
    modeller = to_openmm(top, openmm_object='modeller')
    n_modeller_atoms = len([i for i in modeller.topology.atoms()])

    assert n_topology_sites == n_modeller_atoms

def test_box_dims():
    top = Topology()
    top.box = Box(lengths=[1,1,1])
    H = Element(name='H', symbol='H', mass=1)
    site1 = Site(name='site1',
            element=H,
            atom_type=AtomType(name="at1",
                               mass=H.mass)
            )
    for i in range(10):
        top.add_site(site1)

    n_topology_sites = len(top.site_list)
    omm_top = to_openmm(top)
    topology_lengths = top.box.lengths
    omm_lengths = omm_top.getUnitCellDimensions()

    assert (topology_lengths.value == omm_lengths._value).all()

def test_particle_positions():
    top = Topology()
    top.box = Box(lengths=[2,2,2])
    H = Element(name='H', symbol='H', mass=1)
    site1 = Site(name='site1',
            element=H,
            atom_type=AtomType(name="at1",
                               mass=H.mass)
            )
    site1.position=(1,1,1)*u.nanometer
    top.add_site(site1)
    omm_top = to_openmm(top, openmm_object='modeller')
    assert (omm_top.positions._value == top.positions().value).all()

def test_position_units():
    top = Topology()
    top.box = Box(lengths=[1,1,1])
    H = Element(name='H', symbol='H', mass=1)
    site1 = Site(name='site1',
            element=H,
            atom_type=AtomType(name="at1",
                               mass=H.mass)
            )
    for i in range(10):
        top.add_site(site1)

    n_topology_sites = len(top.site_list)
    omm_top = to_openmm(top, openmm_object='modeller')

    assert isinstance(omm_top.positions.unit, type(simtk.unit.nanometer))
