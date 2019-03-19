import numpy as np
import pytest 
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.box import Box
from topology.testing.utils import allclose
from topology.core.atom_type import AtomType
from topology.core.element import Element
from topology.external.convert_openmm import to_modeller


def test_to_modeller():
    top = Topology()
    top.box = Box(lengths=[1,1,1])
    H = Element(name='H', symbol='H', mass=1)
    site1 = Site(name='site1',
            element=H,
            atom_type=AtomType(name="at1",
                               mass=H.mass)
            )
    top.add_site(site1)
    to_modeller(top)

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
    modeller = to_modeller(top)
    n_modeller_atoms = len([i for i in modeller.topology.atoms()])

    assert n_topology_sites == n_modeller_atoms
