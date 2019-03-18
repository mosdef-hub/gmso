import mbuild as mb
import unyt as u
import simtk.unit

from topology.core.topology import Topology
from topology.core.site import Site

from simtk import openmm


def to_openmm(topology):
    # Get into OpenMM form 
    openmm_top = openmm.app.Topology()
    
    # Get topology.positions into OpenMM form
    top_pos = topology.positions()
    value = [i.value for i in topology.positions()]
    openmm_unit = topology.positions().units
    openmm_pos = simtk.unit.Quantity(value=value,
            unit=openmm_unit)

    # Not sure what a 'Chain' is
    chain = openmm_top.addChain()

    # Is a site equivalent to a residue?
    residues = topology.site_list

    for site in topology.site_list:
        openmm_top.addResidue(name=site.name, chain=chain)

    modeller = openmm.app.Modeller(openmm_top, openmm_pos)

    return modeller
