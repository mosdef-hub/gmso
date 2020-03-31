from warnings import warn
import unyt as u
import numpy as np

import gmso
from gmso.utils.io import import_, has_openmm, has_simtk_unit


if has_openmm & has_simtk_unit:
    simtk_unit = import_('simtk.unit')
    from simtk.openmm.app import *
    from simtk.openmm import *

def from_openmm(openmm_object, refer_type=True):
    """Convert an openmm Topology into a gmso.Topology

    Convert an OpenMM object, either Topology or Modeller, to a
    gmso.Topology. User has the option to either refer the types
    information or not.
    Mapping plan:
    OpenMM: Topology - Chains - Residues - Atoms
    GMSO: Topology - SubTopology - SubTopology - Sites

    Parameters
    ----------
    openmm_object : 'simtk.openmm.app.Topology'
        OpenMM Topology object that need to be converted.
    refer_type : bool, optional, default=True
        Whether or not to refer types information.
    Return
    ------
    top : gmso.Topology
        Typed or untyped GMSO Topology object.
    """

def from_openmm_topology(openmm_topology):
    """Convert an openmm Topology to a gmso Topology

    Helper function for the main from_openmm method.
    Specifically handle openmm Topology. Mapping:
    GMSO Topology: Top - Chain_Subtop -

    Parameters
    ----------
    openmm_topology : 'simtk.openmm.app.Topology'
        OpenMM Topology object that need to be converted.

    Return
    ------
    top : gmso.Topology
        Typed or untyped GMSO Topology object.
    """
    msg = 'Given object is not an OpenMM Topology'
    assert isinstance(openmm_topology, openmm.Topology), msg

    # Initialize GMSO Topology object
    top = gmso.Topology()
    top.name = 'Topology'

    # Convert box information
    lengths = openmm_topology.getPeriodicBoxVectors()[0] * u.nm
    angles = openmm_topology.getPeriodicBoxVectors()[1] * u.degree
    top.box = gmso.Box(lengths=lengths, angles=angles)

    # Convert chains-residues-atoms information
    site_map = dict() # mapping atom -> site
    for chain in openmm_topology.chains():
        chain_name = chain.name + chain.id
        gmso_chain = gmso.SubTopology(name=chain_name,
                                parent=top)
        for residue in chain.residues():
            residue_name = residue.name + residue.id
            gmso_residue = gmso.SubTopology(name=res_name,
                                        parent=gmso_chain)
            for atom in residue.atoms():
                site = gmso.Site(name=gmso.name)
                gmso_residue.add_site(site)
                site_map[atom] = site
        top.add_subtop(chain)

    # Convert bonds information
    for bond in openmm_topology.bonds():
        top_connection = gmso.Bond(
        connection_members=[site_map[bond.atom1],
                            site_map[bond.atom2]])
        top.add(top_connection)

    return top

def from_openmm_modeller(openmm_modeller):
    """Convert an openmm Modeller to a gmso Topology

    Helper function for the main from_openmm method.
    Specifically handle openmm Modeller.

    Parameters
    ----------
    openmm_modeller : 'simtk.openmm.app.Modeller'
        OpenMM Modeller object that need to be converted.

    Return
    ------
    top : gmso.Topology
        Typed or untyped GMSO Topology object.
    """
    msg = 'Given object is not an OpenMM Modeller'
    assert isinstance(openmm_modeller, openmm.Modeller), msg

    # Initialize GMSO Topology
    top = gmso.Topology()
    top.name = 'Topology'

    # Separate Top and positions information
    mm_top = openmm_modeller.topology
    mm_pos = openmm_modeller.position

    # Convert box information
    lengths = mm_top.getPeriodicBoxVectors()[0] * u.nm
    angles = mm_top.getPeriodicBoxVectors()[1] * u.degree
    top.box = gmso.Box(lengths=lengths, angles=angles)

    # Convert topology information
    site_map = dict() # mapping atom -> site
    for chain in mm_top.chains():
        chain_name = chain.name + chain.id
        gmso_chain = gmso.SubTopology(name=chain_name,
                                parent=top)
        for residue in chain.residues():
            residue_name = residue.name + residue.id
            gmso_residue = gmso.SubTopology(name=res_name,
                                        parent=gmso_chain)
            for atom in residue.atoms():
                # Assume things are in nm now, will need
                # to actually read from the openmm position
                # itself to determine the unit
                pos = mm_pos.pop(0)._value * u.nm
                site = gmso.Site(name=gmso.name)
                gmso_residue.add_site(site)
                site_map[atom] = site
        top.add_subtop(chain)

    # Convert bonds information
    for bond in openmm_topology.bonds():
        top_connection = gmso.Bond(
        connection_members=[site_map[bond.atom1],
                            site_map[bond.atom2]])
        top.add(top_connection)

    return top

def from_openmm_system(openmm_system, refer_type=True):
    """Convert an openmm System to a gmso Topology

    Helper function for the main from_openmm method.
    Specifically handle openmm System.

    Parameters
    ----------
    openmm_system : 'simtk.openmm.app.System'
        OpenMM Topology object that need to be converted.
    refer_type : bool, optional, default=True
        Whether to refer types information

    Return
    ------
    top : gmso.Topology
        Typed or untyped GMSO Topology object.
    """


def to_openmm(topology, openmm_object='topology'):
    """
    Convert an untyped topology object to an untyped OpenMM modeller or topology. 
    This is useful if it's preferred to atom-type a system within OpenMM. 
    See http://openmm.org for more information.

    Parameters
    ----------
    topology : `Topology` object
        An untyped topology object
    open_mm_object : 'topology' or 'modeller' OpenMM object, default='topology'
        Untyped OpenMM object to convert to

    Returns
    -------
    open_mm_object : Untyped `topology` or `modeller` object

    """
    openmm_top = app.Topology()

    # Get topology.positions into OpenMM form
    openmm_unit = 1 * simtk_unit.nanometer
    topology.positions.convert_to_units(openmm_unit.unit.get_symbol())
    value = [i.value for i in topology.positions]
    openmm_pos = simtk_unit.Quantity(value=value,
            unit=openmm_unit.unit)

    # Convert bonx information
    box = topology.box
    box.lengths.convert_to_units(u.nanometer)
    lengths = box.lengths.value
    openmm_top.setUnitCellDimensions(lengths)

    # Adding a default chain and residue temporarily
    chain = openmm_top.addChain()
    residue = openmm_top.addResidue(name='RES',
                                    chain=chain)
    # Convert toplogy information
    atom_map = dict()
    for site in topology.sites:
        name = site.name
        element = site.element.name if site.element else None
        # addAtom would add an atom to the openmm Topology
        # AND return that atom
        atom_map[site] = openmm_top.addAtom(name=site.name,
                                            element=element,
                                            residue=residue)

    # Convert bonds information
    for bond in topology.bonds:
        openmm_top.addBond(
            atom1=atom_map[bond.connection_members[0]],
            atom2=atom_map[bond.connection_members[1]])

    # TODO: Figure out how to add residues

    if openmm_object == 'topology':

        return openmm_top

    else:
        modeller = app.Modeller(openmm_top, openmm_pos)

        return modeller


def to_system(topology,
              nonbondedMethod=None,
              nonbondedCutoff=0.8*u.nm,
              switchDistance=0.6*u.nm,
              constraints=None,
              rigidWater=True,
              implicitSolvent=None,
              implicitSolventKappa=None,
              implicitSolventSaltConc=0.0*u.Unit('mol/dm**3'),
              temperature=300*u.Kelvin,
              soluteDielectric=1.0,
              solventDielectric=78,
              removeCMMotion=True,
              hydrogenMass=None,
              ewaldErrorTolerance=0.0005,
              flexibleConstraints=True,
              verbose=False,
              splitDihedrals=False):
    """
    Convert a typed topology object to a typed OpenMM System.  See http://openmm.org for more information.

    Parameters
    ----------
    topology : `Topology` object, default=None
        An untyped topology object.
    nonbondedMethod : cutoff method, optional, default=None
        Cutoff method specified for OpenMM system.  
        Options supported are 'NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'PME', or Ewald objects from simtk.openmm.app.
    nonbondedCutoff : unyt array or float, default=0.8*u.nm
        The nonbonded cutoff must either be a float or a unyt array.  
        Float interpreted in units of nm.
    switchDistance : unyt array or float, default=0.6*u.nm
        The distance at which the switching function is turned on for van der waals interactions.  
        This is ignored when no cutoff is used, and no switch is used if switchDistance is 0, negative, or greater than the cutoff. 
        Float point interpreted in units of nm.
    constraints : 'None', 'app.HBonds', 'app.HAngles', or 'app.AllBonds'
        Type of constraints to add to the System (e.g., SHAKE).
    rigidWater : boolean, optional, default=True
        If True, water is kept rigid regardless of constraint values.
        False value is overriden if constraints is not None
    implicitSolvent : 'None', 'app.HCT', 'app.OBC1', 'app.OBC2', 'app.GBn', 'app.GBn2'. Default=None
        The Generalized Born implicit solvent model to use.
    implicitSolventKappa : 1/distance unyt array or float, default=None
        Debye kappa property related to modeling saltware conditions in GB.  
        It should have units of 1/distance (interpreted as 1/nanometers if no units reported).  
        A value of 'None' means that kappa will be calculated from implicitSolventSaltConc.
    implicitSolventSaltConc : amount/volume unyt array or float, default=0 moles/Liter (moles/dm^3)
        if implicitSolventKappa is 'None', the kappa will be computed from salt concentration.  
        Units should be compatible with mol/L.
    temperature : temperature unyt array or flow, default=300*u.K 
        This is only used to compute kappa from implicitySolventSaltConc.  
        If not unit given, temperature will be interpreted in units of Kelvin.
    soluteDielectric : float, default=1.0
        The dielectric constant of protein interior used in GB.
    solventDielectric : float, default=78.5
        The dielectric constant of water used in GB
    useSASA : boolean, default=False
        If True, use the ACE non-polar solvation model.  
        Otherwise, no SASA-based nonpolar solvation model is used.
    removeCMMotion : boolean, default=True
        If True, the center-of-mass motion will be removed periodically during the simulation.  
        If False, it will not.
    hydrogenMass : mass unyt array or float, default=None
        If not None, hydrogen masses will be changed to this mass and the difference subtracted from the attached heavy atom (hydrogen mass repartitioning).
    ewaldErrorTolerance : float, default=0.0005
        When using PME or Ewald, the Ewald parameters will be calculated from this value.
    flexibleConstraints : boolean, optional, default=True
        If False, the energies and forces from the constrained degrees of freedom will NOT be computed.  
        If True, they will but those degrees of freedom will *still* be constrained).
    verbose : boolean, optional, default=False
        If True, the progress of this subroutine will be printed to stdout.
    splitDihedrals : boolean, optional, default=False
        If True, the dihedrals will be split into two forces -- propers and impropers.  
        This is primarily useful for debugging torsion parameter assignments.

    Returns
    -------
    openmm_system : A typed OpenMM System object

    """

    # TODO: Everything
