from warnings import warn
import unyt as u
import numpy as np

import gmso
from gmso.utils.io import import_, has_openmm, has_simtk_unit


if has_openmm & has_simtk_unit:
    simtk_unit = import_('simtk.unit')
    from simtk import openmm as mm
    from simtk.openmm.app import *
    from simtk.openmm import *

def from_openmm(openmm_object, system=None):
    """Convert an openmm topology object into a gmso.Topology

    Convert an OpenMM object, either Topology or Modeller,
    to a gmso.Topology. User has the option to refer types
    information from a openmm System.
    Mapping plan:
    OpenMM: Topology - Chains - Residues - Atoms
    GMSO: Topology - SubTopology - SubTopology - Sites

    Parameters
    ----------
    openmm_object : simtk.openmm.app.Topology or
                    simtk.openmm.app.Modeller
        OpenMM Topology or Modeller object
        that need to be converted.
    system : 'simtk.openmm.System', optional, default=None
        Referenced OpenMM system, if none is supplied,
        the topology will be returned untyped.

    Return
    ------
    top : gmso.Topology
        Typed or untyped GMSO Topology object.
    """
    msg = 'Given object is not an OpenMM Topology or \
           OpenMM Modeller'

    # Separate Top and positions information
    # Assuming positions unit to be nm, may need to actually
    # parse the unit from OpenMM Modeller
    if isinstance(openmm_object, openmm.Modeller):
        mm_top = openmm_object.topology
        mm_pos = openmm_object.positions._value * u.nm
    elif isinstance(openmm_object, openmm.Topology):
        # If this is a Topology, all position is set to 0
        mm_top = openmm_object
        mm_pos = [np.zeros(3) for i in
                  range(len(mm_top.getNumAtoms))]
    else:
        raise TypeError(msg)


    # Initialize GMSO Topology
    top = gmso.Topology()
    top.name = 'Topology'

    # Convert box information
    if not system:
        lengths = mm_top.getPeriodicBoxVectors()[0] * u.nm
        angles = mm_top.getPeriodicBoxVectors()[1] * u.degree
    else:
        lenghts = # from system
        angles = # from system
        # Give a warning about the system override the
        # box informatoin from Topology/Modeller
    top.box = gmso.Box(lengths=lengths, angles=angles)

    # Convert topology information
    site_map = dict() # mapping atom id-> site
    for chain in mm_top.chains():
        chain_name = chain.name + chain.id
        gmso_chain = gmso.SubTopology(name=chain_name,
                                parent=top)
        for residue in chain.residues():
            residue_name = residue.name + residue.id
            gmso_residue = gmso.SubTopology(name=res_name,
                                        parent=gmso_chain)
            for atom in residue.atoms():
                pos = mm_pos.pop(0)._value
                name = atom.name + '_' + atom.id
                site = gmso.Site(name=name, pos= pos)
                gmso_residue.add_site(site)
                site_map[atom.id] = site
        top.add_subtop(chain)

    # Convert bonds information
    for bond in openmm_topology.bonds():
        top_connection = gmso.Bond(
        connection_members=[site_map[bond.atom1.id],
                            site_map[bond.atom2.id]])
        top.add(top_connection)

    # Checkpoint for barebone, untyped GMSO Top
    if not system:
        pass
    else:
        # Call helper function to apply forces system to top
        top = apply_system(top, system, site_map)

    return top

def apply_system(top, system, site_map=None):
    """ Helper function to apply OpenMM System to GMSO Topology

    Applying System information to relevant Topology object.
    Only Non-bonded, Bond, Angle, and Dihedral Forces will
    be considered/translated, why information about
    Thermometers and Barometers in OpenMM System will be
    discarded. The Thermometers and Barometers are still
    stored in a set, so this can still be translate to
    relevant GMSO Topology variable in the future.

    Paramters
    ---------
    top : gmso.Topology
        Host GMSO Topology object. If site_dict is not
        provided, site's name must follow specific convention,
        `element_id` (element name separated by OpenMM id).
    system : 'simtk.openmm.System'
        The System from which we want to translated the
        forces information.
    site_map : dict, optional, default=None
        Dictionary of sites with key is OpenMM atom index.
    Return
    ------
    top : gmso.Topology
        Typed GMSO Topology object
    """
    # Sanity checks
    # May add method to load system from file
    # (Meaning system can be str=path/to/file)
    msg1 = 'Given system is not an OpenMM System'
    assert isinstance(system, mm.System), msg
    msg2 = 'Topology and System have different \
            number of atoms'
    assert len(struc.atomts) == system.getNumAtoms(), msg3

    # At this point, GMSO Topology only support potential
    # forces for atoms (non-bonded force), bond, angle,
    # and dihedrals (torsion). Other forces of the
    # OpenMM System will be disregarded.
    # Forces include information about bonds forces,
    # angle forces, dihedral forces, and non-bonded forces
    atom_forces = list()
    bond_forces = list()
    angle_forces = list()
    dihedral_forces = list()

    for force in system.getForces():
        if 'Nonbonded' in str(force):
            atom_forces.append(force):
        elif 'Bond' in str(force):
            bond_forces.append(force)
        elif 'Angle' in str(force):
            angle_forces.append(force)
        elif 'Torsion' in str(force):
            dihedral_forces.append(force)
        else:
            warn('OpenMM System {} will be \
                  disregarded').format(str(force))

    # Rebuild site_map
    if site_map:
        pass
    else:
        site_map = dict()
        for site in top.sites:
            id = int(site.name.split('_')[1)]
            # Probably need to add a check here
            # and gives better error if id is not int
            site_map[id] = site

    for atom_force in atom_forces:
        # To do
    for bond_force in bond_forces:
        # To do
    for angle_force in angle_forces:
        # To do
    for dihedral_force in dihedral_forces:
        # To do, to do, to do, to doooooo
        '''
                                  .--.            .--.
                                 ( (`\\."--``--".//`) )
                                  '-.   __   __    .-'
                                   /   /__\ /__\   \
                                  |    \ 0/ \ 0/    |
                                  \     `/   \`     /
                                   `-.  /-"""-\  .-`
                                     /  '.___.'  \
                                     \     I     /
                                      `;--'`'--;`
                                jgs     '.___.'
        '''
    return top

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
