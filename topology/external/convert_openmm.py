import mbuild as mb
import unyt as u
import simtk.unit

from topology.core.topology import Topology
from topology.core.site import Site
from topology.utils.sorting import natural_sort

from simtk import openmm


def to_modeller(topology):
    """
    Convert an untyped topology object to an untyped OpenMM modeller.
    This is useful if it's preferred to atom-type a system within
    OpenMM.

    Parameters
    ----------
    topology: topology object
        An untyped topology object
    """
    openmm_top = openmm.app.Topology()
    
    # Get topology.positions into OpenMM form
    top_pos = topology.positions()
    value = [i.value for i in topology.positions()]
    openmm_unit = topology.positions().units
    openmm_pos = simtk.unit.Quantity(value=value,
            unit=openmm_unit)

    # Adding a default chain and residue temporarily
    chain = openmm_top.addChain()
    residue = openmm_top.addResidue(name='RES',
                                    chain=chain)

    #types = list()
    for site in topology.site_list:
        #types.append(site.atom_type.name)
        openmm_top.addAtom(name=site.name,
                           element=site.element.name,
                           residue=residue)

    #unique_types = list(set(types))
    #unique_types.sort(key=natural_sort)
    #for types in unique_types:
    #    res = openmm_top.addResidue(name=types,
    #                          chain=chain)
    # TODO: Figure out how to add residues
    # TODO: Convert connections to OpenMM Bonds


    modeller = openmm.app.Modeller(openmm_top, openmm_pos)

    return modeller


def to_system(topology,
              nonbondedMethod=None,
              nonbondedCutoff=0.8,
              switchDistance=0.6,
              constraints=None,
              rigidWater=True,
              implicitSolvent=None,
              implicitSolventKappa=None,
              implicitSolventSaltConc=0.0,
              temperature=300,
              soluteDielectric=1.0,
              solventDielectric=78,
              removeCMMotion=True,
              hydrogenMass=None,
              ewaldErrorTolerance=0.0005,
              flexibleConstraints=True,
              verbose=False,
              splitDihedrals=False):
    """
    Convert an typed topology object to a typed OpenMM System.

    Parameters
    ----------
    topology: topology object, default=None
        An untyped topology object.
    nonbondedMethod: cutoff method
        Cutoff method specified for OpenMM system.  Options supported
        are 'NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'PME',
        or Ewald objects from simtk.openmm.app.
    nonbondedCutoff: float or distance Quantity
        The nonbonded cutoff must either be a floating point number
        (interpreted as nanometers) or an OpenMM Quantity with
        attached units.
    switchDistance: float or distance Quantity
        The distance at which the switching function is turned on for
        van der waals interactions.  This is ignored when no cutoff is
        used, and no switch is used if switchDistance is 0, negative,
        or greater than the cutoff.
    constraints: 'None', 'app.HBonds', 'app.HAngles', or 'app.AllBonds'
        Type of constraints to add to the System (e.g., SHAKE).
    rigidWater: bool=True
        If True, water is kept rigid regardless of constraint values.
        False value is overriden if constarints is not 'None'.
    implicitSolvent: 'None', 'app.HCT', 'app.OBC1', 'app.OBC2', 'app.GBn', 'app.GBn2'
        The Generalized Born implicit solvent model to use. Default=None
    implicitSolventKappa: float or 1/distance Quantity, default=None
        Debye kappa property related to modeling saltware conditions
        in GB.  It should have units of 1/distance (interpreted as
        1/nanometers if no units reported).  A value of 'None' means
        that kapps will be calculated from implicitSolventSaltConc.
    implicitSolventSaltConc: float or amount/volume Quantity, default=0 moles/Liter
        if implicitSolventKappa is 'None', the kappa will be computed
        from salt concentration.  Units should be compatible with
        mol/L.
    temperature: float or temperature Quantity, default=300 Kelvin
        This is only used to compute kappa from
        implicitySolventSaltConc.
    soluteDielectric: float, default=1.0
        The dielectric constant of protein interior used in GB.
    solventDielectric: float, default=78.5
        The dielectric constant of water used in GB
    useSASA: bool, default=False
        If True, use the ACE non-polar solvation model.  Otherwise, no
        SASA-based nonpolar solvation model is used.
    removeCMMotion: bool, default=True
        If True, the center-of-mass motion will be removed
        periodically during the simulation.  If False, it will not.
    hydrogenMass: float or mass quantity, default=None
        If not None, hydrogen masses will be changed to this mass and
        the difference subtracted from the attached heavy atom
        (hydrogen mass repartitioning).
    ewaldErrorTolerance: float, default=0.0005
        When using PME or Ewald, the Ewald parameters will be
        calculated from this value.
    flexibleConstraints: bool, default=True
        If False, the energies and forces from the constrained degrees
        of freedom will NOT be computed.  If True, they will but those
        degrees of freedom will *still* be constrained).
    verbose: bool, default=False
        If True, the progress of this subroutine will be printed to
        stdout.
    splitDihedrals: bool, default=False
        If True, the dihedrals will be split into two forces --
        propers and impropers.  This is primarily useful for debugging
        torsion parameter assignments.
    """

    # TODO: Everything
