"""Convert to and from an OpenMM Topology or System object."""

import unyt as u

from gmso.utils.io import has_openmm, has_openmm_unit, import_

if has_openmm & has_openmm_unit:
    openmm_unit = import_("openmm.unit")
    from openmm import *
    from openmm.app import *


def to_openmm(topology, openmm_object="topology"):
    """Convert an untyped topology object to an untyped OpenMM modeller or topology.

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
    topology.positions.convert_to_units(u.nm)
    value = [i.value for i in topology.positions]
    openmm_pos = openmm_unit.Quantity(value=value, unit=openmm_unit.nanometer)

    # Adding a default chain and residue temporarily
    chain = openmm_top.addChain()
    residue = openmm_top.addResidue(name="RES", chain=chain)

    for site in topology.sites:
        openmm_top.addAtom(
            name=site.name, element=site.element.name, residue=residue
        )

    # Set box
    box = topology.box
    box.lengths.convert_to_units(u.nanometer)
    lengths = box.lengths.value
    openmm_top.setUnitCellDimensions(lengths)

    # TODO: Figure out how to add residues
    # TODO: Convert connections to OpenMM Bonds

    if openmm_object == "topology":
        return openmm_top

    else:
        modeller = app.Modeller(openmm_top, openmm_pos)

        return modeller


def to_system(
    topology,
    nonbondedMethod=None,
    nonbondedCutoff=0.8 * u.nm,
    switchDistance=0.6 * u.nm,
    constraints=None,
    rigidWater=True,
    implicitSolvent=None,
    implicitSolventKappa=None,
    implicitSolventSaltConc=0.0 * u.Unit("mol/dm**3"),
    temperature=300 * u.Kelvin,
    soluteDielectric=1.0,
    solventDielectric=78,
    removeCMMotion=True,
    hydrogenMass=None,
    ewaldErrorTolerance=0.0005,
    flexibleConstraints=True,
    verbose=False,
    splitDihedrals=False,
):
    """
    Convert a typed topology object to a typed OpenMM System.  See http://openmm.org for more information.

    Parameters
    ----------
    topology : `Topology` object, default=None
        An untyped topology object.
    nonbondedMethod : cutoff method, optional, default=None
        Cutoff method specified for OpenMM system.
        Options supported are 'NoCutoff', 'CutoffNonPeriodic', 'CutoffPeriodic', 'PME', or Ewald objects from openmm.app.
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
