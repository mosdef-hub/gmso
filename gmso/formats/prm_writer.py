"""CHARMM Par format."""

from gmso.core.views import PotentialFilters
from gmso.formats.formats_registry import saves_as
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.compatibility import check_compatibility
from gmso.utils.units import LAMMPS_UnitSystems


def _validate_potential_compatibility(top):
    """Check compatability of topology object potentials with LAMMPSDATA format."""
    pfilter = PotentialFilters.UNIQUE_EXPRESSION
    pot_types = check_compatibility(
        top, _accepted_potentials(), site_pfilter=pfilter, conn_pfilter=pfilter
    )
    return pot_types


def _accepted_potentials():
    """List of accepted potentials that LAMMPS can support."""
    templates = PotentialTemplateLibrary()
    lennard_jones_potential = templates["LennardJonesPotential"]
    lennard_jones_potential.expression /= 4 # no 4*epsilon term
    harmonic_bond_potential = templates["LAMMPSHarmonicBondPotential"]
    harmonic_angle_potential = templates["LAMMPSHarmonicAnglePotential"]
    ub_angle_potential = templates["UreyBradleyAnglePotential"]
    periodic_torsion_potential = templates["PeriodicTorsionPotential"]
    harmonic_improper_potential = templates["LAMMPSHarmonicImproperPotential"]
    accepted_potentialsList = [
        lennard_jones_potential,
        harmonic_bond_potential,
        harmonic_angle_potential,
        ub_angle_potential,
        periodic_torsion_potential,
        harmonic_improper_potential,
    ]
    return accepted_potentialsList


@saves_as(".prm", ".par")
def write_prm(topology, filename, strict_potentials=False):
    """Write CHARMM Parameter .prm file given a parametrized structure.

    Notes
    -----
    Follows format according to
    https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/
    node25.html
    Furthermore, ParmEd can support writing CHARMM par, rtf, str files
    by converting the parmed.Structure into parmed.CharmmParameterSet

    Parmed stores rmin/2 in "rmin"
    """
    # Validation
    # TODO: Use strict_x, (e.g. x=bonds) to validate what topology attrs   to convert
    if not strict_potentials:
        default_parameterMaps = {  # TODO: sites are not checked currently because gmso
            # doesn't store pair potential eqn the same way as the connections.
            "impropers": ["LAMMPSHarmonicImproperPotential"],
            "dihedrals": ["PeriodicTorsionPotential"],
            "angles": ["LAMMPSHarmonicAnglePotential", "UreyBradleyAnglePotential"],
            "bonds": ["LAMMPSHarmonicBondPotential"],
            "sites":["LennardJonesPotential"],
            # "sites":"CoulombicPotential"
        }
        _try_default_potential_conversions(topology, default_parameterMaps)
    _validate_potential_compatibility(topology)

    unit_system = LAMMPS_UnitSystems("real")  # ang, kcal/mol, amu
    # ATOMS
    with open(filename, "w") as f:
        f.write("ATOMS\n")
        unique_atoms = set()
        for site in topology.sites:
            unique_atoms.add(
                (
                    site.atom_type.name,
                    unit_system.convert_parameter(
                        site.atom_type.mass, n_decimals=3, name="mass"
                    ),
                )
            )
        for atom in unique_atoms:
            f.write("MASS -1 {:8s} {:8s}\n".format(atom[0], atom[1]))

        # TODO: Works for harmonic bonds
        f.write("\nBONDS\n")
        for bond in topology.bond_types:
            atom1, atom2 = bond.member_types
            f.write(
                "{:8s} {:8s} {:7s} {:7s}\n".format(
                    atom1,
                    atom2,
                    unit_system.convert_parameter(bond.parameters["k"], n_decimals=3),
                    unit_system.convert_parameter(
                        bond.parameters["r_eq"], n_decimals=3
                    ),
                )
            )

        f.write("\nANGLES\n")
        for angle in topology.angle_types:
            atom1, atom2, atom3 = angle.member_types
            if angle.name == "UreyBradleyAnglePotential":
                f.write(
                    "{:8s} {:8s} {:8s} {:7s} {:7s} {:7s} {:7s}\n".format(
                        atom1,
                        atom2,
                        atom3,
                        unit_system.convert_parameter(
                            angle.parameters["k"], n_decimals=3
                        ),
                        unit_system.convert_parameter(
                            angle.parameters["theta_eq"],
                            n_decimals=3,
                            name="theta_eq",  # necessary for conversion to degrees not radians
                        ),
                        unit_system.convert_parameter(
                            angle.parameters["kub"], n_decimals=3
                        ),
                        unit_system.convert_parameter(
                            angle.parameters["r_eq"], n_decimals=3
                        ),
                    )
                )
            else:  # assume harmonic style:
                f.write(
                    "{:8s} {:8s} {:8s} {:7s} {:7s}\n".format(
                        atom1,
                        atom2,
                        atom3,
                        unit_system.convert_parameter(
                            angle.parameters["k"], n_decimals=3
                        ),
                        unit_system.convert_parameter(
                            angle.parameters["theta_eq"], n_decimals=3, name="theta_eq"
                        ),
                    )
                )

        # These dihedrals need to be PeriodicTorsion Style (Charmm style)
        f.write("\nDIHEDRALS\n")
        for dihedral in topology.dihedral_types:
            # works for PeriodicTorsion
            atom1, atom2, atom3, atom4 = dihedral.member_types
            variable_dtypes = ["k", "n", "phi_eq"]
            zipped_params = (dihedral.parameters[x] for x in variable_dtypes)
            for k, n, phi_eq in zip(*zipped_params):
                f.write(
                    "{:8s} {:8s} {:8s} {:8s} {:7s} {:7s} {:7s}\n".format(
                        atom1,
                        atom2,
                        atom3,
                        atom4,
                        unit_system.convert_parameter(k, n_decimals=3),
                        unit_system.convert_parameter(n, n_decimals=3),
                        unit_system.convert_parameter(
                            phi_eq, n_decimals=3, name="phi_eq"
                        ),
                    )
                )

        f.write("\nIMPROPER\n")
        for improper in topology.improper_types:
            atom1, atom2, atom3, atom4 = improper.member_types
            f.write(
                "{:8s} {:8s} {:8s} {:8s} {:.5s} {:5.3f} {:.5s}\n".format(
                    atom1,
                    atom2,
                    atom3,
                    atom4,
                    unit_system.convert_parameter(
                        improper.parameters["k"], n_decimals=3
                    ),
                    0.0,
                    unit_system.convert_parameter(
                        improper.parameters["phi_eq"], n_decimals=3
                    ),
                )
            )

        f.write("\nNONBONDED\n")
        # nonbonded14 = topology.scaling_factors[0][2]

        for atype in topology.atom_types:
            # atype, 0.0, epsilon, rmin/2, 0.0, epsilon(1-4), rmin/2 (1-4)
            f.write(
                "{:8s} {:8.3f} {:8s} {:8s} {:8.3f} {:8s} {:8s}\n".format(
                    atype.name,
                    0.0,  # ignored,
                    unit_system.convert_parameter(
                        atype.parameters["epsilon"], n_decimals=3
                    ),
                    unit_system.convert_parameter(
                        atype.parameters["sigma"], n_decimals=3
                    ),
                    0.0,
                    unit_system.convert_parameter(
                        atype.parameters["epsilon"], n_decimals=3
                    ),
                    unit_system.convert_parameter(
                        atype.parameters["sigma"], n_decimals=3
                    ),
                )
            )
        # TODO: Add NONBONDED section: NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
        # !adm jr., 5/08/91, suggested cutoff scheme
        f.write("\nEND")


def _try_default_potential_conversions(top, potentialsDict):
    """Take a topology a convert all potentials to the style in potentialDict.""" ""
    for pot_container in potentialsDict:
        containerStr = pot_container[:-1] + "_types"
        if getattr(top, containerStr):
            for potential in potentialsDict[pot_container]:
                top.convert_potential_styles({pot_container: potential})
        elif getattr(top, pot_container):
            raise AttributeError(
                f"Missing parameters in {pot_container} for {top.get_untyped(pot_container)}"
            )
