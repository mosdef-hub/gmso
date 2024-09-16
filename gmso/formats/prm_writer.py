"""CHARMM Par format."""

from gmso.formats.formats_registry import saves_as
from gmso.utils.units import LAMMPS_UnitSystems


@saves_as(".prm", ".par")
def write_prm(topology, filename):
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
                        site.atom_type.mass, n_decimals=6, name="mass"
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
                "{:8s} {:8s} {:.5f} {:.5f}\n".format(
                    atom1, atom2, bond.parameters["k"], bond.parameters["r_eq"]
                )
            )

        f.write("\nANGLES\n")
        for angle in topology.angle_types:
            atom1, atom2, atom3 = angle.member_types
            if angle.name == "UreyBradleyAnglePotential":
                f.write(
                    "{:8s} {:8s} {:8s} {:.5f} {:.5f} {:.5f} {:.5f}\n".format(
                        atom1,
                        atom2,
                        atom3,
                        angle.parameters["k"],
                        angle.parameters["theta_eq"],
                        angle.parameters["kub"],
                        angle.parameters["r_eq"],
                    )
                )
            else:  # assume harmonic style:
                f.write(
                    "{:8s} {:8s} {:8s} {:.5f} {:.5f}\n".format(
                        atom1,
                        atom2,
                        atom3,
                        angle.parameters["k"],
                        angle.parameters["theta_eq"],
                    )
                )

        # These dihedrals need to be PeriodicTorsion Style (Charmm style)
        f.write("\nDIHEDRALS\n")
        for dihedral in topology.dihedral_types:
            # works for PeriodicTorsion
            atom1, atom2, atom3, atom4 = dihedral.member_types
            f.write(
                "{:8s} {:8s} {:8s} {:8s} {:.5f} {:5f} {:.5f}\n".format(
                    atom1,
                    atom2,
                    atom3,
                    atom4,
                    dihedral.parameters["k"][0],
                    dihedral.parameters["n"][0],
                    dihedral.parameters["phi_eq"][0],
                )
            )

        # TODO: No support for harmonic impropers

        f.write("\nIMPROPER\n")
        for improper in topology.improper_types:
            atom1, atom2, atom3, atom4 = improper.member_types
            f.write(
                "{:8s} {:8s} {:8s} {:8s} {:.5f} {:5f} {:.5f}\n".format(
                    atom1,
                    atom2,
                    atom3,
                    atom4,
                    improper.parameters["phi_k"][0],
                    improper.parameters["n"][0],
                    improper.parameters["phi_eq"][0],
                )
            )

        f.write("\nNONBONDED\n")
        nonbonded14 = topology.scaling_factors[0][2]

        for atype in topology.atom_types:
            # atype, 0.0, epsilon, rmin/2, 0.0, epsilon(1-4), rmin/2 (1-4)
            f.write(
                "{:8s} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(
                    atype.name,
                    0.0,  # ignored,
                    atype.parameters["epsilon"],
                    atype.parameters["sigma"],
                    0.0,
                    atype.parameters["epsilon"] * nonbonded14,
                    atype.parameters["sigma"] * nonbonded14,
                )
            )
        f.write("\nEND")
