"""CHARMM Par format."""

import warnings
from gmso.formats.formats_registry import saves_as

__all__ = ["write_prm"]

@saves_as("prm", "par")
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
    # ATOMS
    with open(filename, "w") as f:
        f.write("ATOMS\n")
        unique_atoms = set()
        for site in topology.sites:
            unique_atoms.add((site.atom_type.name, site.mass))
        for atom in unique_atoms:
            f.write("MASS -1 {:8s} {:8.4f}\n".format(atom[0], atom[1]))

        f.write("\nBONDS\n")
        unique_bonds = set()
        for bond in topology.bonds:
            atom1, atom2 = bond.connection_members
            unique_bonds.add(
                (
                    atom1.atom_type.name,
                    atom2.atom_type.name,
                    bond.bond_type,
                )
            )

        for bond in unique_bonds:
            # TODO: only works for harmonic bonds
            f.write(
                "{:8s} {:8s} {:.5f} {:.5f}\n".format(
                    bond[0], bond[1], bond[2].k, bond[2].req
                )
            )

        f.write("\nANGLES\n")
        unique_harmonic_angles = set()
        unique_ubs = set()
        for angle in topology.angles:
            atom1. atom2, atom3 = angle.connection_members
            unique_ubs.add(
                atom1.atom_type.name,
                atom2.atom_type.name,
                atom3.atom_type.name,
                angle.angle_type,
                angle.angle_type.name == "UreyBradleyAnglePotential"
            )

        for angle in unique_angles:
            if angle[4]: # urey_bradley flag
                f.write(
                    "{:8s} {:8s} {:8s} {:.5f} {:.5f} {:.5f} {:.5f}\n".format(
                        angle[0],
                        angle[1],
                        angle[2],
                        angle[3].k,
                        angle[3].theteq,
                        angle[3].kub,
                        angle[3].r_eq,
                    )
            else: # assume harmonic angle
                f.write(
                    "{:8s} {:8s} {:8s} {:.5f} {:.5f}\n".format(
                        angle[0],
                        angle[1],
                        angle[2],
                        angle[3].k,
                        angle[3].theteq,
                    )


        # These dihedrals need to be PeriodicTorsion Style (Charmm style)
        if len(topology.rb_torsions) > 0:
            warnings.warn("RB Torsions detected, but unsupported in par writer")
        f.write("\nDIHEDRALS\n")
        unique_dihedrals = set()
        for dihedral in topology.dihedrals:
            # works for PeriodicTorsion
            atom1, atom2, atom3, atom4 = dihedral.connection_members
            unique_dihedrals.add(
                (
                    atom1.atom_type.name,
                    atom2.atom_type.name,
                    atom3.atom_type.name,
                    atom4.atom_type.name,
                    dihedral.dihedral_type,
                )
            )

        for dihedral in unique_dihedrals:
            f.write(
                "{:8s} {:8s} {:8s} {:8s} {:.5f} {:5d} {:.5f}\n".format(
                    dihedral[0],
                    dihedral[1],
                    dihedral[2],
                    dihedral[3],
                    dihedral[4].k,
                    dihedral[4].n,
                    dihedral[4].phi_eq,
                )
            )

        # TODO: No support for harmonic impropers

        f.write("\nIMPROPER\n")
        unique_impropers = set()
        for improper in topology.impropers:
            atom1, atom2, atom3, atom4 = improper.connection_members
            unique_impropers.add(
                (
                    atom1.atom_type.name,
                    atom2.atom_type.name,
                    atom3.atom_type.name,
                    atom4.atom_type.name,
                    improper.type,
                )
            )
        for improper in unique_impropers:
            # order of impropers goes
            f.write(
                "{:8s} {:8s} {:8s} {:8s} {:.5f} {:5d} {:.5f}\n".format(
                    improper[2],
                    improper[0],
                    improper[1],
                    improper[3],
                    improper[4].psi_k,
                    0,
                    improper[4].psi_eq,
                )
            )

        sc_nb = [a for a in scnb]
        if len(sc_nb) > 1:
            warnings.warn(
                "Multiple 1-4 LJ scalings were detected, "
                "defaulting to first LJ scaling detected, {}".format(sc_nb[0])
            )
            sc_nb = sc_nb[0]
        elif len(sc_nb) == 1:
            sc_nb = sc_nb[0]
        elif len(sc_nb) == 0:
            warnings.warn("No 1-4 LJ scaling was detected, defaulting 1")
            sc_nb = 1.0

        f.write("\nNONBONDED\n")
        unique_atypes = set()
        for atom in topology.atoms:
            unique_atypes.add(atom.atom_type)
        for atype in unique_atypes:
            # atype, 0.0, epsilon, rmin/2, 0.0, epsilon(1-4), rmin/2 (1-4)
            f.write(
                "{:8s} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(
                    atype.name,
                    0.0,
                    -1 * atype.epsilon,
                    atype.rmin,
                    0.0,
                    -1 * sc_nb * atype.epsilon,
                    atype.rmin,
                )
            )

        if topology.has_NBFIX():
            warnings.warn("NBFixes detected but unsupported in par writer")

        f.write("\nEND")
