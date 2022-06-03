"""Generic utilities for parameterizing a gmso Topology."""

from gmso.core.angle import Angle
from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.dihedral import Dihedral
from gmso.core.improper import Improper

POTENTIAL_GROUPS = {
    Bond: "bond_type",
    Angle: "angle_type",
    Dihedral: "dihedral_type",
    Improper: "improper_type",
    Atom: "atom_type",
}
