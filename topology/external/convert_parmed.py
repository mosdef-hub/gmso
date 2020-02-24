import numpy as np
import unyt as u

import topology as topo
from topology.utils.io import import_, has_parmed


if has_parmed:
    pmd = import_('parmed')

def from_parmed(structure):
    """Convert a parmed.Structure to a topology.Topology

    Parameters
    ----------
    structure : parmed.Structure
        parmed.Structure instance that need to be converted

    Returns
    -------
    top : topology.Topology
    """
    msg = ("Provided argument that is not a Parmed Structure")
    assert isinstance(structure, pmd.Structure), msg

    top = topo.Topology(name=structure.title)
    site_map = dict()
    for atom in structure.atoms:
        if isinstance(atom.atom_type, pmd.AtomType):
            atom_type = topo.AtomType(
                name=atom.atom_type.name,
                charge=atom.atom_type.charge * u.elementary_charge,
                parameters={
                    'sigma': (atom.sigma * u.angstrom).in_units(u.nm),
                    'epsilon': atom.epsilon * u.Unit('kcal / mol')
                })
            site = topo.Site(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=([atom.xx, atom.xy, atom.xz] * u.angstrom).in_units(
                    u.nm),
                atom_type=atom_type)
        else:
            site = topo.Site(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=([atom.xx, atom.xy, atom.xz] * u.angstrom).in_units(
                    u.nm),
                atom_type=None)
        site_map[atom] = site
        top.add_site(site)
    top.update_topology()

    if np.all(structure.box):
        # This is if we choose for topology to have abox
        top.box = topo.Box(
            (structure.box[0:3] * u.angstrom).in_units(u.nm),
            angles=u.degree * structure.box[3:6])

    for bond in structure.bonds:
        # Generate bond parameters for BondType that gets passed
        # to Bond
        if isinstance(bond.type, pmd.BondType):
            bond_params = {
                'k': (2 * bond.type.k * u.Unit('kcal / (nm**2 * mol)')),
                'r_eq': (bond.type.req * u.angstrom).in_units(u.nm)
            }
            new_connection_type = topo.BondType(parameters=bond_params,
                    member_types=[bond.atom1.type, bond.atom2.type])
            top_connection = topo.Bond(connection_members=[site_map[bond.atom1],
                site_map[bond.atom2]],
                connection_type=new_connection_type)

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = topo.Bond(connection_members=[site_map[bond.atom1],
                site_map[bond.atom2]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)
    top.update_topology()
    print(top.n_bonds)

    for angle in structure.angles:
        # Generate angle parameters for AngleType that gets passed
        # to Angle
        if isinstance(angle.type, pmd.AngleType):
            angle_params = {
                'k': (2 * angle.type.k * u.Unit('kcal / (rad**2 * mol)')),
                'theta_eq': (angle.type.theteq * u.degree)
            }
            new_connection_type = topo.AngleType(parameters=angle_params,
                    member_types=[angle.atom1.type, angle.atom2.type, 
                                angle.atom3.type])
            top_connection = topo.Angle(connection_members=[site_map[angle.atom1],
                site_map[angle.atom2], site_map[angle.atom3]],
                connection_type=new_connection_type)

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = topo.Angle(connection_members=[site_map[angle.atom1],
                site_map[angle.atom2], site_map[angle.atom3]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)
    top.update_topology()

    for dihedral in structure.dihedrals:
        # Generate dihedral parameters for DihedralType that gets passed
        # to Dihedral
        # These all follow periodic torsions functions
        # (even if they are improper dihedrals)
        # Which are the default expression in top.DihedralType
        # These periodic torsion dihedrals get stored in top.dihedrals
        if dihedral.improper:
            warnings.warn("ParmEd improper dihedral {} ".format(dihedral) +
                    "following periodic torsion " +
                    "expression detected, currently accounted for as " +
                    "topology.Dihedral with a periodic torsion expression")
        if isinstance(dihedral.type, pmd.DihedralType):
            dihedral_params = {
                'k': (dihedral.type.phi_k * u.Unit('kcal / mol')),
                'phi_eq': (dihedral.type.phase * u.degree),
                'n': dihedral.type.per * u.dimensionless
            }
            new_connection_type = topo.DihedralType(parameters=dihedral_params,
                    member_types=[dihedral.atom1.type, dihedral.atom2.type,
                                dihedral.atom3.type, dihedral.atom4.type])
            top_connection = topo.Dihedral(connection_members=
                    [site_map[dihedral.atom1], site_map[dihedral.atom2], 
                        site_map[dihedral.atom3], site_map[dihedral.atom4]],
                connection_type=new_connection_type)

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = topo.Dihedral(connection_members=
                    [site_map[dihedral.atom1], site_map[dihedral.atom2], 
                        site_map[dihedral.atom3], site_map[dihedral.atom4]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)

    for rb_torsion in structure.rb_torsions:
        # Generate dihedral parameters for DihedralType that gets passed
        # to Dihedral
        # These all follow RB torsion functions
        # These RB torsion dihedrals get stored in top.dihedrals
        if rb_torsion.improper:
            warnings.warn("ParmEd improper dihedral {} ".format(rb_torsion) +
                    "following RB torsion " +
                    "expression detected, currently accounted for as " +
                    "topology.Dihedral with a RB torsion expression")
        if isinstance(rb_torsion.type, pmd.RBTorsionType):
            dihedral_params = {
                'c0': (rb_torsion.type.c0 * u.Unit('kcal/mol')),
                'c1': (rb_torsion.type.c1 * u.Unit('kcal/mol')),
                'c2': (rb_torsion.type.c2 * u.Unit('kcal/mol')),
                'c3': (rb_torsion.type.c3 * u.Unit('kcal/mol')),
                'c4': (rb_torsion.type.c4 * u.Unit('kcal/mol')),
                'c5': (rb_torsion.type.c5 * u.Unit('kcal/mol')),
            }
            new_connection_type = topo.DihedralType(parameters=dihedral_params,
                    expression='c0 * cos(phi)**0 + c1 * cos(phi)**1 + ' +
                    'c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + ' +
                    'c5 * cos(phi)**5',
                    independent_variables='phi',
                    member_types=[rb_torsion.atom1.type, rb_torsion.atom2.type,
                                rb_torsion.atom3.type, rb_torsion.atom4.type])
            top_connection = topo.Dihedral(connection_members=
                    [site_map[rb_torsion.atom1], site_map[rb_torsion.atom2], 
                        site_map[rb_torsion.atom3], site_map[rb_torsion.atom4]],
                connection_type=new_connection_type)

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = topo.Dihedral(connection_members=
                    [site_map[rb_torsion.atom1], site_map[rb_torsion.atom2], 
                        site_map[rb_torsion.atom3], site_map[rb_torsion.atom4]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)

    top.update_top()


    return top
