from __future__ import division

import warnings
import numpy as np
import unyt as u
import datetime

from gmso.utils.testing import allclose

def write_lammpsdata(topology, filename, atom_style='full'):
    """Output a LAMMPS data file.

    Outputs a LAMMPS data file in the 'full' atom style format.
    Assumes use of 'real' units.
    See http://lammps.sandia.gov/doc/atom_style.html for more information on atom styles.

    Parameters
    ----------
    Topology : `Topology`
        A Topology Object
    filename : str
        Path of the output file
    atom_style : str, optional, default='full'
        Defines the style of atoms to be saved in a LAMMPS data file.
        The following atom styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more information on atom styles.

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description of the LAMMPS data format.  
    This is a work in progress, as only atoms, masses, and atom_type information can be written out.

    Some of this function has been adopted from `mdtraj`'s support of the LAMMPSTRJ trajectory format. 
    See https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/lammpstrj.py for details.

    """
    if atom_style not in ['atomic', 'charge', 'molecular', 'full']:
        raise ValueError('Atom style "{}" is invalid or is not currently supported'.format(atom_style))

    # TODO: Support various unit styles

    box = topology.box

    # TODO: charges

    with open(filename, 'w') as data:
        data.write('{} written by topology at {}\n\n'.format(
            topology.name if topology.name is not None else '',
            str(datetime.datetime.now())))
        data.write('{:d} atoms\n'.format(len(topology.sites)))
        if atom_style in ['full', 'molecular']:
            if len(topology.bonds) != 0:
                data.write('{:d} bonds\n'.format(len(topology.bonds)))
            else:
                data.write('0 bonds\n')
            if len(topology.angles) != 0:
                data.write('{:d} angles\n'.format(len(topology.angles)))
            else:
                data.write('0 angles\n')
            if topology.dihedrals != 0:
                data.write('{:d} dihedrals\n\n'.format(len(topology.dihedrals)))
            else:
                data.write('0 dihedrals\n')

        data.write('{:d} atom types\n'.format(len(topology.atom_types)))

        data.write('\n')

        # Box data
        if allclose(box.angles, u.unyt_array([90,90,90],'degree')):
            warnings.warn("Orthorhombic box detected")
            box.lengths.convert_to_units(u.angstrom)
            for i,dim in enumerate(['x', 'y', 'z']):
                data.write('{0:.6f} {1:.6f} {2}lo {2}hi\n'.format(
                    0, box.lengths.value[i],dim))
        else:
            warnings.warn("Non-orthorhombic box detected")
            box.lengths.convert_to_units(u.angstrom)
            box.angles.convert_to_units(u.radian)
            vectors = box.get_vectors()
            a, b, c = box.lengths
            alpha, beta, gamma = box.angles

            lx = a
            xy = b * np.cos(gamma)
            xz = c * np.cos(beta)
            ly = np.sqrt(b**2 - xy**2)
            yz = (b*c*np.cos(alpha) - xy*xz) / ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

            xhi = vectors[0][0]
            yhi = vectors[1][1]
            zhi = vectors[2][2]
            xy = vectors[1][0]
            xz = vectors[2][0]
            yz = vectors[2][1]
            xlo = u.unyt_array(0, xy.units)
            ylo = u.unyt_array(0, xy.units)
            zlo = u.unyt_array(0, xy.units)

            xlo_bound = xlo + u.unyt_array(np.min([0.0, xy, xz, xy+xz]),
                    xy.units)
            xhi_bound = xhi + u.unyt_array(np.max([0.0, xy, xz, xy+xz]),
                    xy.units)
            ylo_bound = ylo + u.unyt_array(np.min([0.0, yz]), xy.units)
            yhi_bound = yhi + u.unyt_array(np.max([0.0, yz]), xy.units)
            zlo_bound = zlo
            zhi_bound = zhi

            data.write('{0:.6f} {1:.6f} xlo xhi\n'.format(
                xlo_bound.value, xhi_bound.value))
            data.write('{0:.6f} {1:.6f} ylo yhi\n'.format(
                ylo_bound.value, yhi_bound.value))
            data.write('{0:.6f} {1:.6f} zlo zhi\n'.format(
                zlo_bound.value, zhi_bound.value))
            data.write('{0:.6f} {1:.6f} {2:.6f} xy xz yz\n'.format(
                xy.value, xz.value, yz.value))

        # Write out mass data
        data.write('\nMasses\n\n')
        for atom_type in topology.atom_types:
            data.write('{:d}\t{:.6f}\t# {}\n'.format(
                topology.atom_types.index(atom_type)+1,
                atom_type.mass.in_units(u.g/u.mol).value,
                atom_type.name
                ))

        # TODO: Get a dictionary of indices and atom types
        if topology.is_typed():
            # TODO: Modified cross-interactions
            # Pair coefficients
            data.write('\nPair Coeffs # lj\n\n')
            for idx, param in enumerate(topology.atom_types):
                data.write('{}\t{:.5f}\t{:.5f}\n'.format(
                    idx+1,
                    param.parameters['epsilon'].in_units(u.Unit('kcal/mol')).value,
                    param.parameters['sigma'].in_units(u.angstrom).value
                    ))

            data.write('\nBond Coeffs\n\n')
            for idx, bond_type in enumerate(topology.bond_types):
                data.write('{}\t{:.5f}\t{:.5f}\n'.format(
                    idx+1,
                    bond_type.parameters['k'].in_units(u.Unit('kcal/mol/angstrom**2')).value,
                    bond_type.parameters['r_eq'].in_units(u.Unit('angstrom')).value
                    ))

            data.write('\nAngle Coeffs\n\n')
            for idx, angle_type in enumerate(topology.angle_types):
                data.write('{}\t{:.5f}\t{:.5f}\n'.format(
                    idx+1,
                    angle_type.parameters['k'].in_units(u.Unit('kcal/mol/degree**2')).value,
                    angle_type.parameters['theta_eq'].in_units(u.Unit('degree')).value
                    ))

            data.write('\nDihedral Coeffs\n')
            for idx, dihedral_type in enumerate(topology.dihedral_types):
                data.write('{}\t{:.5f}\t{:5f}\t{:5f}\t{:.5f}\t{:d}\t{:d}\n'.format(
                    idx+1,
                    dihedral_type.parameters['c0'],
                    dihedral_type.parameters['c1'],
                    dihedral_type.parameters['c2'],
                    dihedral_type.parameters['c3'],
                    0,
                    0
                    ))


        # Atom data
        data.write('\nAtoms\n\n')
        if atom_style == 'atomic':
            atom_line = '{index:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'charge':
            atom_line = '{index:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'molecular':
            atom_line = '{index:d}\t{zero:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'full':
            atom_line ='{index:d}\t{zero:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'

        for i, site in enumerate(topology.sites):
            data.write(atom_line.format(
                index=topology.sites.index(site)+1,
                type_index=topology.atom_types.index(site.atom_type)+1,
                zero=0,charge=0, # TODO: handle charges from atomtype and/or site
                x=site.position[0].in_units(u.angstrom).value,
                y=site.position[1].in_units(u.angstrom).value,
                z=site.position[2].in_units(u.angstrom).value))

        if topology.bonds:
            data.write('\nBonds\n\n')
            for i, bond in enumerate(topology.bonds):
                data.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                i+1,
                topology.bond_types.index(bond.connection_type)+1,
                topology.sites.index(bond.connection_members[0]),
                topology.sites.index(bond.connection_members[1])
                ))

        if topology.angles:
            data.write('\nAngles\n\n')
            for i, angle in enumerate(topology.angles):
                data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                i+1,
                topology.angle_types.index(angle.connection_type)+1,
                topology.sites.index(angle.connection_members[0]),
                topology.sites.index(angle.connection_members[1]),
                topology.sites.index(angle.connection_members[2])
                ))

        if topology.dihedrals:
            data.write('\nDihedrals\n\n')
            for i, dihedral in enumerate(topology.dihedrals):
                data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                i+1,
                topology.dihedral_types.index(dihedral.connection_type)+1,
                topology.sites.index(dihedral.connection_members[0]),
                topology.sites.index(dihedral.connection_members[1]),
                topology.sites.index(dihedral.connection_members[2]),
                topology.sites.index(dihedral.connection_members[3])
                ))
