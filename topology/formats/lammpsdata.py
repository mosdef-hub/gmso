from __future__ import division

import warnings
import numpy as np
import unyt as u
import datetime

from topology.utils.sorting import natural_sort
from topology.core.site import Site
from topology.core.atom_type import AtomType
from topology.utils.testing import allclose
from topology.core.topology import Topology
from topology.core.box import Box


def read_lammpsdata(filename, atom_style='full'):
    top = Topology()
    u.define_unit("kcal", u.cal * 1000)

    with open(filename, 'r') as lammps_file:
        top.name = str(lammps_file.readline().strip())
        lammps_file.readline()
        n_atoms = int(lammps_file.readline().split()[0])
        coords = u.nm * np.zeros(shape=(n_atoms, 3))
        n_bonds = int(lammps_file.readline().split()[0])
        n_angles = int(lammps_file.readline().split()[0])
        n_dihedrals = int(lammps_file.readline().split()[0])
        lammps_file.readline()
        n_atomtypes = int(lammps_file.readline().split()[0])

        for i in range(4):
            n_gen_types = lammps_file.readline()
            if n_gen_types == '\n':
                print("No more types to track")
                break
            elif n_gen_types.split()[1] == 'bond':
                n_bondtypes = int(n_gen_types.split()[0])
            elif n_gen_types.split()[1] == 'angle':
                n_angletypes = int(n_gen_types.split()[0])
            elif n_gen_types.split()[1] == 'dihedral':
                n_dihedraltypes = int(n_gen_types.split()[0])

        x_line = lammps_file.readline().split()
        y_line = lammps_file.readline().split()
        z_line = lammps_file.readline().split()

        x =  float(x_line[1])-float(x_line[0])
        y =  float(y_line[1])-float(y_line[0])
        z =  float(z_line[1])-float(z_line[0])

        # Box Information
        top.box = Box([x*u.angstrom,y*u.angstrom,z*u.angstrom])


    unique_types = _get_masses(filename,n_atomtypes)
    #atoms = _get_atoms(filename, n_atoms, coords)
    charge_dict, coords_dict, type_dict = _get_atoms(filename, n_atoms, coords)
    sigma_dict, epsilon_dict = _get_pairs(filename, n_atomtypes)

    #type_ids = [typ for typ in atoms[0]]
    #charges = [charge for charge in atoms[1]]
    #coordinates = [coord for coord in atoms[2]]

    #for type_id, charge, coord, ep, sigma in zip(type_ids, charges,coordinates,epsilons,sigmas):
    for k, v in type_dict.items():
        atomtype = AtomType(name=k,
                mass=unique_types[k],
                charge=charge_dict[k],
                parameters={
                    'sigma': 2*u.angstrom,
                    'epsilon':2*u.kcal/u.mol})
                #parameters={
                #    'sigma': sigma_dict[k],
                #    'epsilon': epsilon_dict[k]}
                #)
        for i in range(v):
            site = Site(name="atom",
                position=coords_dict[k][i],
                atom_type=atomtype
                )
            top.add_site(site, update_types=False)
            print('{}:{}'.format(k, i))

    return top


def write_lammpsdata(topology, filename, atom_style='full'):
    """Output a LAMMPS data file.

    Outputs a LAMMPS data file in the 'full' atom style format. Assumes use
    of 'real' units. See http://lammps.sandia.gov/doc/atom_style.html for
    more information on atom styles.

    Parameters
    ----------
    Topology : Topology
        Topology Object
    filename : str
        Path of the output file
    atom_style: str
        Defines the style of atoms to be saved in a LAMMPS data file. The following atom
        styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more
        information on atom styles.

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description
    of the LAMMPS data format. Currently the following sections are supported (in
    addition to the header): *Masses*, *Nonbond Coeffs*, *Bond Coeffs*, *Angle
    Coeffs*, *Dihedral Coeffs*, *Atoms*, *Bonds*, *Angles*, *Dihedrals*

    Some of this function has beed adopted from `mdtraj`'s support of the LAMMPSTRJ
    trajectory format. See https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/lammpstrj.py for details.

    """
    if atom_style not in ['atomic', 'charge', 'molecular', 'full']:
        raise ValueError('Atom style "{}" is invalid or is not currently supported'.format(atom_style))

    xyz = list()
    types = list()
    for site in topology.sites:
        xyz.append([site.position[0],site.position[1],site.position[2]])
        types.append(site.atom_type.name)

    forcefield = True
    if topology.sites[0].atom_type.name in ['', None]:
        forcefield = False

    box = topology.box

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    # TODO: charges
    # TODO: bonds
    # TODO: Angles
    # TODO: Dihedrals

    # TODO: Figure out handling bond, angle, and dihedral indices

    # placeholder; change later
    bonds = 0
    angles = 0
    dihedrals = 0

    with open(filename, 'w') as data:
        data.write('{} written by topology at {}\n\n'.format(
            topology.name if topology.name is not None else '',
            str(datetime.datetime.now())))
        data.write('{:d} atoms\n'.format(len(topology.sites)))
        if atom_style in ['full', 'molecular']:
            if bonds != 0:
                data.write('{:} bonds\n'.format(len(bonds)))
            else:
                data.write('0 bonds\n')
            if angles != 0:
                data.write('{:} angles\n'.format(len(angles)))
            else:
                data.write('0 angles\n')
            if dihedrals != 0:
                data.write('{:} dihedrals\n\n'.format(len(dihedrals)))
            else:
                data.write('0 dihedrals\n\n')

        data.write('{:d} atom types\n'.format(len(set(types))))

        # TODO: Write out bonds, angles, and dihedrals

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

        # Mass data
        masses = [site.atom_type.mass for site in topology.sites]
        mass_dict = dict([(unique_types.index(atom_type)+1,mass) for atom_type,mass in zip(types,masses)])

        data.write('\nMasses\n\n')
        for atom_type,mass in mass_dict.items():
            data.write('{:d}\t{:.6f}\t# {}\n'.format(
                atom_type,
                mass.in_units(u.g/u.mol).value,
                unique_types[atom_type-1]))
        if forcefield:
            sigmas = [site.atom_type.parameters['sigma'] for site in topology.sites]
            epsilons = [site.atom_type.parameters['epsilon'] for site in topology.sites]
            sigma_dict = dict([(unique_types.index(atom_type)+1,sigma) for atom_type,sigma in zip(types,sigmas)])
            epsilon_dict = dict([(unique_types.index(atom_type)+1,epsilon) for atom_type,epsilon in zip(types,epsilons)])

            # TODO: Modified cross-interactions
            # Pair coefficients
            data.write('\nPair Coeffs # lj\n\n')
            for idx,epsilon in epsilon_dict.items():
                data.write('{}\t{:.5f}\t{:.5f}\n'.format(
                    idx,
                    epsilon.in_units(u.Unit('kcal')).value,
                    sigma_dict[idx].in_units(u.angstrom).value))

        # TODO: Write out bond coefficients
        # TODO: Write out angle coefficients
        # TODO: Write out dihedral coefficients

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


        # TODO: Add back in correct 'type_index' and 'charge'
        #for i,coords in enumerate(xyz):
        #    data.write(atom_line.format(
        #        index=i+1,type_index=unique_types.index(types[i])+1,
        #        zero=0,charge=charges[i],
        #        x=coords[0],y=coords[1],z=coords[2]))

        for i,coords in enumerate(xyz):
            data.write(atom_line.format(
                index=i+1,type_index=unique_types.index(types[i])+1,
                zero=0,charge=0, # TODO: handle charges from atomtype and/or site
                x=coords[0].in_units(u.angstrom).value,
                y=coords[1].in_units(u.angstrom).value,
                z=coords[2].in_units(u.angstrom).value))

        # TODO: Write out bonds
        # TODO: Write out angles
        # TODO: Write out dihedrals


def _get_masses(filename, n_atomtypes):
    with open(filename, 'r') as lammps_file:
        types = dict()
        for line in lammps_file:
            if 'Masses' in line:
                lammps_file.readline()
                for i in range(n_atomtypes):
                    type_line = lammps_file.readline().split()
                    types[type_line[0]] = float(type_line[1]) * u.g
    
    return types


def _get_atoms(filename, n_atoms, coords):
    #types = list()
    #charges = list()
    #coord_list = list()
    type_dict = dict()
    charge_dict = dict()
    coord_dict = dict()
    with open(filename, 'r') as lammps_file:
        for line in lammps_file:
            if 'Atoms' in line:
                lammps_file.readline()
                for i in range(n_atoms):
                    atom = lammps_file.readline().split()
                    atom_type = atom[2]
                    charge = u.unyt_quantity(float(atom[3]), u.elementary_charge)
                    coords[i] = u.nm * np.array([
                        float(atom[4]),
                        float(atom[5]),
                        float(atom[6])])

                    if atom_type in type_dict:
                        type_dict[atom_type] += 1
                    else:
                        charge_dict[atom_type] = charge
                        coord_dict[atom_type] = coords
                        type_dict[atom_type] = 1

    #return types, charges, coord_list
    return charge_dict, coord_dict, type_dict


def _get_pairs(filename, n_atomtypes):
    #types = list()
    #epsilon_list = list()
    #sigma_list = list()
    sigma_dict = dict()
    epsilon_dict = dict()
    with open(filename, 'r') as lammps_file:
        for line in lammps_file:
            if 'Pair' in line:
                lammps_file.readline()
                for i in range(n_atomtypes):
                    atom_type = lammps_file.readline().split()
                    epsilon = float(atom_type[1]) * u.kcal / u.mol
                    sigma = float(atom_type[2]) * u.angstrom

                    sigma_dict[atom_type[0]] = sigma
                    epsilon_dict[atom_type[0]] = epsilon
                    #types.append(atom_type[0])
                    #epsilon_list.append(epsilon)
                    #sigma_list.append(sigma)

    return sigma_dict, epsilon_dict
