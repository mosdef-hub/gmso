from __future__ import division

import warnings
import numpy as np
import unyt as u
import datetime

from gmso.utils.sorting import natural_sort
from gmso.utils.testing import allclose
from gmso.core.site import Site
from gmso.core.atom_type import AtomType
from gmso.core.bond_type import BondType
from gmso.core.angle_type import AngleType
from gmso.core.bond import Bond
from gmso.core.angle import Angle
from gmso.core.topology import Topology
from gmso.core.box import Box
from gmso.core.element import element_by_mass


def read_lammpsdata(filename, atom_style='full', potential='lj'):
    """
    Read in a lammps data file as a GMSO topology

    Parameters
    ----------
    filename : str
        LAMMPS data file
    atom_style : str, optional, default='full'
        Inferred atom style defined by LAMMPS
    potential: str, optional, default='lj'
        Potential type defined in data file

    Returns
    -------
    top : GMSO Topology
        A GMSO Topology object

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description of the LAMMPS data format.  
    Currently only supporting the 'full' atom style.  
    Currently only supporting LJ potential parameters.
    """
    # TODO: Add argument to ask if user wants to infer bond type
    top = Topology()

    _get_box_coordinates(filename, top)
    top, type_list = _get_ff_information(filename, top)
    _get_atoms(filename, top, type_list)
    _get_connection(filename,top, connection_type='bond')
    _get_connection(filename, top, connection_type='angle')

def _get_connection(filename, topology, connection_type):
    """
    General function to parse connection types
    """
    with open(filename, 'r') as lammps_file:
        for i, line in enumerate(lammps_file):
            if connection_type in line.split():
                n_connection_types = int(line.split()[0])
            if connection_type.capitalize() in line.split():
                break
    connection_type_lines = open(filename, 'r').readlines()[i+2:i+n_connection_types+2]
    connection_type_list = list()
    for line in connection_type_lines:
        if connection_type == 'bond':
            c_type = BondType(name=line.split()[0]
                    )
            c_type.parameters['k']=float(line.split()[1])*u.Unit(
                                     'kcal/mol/angstrom**2')*2
            c_type.parameters['r_eq']=float(line.split()[2])*u.angstrom
        elif connection_type == 'angle':
            c_type = AngleType(name=line.split()[0]
                    )
            c_type.parameters['k']=float(line.split()[1])*u.Unit(
                                     'kcal/mol/radian**2')*2
            c_type.parameters['theta_eq']=float(line.split()[2])*u.radian

        connection_type_list.append(c_type)

    with open(filename, 'r') as lammps_file:
        for i, line in enumerate(lammps_file):
            if connection_type + 's' in line.split():
                n_connections = int(line.split()[0])
            if connection_type.capitalize() + 's' in line.split():
                break
    connection_lines = open(filename, 'r').readlines()[i+2:i+n_connections+2]
    # Determine number of sites to generate
    if connection_type == 'bond':
        n_sites = 2
    elif connection_type == 'angle':
        n_sites = 3
    else:
        n_sites = 4
    for i, line in enumerate(connection_lines):
        site_list = list()
        for j in range(n_sites):
            site = topology.sites[int(line.split()[j+2])-1]
            site_list.append(site)
        if connection_type == 'bond':
            connection = Bond(
                connection_members=site_list,
                connection_type=connection_type_list[int(line.split()[1])-1],
                    )
        elif connection_type == 'angle':
            connection = Angle(
                connection_members=site_list,
                connection_type=connection_type_list[int(line.split()[1])-1],
                    )
        topology.add_connection(connection)

    return topology

def _get_atoms(filename, topology, type_list):
    with open(filename, 'r') as lammps_file:
        for i, line in enumerate(lammps_file):
            if 'atoms' in line.split():
                n_atoms = int(line.split()[0])
            if 'Atoms' in line.split():
                break
    atom_lines = open(filename, 'r').readlines()[i+2:i+n_atoms+2]
    for line in atom_lines:
        atom = line.split()
        atom_type = atom[2]
        charge = u.unyt_quantity(float(atom[3]), u.elementary_charge)
        coord = u.angstrom * u.unyt_array([
            float(atom[4]),
            float(atom[5]),
            float(atom[6])])
        site = Site(
            charge=charge,
            position=coord,
            atom_type=type_list[int(atom[2])-1]
            )
        element = element_by_mass(site.atom_type.mass)
        site.name = element.name
        site.element = element
        topology.add_site(site)

    topology.update_sites()

    return topology

def _get_box_coordinates(filename, topology):
    with open(filename, 'r') as lammps_file:
        for line in lammps_file:
            if 'xlo' in line.split():
                break
        x_line = line.split()
        y_line = lammps_file.readline().split()
        z_line = lammps_file.readline().split()

        x =  float(x_line[1])-float(x_line[0])
        y =  float(y_line[1])-float(y_line[0])
        z =  float(z_line[1])-float(z_line[0])

        # Box Information
        lengths = u.unyt_array([x,y,z], u.angstrom)
        topology.box = Box(lengths)

        return topology

def _get_ff_information(filename, topology):
    with open(filename, 'r') as lammps_file:
        for i, line in enumerate(lammps_file):
            if 'atom' in line:
                n_atomtypes = int(line.split()[0])
            elif 'Masses' in line:
                break
    mass_lines = open(filename, 'r').readlines()[i+2:i+n_atomtypes+2]
    type_list = list()
    for line in mass_lines:
        atom_type = AtomType(name=line.split()[0],
                             mass=float(line.split()[1])
                             )
        type_list.append(atom_type)

    with open(filename, 'r') as lammps_file:
        for i, line in enumerate(lammps_file):
            if 'Pair' in line:
                break
    # Need to figure out if we're going have mixing rules printed out
    # Currently only reading in LJ params
    pair_lines = open(filename, 'r').readlines()[i+2:i+n_atomtypes+2]
    for i, pair in enumerate(pair_lines):
        if len(pair.split()) == 3:
            type_list[i].parameters['sigma'] = float(
                    pair.split()[2]) * u.angstrom
            type_list[i].parameters['epsilon'] = float(
                    pair.split()[1]) * (u.kcal/u.mol)
        elif len(pair.split()) == 4:
            raise warnings.warn('Currently not reading in mixing rules')

    return topology, type_list

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
                    epsilon.in_units(u.Unit('kcal/mol')).value,
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
