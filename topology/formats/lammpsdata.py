from __future__ import division

from collections import OrderedDict
from warnings import warn
import itertools as it

import numpy as np
import unyt as u

from topology.core.topology import Topology
from topology.core.box import Box
from topology.core.site import Site
from topology.utils.sorting import natural_sort
from topology.testing.utils import allclose


def write_lammpsdata(topology, filename, atom_style='full',
        nbfix_in_data_file=True):
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

    xyz = np.array([[site.position[0],site.position[1],site.position[2]] for site in
        topology.site_list])

    forcefield = True
    if topology.site_list[0].atom_type.name in ['', None]:
        forcefield = False

    box = topology.box

    if forcefield:
        types = [site.atom_type.name for site in topology.site_list]
        unique_types = list(set(types))
        unique_types.sort(key=natural_sort)
    else:
        # TODO: I think we want to be able to write out untyped
        # topologies
        # types = [site.name for site in topology.site_list]
        pass

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    #charges = [site.charge for site in topology.site_list]
    #site_with_index = [i for i in topology.site_list]
    #for i,site in enumerate(site_with_index):
    #    site.index = i

    # bonds = [[bond.site1, bond.site2] for bond in topology.connection_list]
    # [(i, j.position) for i,j in enumerate(top.site_list)]
    # TODO: Angles
    # angles 
    # TODO: Dihedrals
    # dihedrals = 

    # TODO: Figure out handling bond, angle, and dihedral indices
    #if topology.n_connections == 0:
    #    bond_types = np.ones(len(bonds), dtype=int)
    #else:
    #    pass
    #    unique_bond_types = 

    # placeholder; change later
    bonds = 0
    angles = 0
    dihedrals = 0

    with open(filename, 'w') as data:
        data.write(filename+' - created by MoSDeF\n\n')
        data.write('{:d} atoms\n'.format(len(topology.site_list)))
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
                data.write('0 dihedrals\n')

        data.write('{:d} atom types\n'.format(len(set(types))))
        #if atom_style in ['full', 'molecular']:
        #    if bonds != 0:
        #        data.write('{:d} bond types\n'.format(len(set(bond_types))))
        #    if angles != 0:
        #        data.write('{:d} angle types\n'.format(len(set(angle_types))))
        #    if dihedrals != 0:
        #        data.write('{:d} dihedral types\n'.format(len(set(dihedral_types))))

        data.write('\n')
        # Box data
        box.lengths.convert_to_units(u.angstrom)
        box.angles.convert_to_units(u.radian)
        if allclose(box.angles, u.unyt_array([90,90,90],'degree')):
            for i,dim in enumerate(['x', 'y', 'z']):
                data.write('{0:.6f} {1:.6f} {2}lo {2}hi\n'.format(
                    0, box.lengths.value[i],dim))
        else:
            a, b, c = box.lengths.value[i]
            alpha, beta, gamma = box.angles(u.radian)

            lx = a
            xy = b * np.cos(gamma)
            xz = c * np.cos(beta)
            ly = np.sqrt(b**2 - xy**2)
            yz = (b*c*np.cos(alpha) - xy*xz) / ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

            xlo, ylo, zlo = 0 # I believe Topology box is always starting at 0
            xhi = xlo + lx
            yhi = ylo + ly
            zhi = zlo + lz

            xlo_bound = xlo + np.min([0.0, xy, xz, xy+xz])
            xhi_bound = xhi + np.max([0.0, xy, xz, xy+xz])
            ylo_bound = ylo + np.min([0.0, yz])
            yhi_bound = yhi + np.max([0.0, yz])
            zlo_bound = zlo
            zhi_bound = zhi

            data.write('{0} {1} {2}\n'.format(xlo_bound, xhi_bound, xy))
            data.write('{0} {1} {2}\n'.format(ylo_bound, yhi_bound, xz))
            data.write('{0} {1} {2}\n'.format(zlo_bound, zhi_bound, yz))

        # Mass data
        # TODO: Determine where masses are coming from
        #masses = [site.atom_type.mass for site in topology.site_list]
        masses = [0,0] # hardcoded to test right now
        mass_dict = dict([(unique_types.index(atom_type)+1,mass) for atom_type,mass in zip(types,masses)])

        data.write('\nMasses\n\n')
        for atom_type,mass in mass_dict.items():
            data.write('{:d}\t{:.6f}\t# {}\n'.format(atom_type,mass,unique_types[atom_type-1]))
        if forcefield:
            sigmas = [site.atom_type.parameters['sigma'] for site in topology.site_list]
            epsilons = [site.atom_type.parameters['epsilon'] for site in topology.site_list]
            sigma_dict = dict([(unique_types.index(atom_type)+1,sigma) for atom_type,sigma in zip(types,sigmas)])
            epsilon_dict = dict([(unique_types.index(atom_type)+1,epsilon) for atom_type,epsilon in zip(types,epsilons)])

        #    # Modified cross-interactions
            if topology.site_list[0].atom_type.nb_function:
                # Temporary -->
                data.write('\nPair Coeffs # lj\n\n')
                for idx,epsilon in epsilon_dict.items():
                    data.write('{}\t{:.5f}\t{:.5f}\n'.format(idx,epsilon,sigma_dict[idx]))

        #        params = ParameterSet.from_structure(structure)
        #        # Sort keys (maybe they should be sorted in ParmEd)
        #        new_nbfix_types = OrderedDict()
        #        for key, val in params.nbfix_types.items():
        #            sorted_key = tuple(sorted(key))
        #            if sorted_key in new_nbfix_types:
        #                warn('Sorted key matches an existing key')
        #                if new_nbfix_types[sorted_key]:
        #                    warn('nbfixes are not symmetric, overwriting old nbfix')
        #            new_nbfix_types[sorted_key] = params.nbfix_types[key]
        #        params.nbfix_types = new_nbfix_types
        #        warn('Explicitly writing cross interactions using mixing rule: {}'.format(
        #            structure.combining_rule))
        #        coeffs = OrderedDict()
        #        for combo in it.combinations_with_replacement(unique_types, 2):
        #            # Attempt to find pair coeffis in nbfixes
        #            if combo in params.nbfix_types:
        #                type1 = unique_types.index(combo[0])+1
        #                type2 = unique_types.index(combo[1])+1
        #                rmin = params.nbfix_types[combo][0] # Angstrom
        #                epsilon = params.nbfix_types[combo][1] # kcal
        #                sigma = rmin/2**(1/6)
        #                coeffs[(type1, type2)] = (round(sigma, 8), round(epsilon, 8))
        #            else:
        #                type1 = unique_types.index(combo[0]) + 1
        #                type2 = unique_types.index(combo[1]) + 1
        #                # Might not be necessary to be this explicit
        #                if type1 == type2:
        #                    sigma = sigma_dict[type1]
        #                    epsilon = epsilon_dict[type1]
        #                else:
        #                    if structure.combining_rule == 'lorentz':
        #                        sigma = (sigma_dict[type1]+sigma_dict[type2])*0.5
        #                    elif structure.combining_rule == 'geometric':
        #                        sigma = (sigma_dict[type1]*sigma_dict[type2])**0.5
        #                    else:
        #                        raise ValueError('Only lorentz and geometric combining rules are supported')
        #                    epsilon = (epsilon_dict[type1]*epsilon_dict[type2])**0.5
        #                coeffs[(type1, type2)] = (round(sigma, 8), round(epsilon, 8))
        #        if nbfix_in_data_file:
        #            data.write('\nPairIJ Coeffs # modified lj\n\n')
        #            for (type1, type2), (sigma, epsilon) in coeffs.items():
        #                data.write('{0} {1} {2} {3}\n'.format(
        #                    type1, type2, epsilon, sigma))
        #        else:
        #            data.write('\nPair Coeffs # lj\n\n')
        #            print('Copy these commands into your input script:\n')
        #            for (type1, type2), (sigma, epsilon) in coeffs.items():
        #                if type1 == type2:
        #                    data.write('{}\t{:.5f}\t{:.5f}\n'.format(
        #                        type1,epsilon_dict[type1],sigma_dict[type1]))
        #                else:
        #                    print('pair_coeff\t{0} {1} {2} {3}'.format(
        #                        type1, type2, epsilon, sigma))

            # Pair coefficients
            else:
                data.write('\nPair Coeffs # lj\n\n')
                for idx,epsilon in epsilon_dict.items():
                    data.write('{}\t{:.5f}\t{:.5f}\n'.format(idx,epsilon,sigma_dict[idx]))

        #    # Bond coefficients
        #    if bonds:
        #        data.write('\nBond Coeffs # harmonic\n\n')
        #        for params,idx in unique_bond_types.items():
        #            data.write('{}\t{}\t{}\n'.format(idx,*params))

        #    # Angle coefficients
        #    if angles:
        #        data.write('\nAngle Coeffs # harmonic\n\n')
        #        for params,idx in unique_angle_types.items():
        #            data.write('{}\t{}\t{:.5f}\n'.format(idx,*params))

        #    # Dihedral coefficients
        #    if dihedrals:
        #        data.write('\nDihedral Coeffs # opls\n\n')
        #        for params,idx in unique_dihedral_types.items():
        #            opls_coeffs = RB_to_OPLS(params[0],
        #                                     params[1],
        #                                     params[2],
        #                                     params[3],
        #                                     params[4],
        #                                     params[5])
        #            data.write('{}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(idx,*opls_coeffs))


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
                index=i+1,type_index=0,
                zero=0,charge=0,
                x=coords[0],y=coords[1],z=coords[2]))

        #if atom_style in ['full', 'molecular']:
        #    # Bond data
        #    if bonds:
        #        data.write('\nBonds\n\n')
        #        for i,bond in enumerate(bonds):
        #            data.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(
        #                i+1,bond_types[i],bond[0],bond[1]))

        #    # Angle data
        #    if angles:
        #        data.write('\nAngles\n\n')
        #        for i,angle in enumerate(angles):
        #            data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
        #                i+1,angle_types[i],angle[0],angle[1],angle[2]))

        #    # Dihedral data
        #    if dihedrals:
        #        data.write('\nDihedrals\n\n')
        #        for i,dihedral in enumerate(dihedrals):
        #            data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
        #                i+1,dihedral_types[i],dihedral[0],
        #                dihedral[1],dihedral[2],dihedral[3]))
