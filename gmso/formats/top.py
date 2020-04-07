import datetime

import unyt as u

from gmso.core.element import element_by_atom_type
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.compatibility import check_compatibility
from gmso.exceptions import GMSOError


def write_top(top, filename):
    """Write a gmso.core.Topology object to a GROMACS topology (.TOP) file"""

    _validate_compatibility(top)
    top_vars = _get_top_vars()

    with open(filename, 'w') as out_file:
        out_file.write(
            '; File {} written by topology at {}\n\n'.format(
                top.name if top.name is not None else '',
                str(datetime.datetime.now())
            )
        )
        out_file.write(
            '[ defaults ]\n'
            '; nbfunc\t'
            'comb-rule\t'
            'gen-pairs\t'
            'fudgeLJ\t'
            'fudgeQQ\n'
        )
        out_file.write(
            '{0}\t\t\t'
            '{1}\t\t\t'
            '{2}\t\t\t'
            '{3}\t\t'
            '{4}\n\n'.format(
                top_vars['nbfunc'],
                top_vars['comb_rule'],
                top_vars['gen_pairs'],
                top_vars['fudgeLJ'],
                top_vars['fudgeQQ'],
            )
        )

        out_file.write(
            '[ atomtypes ]\n'
            '; name\t\t'
            'at.num\t\t'
            'mass\t\t'
            'charge\t\t'
            'ptype\t\t'
            'sigma\t\t'
            'epsilon\n'
        )
        for atom_type in top.atom_types:
            out_file.write(
                '{0}\t\t\t'
                '{1}\t\t\t'
                '{2}\t\t'
                '{3}\t\t'
                '{4}\t\t\t'
                '{5}\t\t\t'
                '{6}\n'.format(
                    atom_type.name,
                    _lookup_atomic_number(atom_type),
                    atom_type.mass.in_units(u.amu).value,
                    atom_type.charge.in_units(u.charge_electron).value,
                    'A',
                    atom_type.parameters['sigma'].in_units(u.nanometer).value,
                    atom_type.parameters['epsilon'].in_units(u.Unit('kJ/mol')).value,
                )
            )


        out_file.write(
            '\n[ moleculetype ]\n'
            '; name\t\tnrexcl\n'
        )

        # TODO: Better parsing of subtops into residues/molecules
        n_unique_subtops = len(set([s.name for s in top.subtops]))
        if n_unique_subtops > 1:
            raise NotImplementedError
        # Treat top without subtops as one residue-like "molecule"
        elif n_unique_subtops == 0:
            out_file.write(
                '{0}\t\t\t'
                '{1}\n\n'.format(
                    top.name,
                    3, # Typically exclude 3 nearest neighbors
                )
            )
        # TODO: Lookup and join nrexcl from each subtop object
        elif n_unique_subtops == 1:
            out_file.write(
                '{0}\t\t\t'
                '{1}\n\n'.format(
                    top.name,
                    3
                )
            )

        out_file.write(
            '[ atoms ]\n'
            '; nr\t\ttype\tresnr\tresidue\t\tatom\tcgnr\tcharge\t\tmass\n'
        )
        for idx, site in enumerate(top.sites):
            out_file.write(
                '{0}\t\t\t'
                '{1}\t\t'
                '{2}\t\t'
                '{3}\t'
                '{4}\t\t'
                '{5}\t\t'
                '{6}\t\t'
                '{7}\n'.format(
                    idx+1,
                    site.atom_type.name,
                    1, # TODO: subtop idx
                    top.name, # TODO: subtop.name
                    'X', # TODO: establish relationship between atom_type and site ...
                    1, # TODO: care about charge groups
                    site.charge.in_units(u.charge_electron).value,
                    site.atom_type.mass.in_units(u.amu).value,
                )
            )

        out_file.write(
            '\n[ bonds ]\n'
            ';   ai     aj  funct   c0      c1\n'
        )
        for bond in top.bonds:
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}\n'.format(
                    top.sites.index(bond.connection_members[0]),
                    top.sites.index(bond.connection_members[1]),
                    '1',
                    bond.connection_type.parameters['r_eq'].in_units(u.nm).value,
                    bond.connection_type.parameters['k'].in_units(
                        u.Unit('kJ / (mol*nm**2)')).value,
                )
            )

        out_file.write(
            '\n[ angles ]\n'
            ';   ai     aj      ak      funct   c0      c1\n'
        )
        for angle in top.angles:
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}'
                '\t{5}\n'.format(
                    top.sites.index(angle.connection_members[0]),
                    top.sites.index(angle.connection_members[1]),
                    top.sites.index(angle.connection_members[2]),
                    '1',
                    angle.connection_type.parameters['theta_eq'].in_units(u.degree).value,
                    angle.connection_type.parameters['k'].in_units(
                        u.Unit('kJ/(mol*rad**2)')).value,
                )
            )

        out_file.write(
            '\n[ dihedrals ]\n'
            ';   ai     aj      ak      al  funct   c0      c1      c2\n'
        )
        for dihedral in top.dihedrals:
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}'
                '\t{5}'
                '\t{6}'
                '\t{7}'
                '\t{8}'
                '\t{9}'
                '\t{10}\n'.format(
                    top.sites.index(dihedral.connection_members[0]),
                    top.sites.index(dihedral.connection_members[1]),
                    top.sites.index(dihedral.connection_members[2]),
                    top.sites.index(dihedral.connection_members[3]),
                    '3',
                    dihedral.connection_type.parameters['c0'],
                    dihedral.connection_type.parameters['c1'],
                    dihedral.connection_type.parameters['c2'],
                    dihedral.connection_type.parameters['c3'],
                    0,
                    0,
                )
            )

        out_file.write(
            '\n[ system ]\n'
            '; name\n'
            '{0}\n\n'.format(
                top.name
            )
        )

        if len(set([s.name for s in top.subtops])) > 1:
            raise NotImplementedError

        # TODO: Write out atom types for each unique `subtop` in `atoms` section
        # and write out number of molecules in `molecules` section
        #if len(top.subtops) == 0:
        out_file.write(
            '[ molecules ]\n'
            '; molecule\tnmols\n'
            '{0}\t\t{1}'.format(top.name, 1)
        )
        #elif len(top.subtops) > 0:
        #    out_file.write(
        #        '[ molecules ]\n'
        #        '; molecule\tnmols\n'
        #        '{0}\t\t{1}'.format(top.subtops[0].name, top.n_subtops)
        #    )


def _validate_compatibility(top):
    """Check compatability of topology object with GROMACS TOP format"""
    templates = PotentialTemplateLibrary()
    lennard_jones_potential = templates['LennardJonesPotential']
    harmonic_bond_potential = templates['HarmonicBondPotential']
    harmonic_angle_potential = templates['HarmonicAnglePotential']
    accepted_potentials = [lennard_jones_potential, harmonic_bond_potential, harmonic_angle_potential]
    check_compatibility(top, accepted_potentials)


def _get_top_vars():
    """Generate a dictionary of values for the defaults directive."""
    top_vars = dict()

    top_vars['nbfunc'] = 1
    top_vars['comb_rule'] = 2
    top_vars['gen_pairs'] = 'no'
    top_vars['fudgeLJ'] = 1
    top_vars['fudgeQQ'] = 1

    return top_vars


def _lookup_atomic_number(atom_type):
    """Attempt to look up an atomic_number based on atom type information, 0 if non-element type"""
    try:
        element = element_by_atom_type(atom_type)
        return element.atomic_number
    except GMSOError:
        return 0
