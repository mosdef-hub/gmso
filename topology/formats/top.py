import datetime

import unyt as u
from foyer.smarts import SMARTS

from topology.core.element import element_by_mass, element_by_symbol


PARSER = SMARTS()


def write_top(top, filename):
    """Write a topology to a GROMACS .TOP file"""

    _validate_compatibility(top)
    top_vars = _get_top_vars(top)
    _assign_indices(top)

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
            '{0}\t\t'
            '{1}\t\t'
            '{2}\t\t'
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
                '{0}\t\t'
                '{1}\t\t'
                '{2}\t\t'
                '{3}\t\t'
                '{4}\t\t'
                '{5}\t\t'
                '{6}\n'.format(
                    atom_type.name,
                    _lookup_element(atom_type),
                    atom_type.mass.in_units(u.amu).value,
                    atom_type.charge.in_units(u.charge_electron).value,
                    'A',
                    atom_type.parameters['sigma'].in_units(u.nanometer).value,
                    atom_type.parameters['epsilon'].in_units(u.Unit('kJ/mol')).value,
                )
            )


        out_file.write(
            '[ moleculetype ]\n'
            '; name\t\t\tnrexcl\n'
        )

        # TODO: Lookup and join nrexcl from each subtop object
        for subtop_name in set([s.name for s in top.subtops]):
            out_file.write(
                '{0}\t'
                '{1}\n\n'.format(
                    subtop_name,
                    3
                )
            )

        out_file.write(
            '[ atoms ]\n'
            ';   nr       type  resnr residue  atom   cgnr    charge       mass\n'
        )
        for idx, atom_type in enumerate(top.atom_types):
            out_file.write(
                '{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}'
                '\t{5}'
                '\t{6}'
                '\t{7}\n'.format(
                    idx+1,
                    atom_type.name,
                    1, # TODO: subtop idx
                    top.name, # TODO: subtop.name
                    'X', # TODO: establish relationship between atom_type and site ...
                    1, # TODO: care about charge groups
                    atom_type.charge.in_units(u.charge_electron).value,
                    atom_type.mass.in_units(u.amu).value,
                )
            )

        out_file.write(
            '\n[ bonds ]\n'
            ';   nr       type  resnr residue  atom   cgnr    charge       mass\n'
            ';   ai     aj  funct   c0      c1      c2      c3\n'
        )
        for bond_idx, bond in enumerate(top.bonds):
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}\n'.format(
                    bond.connection_members[0].idx,
                    bond.connection_members[1].idx,
                    '1',
                    bond.connection_type.parameters['r_eq'].in_units(u.nm).value,
                    bond.connection_type.parameters['k'].in_units(
                        u.Unit('kJ / (mol*nm**2)')).value,
                )
            )

        out_file.write(
            '\n[ angles ]\n'
            ';   nr       type  resnr residue  atom   cgnr    charge       mass\n'
            ';   ai     aj      ak      funct   c0      c1      c2      c3\n'
        )
        for angle_idx, angle in enumerate(top.angles):
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}'
                '\t{5}\n'.format(
                    angle.connection_members[0].idx,
                    angle.connection_members[1].idx,
                    angle.connection_members[2].idx,
                    '1',
                    angle.connection_type.parameters['theta_eq'].in_units(u.degree).value,
                    angle.connection_type.parameters['k'].in_units(
                        u.Unit('kJ/(mol*rad**2)')).value,
                )
            )

        for dihedral_idx, dihedral in enumerate(top.dihedrals):
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
                    dihedral.connection_members[0].idx,
                    dihderal.connection_members[1].idx,
                    dihedral.connection_members[2].idx,
                    dihedral.connection_members[3].idx,
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
            '[ system ]\n'
            '; name\n'
            '{0}\n\n'.format(
                top.name
            )
        )

        if len(set([s.name for s in top.subtops])) > 1:
            raise NotImplementedError

        out_file.write(
            '[ molecules ]\n'
            '; molecule     nmols\n'
            '{0}\t{1}'.format(top.subtops[0].name, top.n_subtops)
        )


def _validate_compatibility(top):
    """Check compatability of topology object with GROMACS TOP format"""
    pass

def _get_top_vars(top):
    """Generate a dictionary of values for the defaults directive."""
    top_vars = dict()

    top_vars['nbfunc'] = 1
    top_vars['comb_rule'] = 2
    top_vars['gen_pairs'] = 'no'
    top_vars['fudgeLJ'] = 1
    top_vars['fudgeQQ'] = 1

    return top_vars


def _assign_indices(top):
    for idx, site in enumerate(top.sites):
        site.idx = idx + 1


def _lookup_element(atom_type):
    """Attempt to look up an element based on atom type information"""
    elem = None
    while elem is None:
        if atom_type.mass is not None:
            elem = element_by_mass(atom_type.mass)
        if atom_type.name is not None:
            elem = element_by_symbol(atom_type.name)
            elem = element_by_symbol(atom_type.name)
        if atom_type.definition is not None:
            elem = _element_from_smarts_string(atom_type.definition)
        elem = 'X'
    return elem


def _element_from_smarts_string(smarts_string):
    symbol = next(PARSER.parse(smarts_string).find_data('atom_symbol')).children[0]
    return element_by_symbol(symbol)
