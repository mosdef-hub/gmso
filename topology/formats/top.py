import datetime

import unyt as u


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
            '; nbfunc\t\t'
            'comb-rule\t\t'
            'gen-pairs\t\t'
            'fudgeLJ\t\t'
            'fudgeQQ\n'
        )
        out_file.write(
            '{0}\t\t\t\t'
            '{1}\t\t\t\t'
            '{2}\t\t\t\t'
            '{3}\t\t\t\t'
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
                '{0}\t\t\t\t'
                '{1}\t\t\t\t'
                '{2}\t\t\t\t'
                '{3}\t\t\t\t'
                '{4}\t\t\t\t'
                '{5}\t\t\t\t'
                '{6}\n'.format(
                    atom_type.name,
                    1.0, # TODO: Use an atomic number here
                    atom_type.mass.in_units(u.amu).value,
                    atom_type.charge.in_units(u.charge_electron).value,
                    'A',
                    atom_type.parameters['sigma'].in_units(u.nanometer).value,
                    atom_type.parameters['epsilon'].in_units(u.Unit('kJ/mol')).value,
                )
            )


        # TODO: for subtop in top.subtops ...
        out_file.write(
            '[ moleculetype ]\n'
            '; name\t\t\tnrexcl'
            '{0}\t'
            '{1}\n\n'.format(
                top.name, # TODO: subtop.name
                3 # TODO: This should be the nrexcl
            )
        )

        out_file.write(
            '[ atoms ]\n'
            ';   nr       type  resnr residue  atom   cgnr    charge       mass\n'
        )
        for site_idx, site in enumerate(top.sites):
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}'
                '\t{5}'
                '\t{6}'
                '\t{7}'.format(
                    site_idx,
                    site.atom_type.name,
                    1, # TODO: subtop idx
                    top.name, # TODO: subtop.name
                    site.element,
                    site_idx, # TODO: care about charge groups
                    site.charge.in_units(u.charge_electron).value,
                    site.mass.in_units(u.amu).value,
                )
            )

        out_file.write(
            '[ bonds ]\n'
            ';   nr       type  resnr residue  atom   cgnr    charge       mass\n'
            ';   ai     aj  funct   c0      c1      c2      c3\n'
        )
        for bond_idx, bond in enumerate(top.bonds):
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}'.format(
                    bond.connection_members[0].idx,
                    bond.connection_members[1].idx,
                    '1',
                    bond.connection_type.parameters['r_eq'].in_units(u.nm).value,
                    bond.connection_type.parameters['k'].in_units(
                        u.Unit('kJ / (mol*nm**2)')).value,
                )
            )

        out_file.write(
            '[ angles ]\n'
            ';   nr       type  resnr residue  atom   cgnr    charge       mass\n'
            ';   ai     aj      ak      funct   c0      c1      c2      c3\n'
        )
        for angle_idx, angle in enumerate(top.angle):
            out_file.write(
                '\t{0}'
                '\t{1}'
                '\t{2}'
                '\t{3}'
                '\t{4}'
                '\t{5}'.format(
                    angle.connection_members[0].angle_idx,
                    angle.connection_members[1].angle_idx,
                    angle.connection_members[2].angle_idx,
                    '1',
                    angle.connection_type.parameters['theta_eq'].in_units(u.degree).value,
                    angle.connection_type.parameters['k'].in_units(
                        u.Unit('kJ/(mol*rad**2)')).value,
                )
            )

        for dihedral_idx, dihedral in enumerate(top.dihedral):
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
                '\t{10}'.format(
                    dihedral.connection_members[0].dihedral_idx,
                    dihderal.connection_members[1].dihedral_idx,
                    dihedral.connection_members[2].dihedral_idx,
                    dihedral.connection_members[3].dihedral_idx,
                    '3',
                    dihedral.dihedral_type.parameters['c0'],
                    dihedral.dihedral_type.parameters['c1'],
                    dihedral.dihedral_type.parameters['c2'],
                    dihedral.dihedral_type.parameters['c3'],
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

        out_file.write(
            '[ molecules ]\n'
            '; molecule     nmols\n'.format(
                'RES            1'
            )
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
        site.idx = idx
