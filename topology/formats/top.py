import datetime


def write_top(top, filename):
    """Write a topology to a GROMACS .TOP file"""

    _validate_compatibility(top)
    top_vars = _get_top_vars(top)

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
                    atom_type.mass,
                    atom_type.charge,
                    'A',
                    atom_type.parameters['sigma'].in_units(u.nanometer),
                    atom_type.parameters['epsilon'].in_units(u.Unit('kJ')),
                )
            )


        # TODO: for subtop in top.subtops ...
        out_file.write(
            '[ moleculetype ]\n'
            '; name\t\t\tnrexcl
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
                    idx,
                    site.atom_type.name,
                    1, # TODO: subtop idx
                    top.name, # TODO: subtop.name
                    site.element,
                    idx, # TODO: care about charge groups
                    site.charge,
                    site.mass,
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
