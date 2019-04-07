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
