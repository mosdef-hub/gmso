import warnings 
import unyt as u
import json
import itertools as it

from topology.exceptions import NotYetImplementedWarning, TopologyError
from topology.lib import potential_templates

# Define the supported functions for HOOMD in a dictionary
# Keys are topology.PotentialTemplate
# Values are the hoomd-metadata keyword
supported_nonbonded = {
        potential_templates.LennardJonesPotential(): 'hoomd.md.pair.lj'}

supported_bondtypes = {
        potential_templates.HarmonicBondPotential(): 'hoomd.md.bond.harmonic'}

supported_angletypes = {
        potential_templates.HarmonicAnglePotential(): 'hoomd.md.angle.harmonic'}


def write_metadata(top, filename, unitsystem={'distance': u.nm, 'energy': u.Unit('kJ/mol'),
                                    'mass': u.amu}):
    """ Routine to write hoomdmeta data given a topology

    Parameters
    ----------
    top : topology.Topology
    filename : output file, string
    unitsystem : dictionary
        Keys are fundamental physical dimension
        Values are the associated unit quantity
        Default is a somewhat common unit system for molecular modeling
        For more information, see 
        https://hoomd-blue.readthedocs.io/en/stable/units.html
        """
    # First update the topology in its entirety
    top.update_top()

    metadata = {'objects':[]}
    # Check for supported functions and print to metadata
    _check_supported_bondtypes(top)
    _check_supported_angletypes(top)
    _check_supported_nonbonded(top)

    _write_defaults(top, metadata)
    _write_bondtypes(top, metadata, unitsystem)
    _write_angletypes(top, metadata, unitsystem)
    _write_nonbonded(top, metadata, unitsystem)
    with open(filename, 'w') as f:
        json.dump(metadata, f, indent=4)


def _check_supported_bondtypes(top):
    """ Verify the topology's bond functions are supported in HOOMD """
    for bond_func in top.bond_type_expressions:
        if any([bond_func == func.expression for func in supported_bondtypes]):
            continue
        else:
            raise NotYetImplementedWarning("Bond function <{}> not implemented "
                "or not supported in HOOMD".format(bond_func))

def _check_supported_angletypes(top):
    """ Verify the topology's angle functions are supported in HOOMD """
    for angle_func in top.angle_type_expressions:
        if any([angle_func == func.expression for func in supported_angletypes]):
            continue
        else:
            raise NotYetImplementedWarning("Angle function <{}> not implemented "
                "or not supported in HOOMD".format(angle_func))

def _check_supported_nonbonded(top):
    """ Verify the topology's nonbonded functions are supported in HOOMD """
    for nb_function in top.atom_type_expressions:
        if any([nb_function == func.expression for func in supported_nonbonded]):
            continue
        else:
            raise NotYetImplementedWarning("Nonbonded function <{}> not implemented "
                "or not supported in HOOMD".format(nb_function))

def _write_defaults(top, metadata):
    metadata['objects'].append({'hoomd.md.nlist.cell':{}})

def _write_bondtypes(top, metadata, unitsystem):
    key = "{}-{}"

    for bondtype in top.bond_types:
        possible_hoomd_funcs = [hoomd_func for top_func, hoomd_func in 
                supported_bondtypes.items() 
                if bondtype.expression == top_func.expression]

        if len(possible_hoomd_funcs) > 1:
            warnings.warn("Multiple possible HOOMD functions "
                    + "found for bondtype {}".format(bondtype))

        if len(possible_hoomd_funcs) == 0:
            raise NotYetImplementedWarning(
                    "Bonded function for Bondtype {}".format(bondtype)
                    + " not implemented or not supported in HOOMD")
        else:
            hoomd_func = possible_hoomd_funcs[0]
            if hoomd_func not in metadata:
                metadata[hoomd_func] = {}

            bond_key = key.format(bondtype.member_types[0], 
                bondtype.member_types[1])
            metadata['objects'].append({
                hoomd_func: {
                    "tracked_fields": {
                        'log': True, 
                        'parameters': {
                            bond_key: {
                                'k': _reduce_units(bondtype.parameters['k'], 
                                    unitsystem),
                                'r0': _reduce_units(bondtype.parameters['r_eq'], 
                                    unitsystem)
                                        }
                                        }
                                    }
                        }
                })

def _write_angletypes(top, metadata, unitsystem):
    key = "{}-{}-{}"

    for angletype in top.angle_types:
        possible_hoomd_funcs = [hoomd_func for top_func, hoomd_func in 
                supported_angletypes.items() 
                if angletype.expression == top_func.expression]

        if len(possible_hoomd_funcs) > 1:
            warnings.warn("Multiple possible HOOMD functions "
                    + "found for angletype {}".format(angletype))

        if len(possible_hoomd_funcs) == 0:
            raise NotYetImplementedWarning(
                    "Bonded function for Angletype {}".format(angletype)
                    + " not implemented or not supported in HOOMD")
        else:
            hoomd_func = possible_hoomd_funcs[0]
            if hoomd_func not in metadata:
                metadata[hoomd_func] = {}

            angle_key = key.format(angletype.member_types[0], 
                angletype.member_types[1], angletype.member_types[2])
            
            metadata['objects'].append({
                hoomd_func: {
                    "tracked_fields": {
                        'log': True, 
                        'parameters': {
                            angle_key: {
                                'k': _reduce_units(angletype.parameters['k'], 
                                    unitsystem),
                                't0': _reduce_units(angletype.parameters['theta_eq'], 
                                    unitsystem)
                                        }
                                        }
                                    }
                        }
                })

def _write_nonbonded(top, metadata, unitsystem, r_cut_default=1.2):
    key = "{},{}"

    for first, second in it.combinations_with_replacement(top.atom_types, 2):
        possible_hoomd_nb_funcs = [hoomd_func for top_func, hoomd_func 
                in supported_nonbonded.items() 
                if (first.expression == top_func.expression and
                    second.expression == top_func.expression)]

        if len(possible_hoomd_nb_funcs) > 1:
            warnings.warn("Multiple possible HOOMD functions found for " 
                        + "atom types {0} and {1}: {2}".format(
                            first, second, possible_hoomd_nb_funcs))

        if len(possible_hoomd_nb_funcs) == 0:
            raise NotYetImplementedWarning(
                "Nonbonded function for atomtypes {0} and {1} ".format(first, second)
                + " not implemented or not supported in HOOMD")
        else:
            hoomd_nb_func = possible_hoomd_nb_funcs[0] 
            if hoomd_nb_func not in metadata:
                metadata[hoomd_nb_func] = {}
            if first != second :
                sigma, epsilon = _apply_combining_rule(
                        top.combining_rule, first, second)
            else:
                sigma, epsilon = (first.parameters['sigma'], 
                        first.parameters['epsilon'])
            nb_key = key.format(first.name, second.name)

            metadata['objects'].append({
                hoomd_nb_func: {
                    "arguments": {
                        'r_cut': r_cut_default,
                        'nlist': "Object #0"
                                },
                    "tracked_fields": {
                        'log': True, 
                        'parameters': {
                            nb_key: {
                                'sigma': _reduce_units(sigma, unitsystem),
                                'epsilon': _reduce_units(epsilon, unitsystem),
                                'r_cut': first.parameters.get('r_cut', r_cut_default)
                                    }
                                    }
                                    }
                        }
                })

def _apply_combining_rule(rule, first, second):
    """ Apply combining rule to atomtypes first and second """
    if rule == 'lorentz':
        sigma = (first.parameters['sigma'] + second.parameters['sigma']) / 2
        epsilon = (first.parameters['epsilon'] * second.parameters['epsilon']) ** 0.5
    else:
        warnings.warn("Combining rule {} not supported".format(rule)
                + " defaulting to lorentz")
        sigma = (first.parameters['sigma'] + second.parameters['sigma']) / 2
        epsilon = (first.parameters['epsilon'] * second.parameters['epsilon']) ** 0.5

    return (sigma, epsilon)


def _reduce_units(quantity, unitsystem):
    """ Reduce the units of quantity according to unitsystem """

    # Okay honestly I don't know the best way to handle this
    # But i DO know we want to parse the quantity,
    # figure out its dimensions (is it distance, distance**3/mass, energy/mass, etc)
    # And then divide out by the associated unitsystem quantities
    # There might be some unyt functionality (unyt.unit_registry)? 

    # For now, just strip out the units, and convert to float
    return float(quantity.value)
