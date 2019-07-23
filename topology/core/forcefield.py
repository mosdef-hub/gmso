import collections
import copy
import glob
import itertools
import os
from tempfile import mktemp, mkstemp
import xml.etree.ElementTree as ET

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
from pkg_resources import resource_filename
import requests
import warnings
import re

import numpy as np
import parmed as pmd
import foyer.element as custom_elem
import topology 
import unyt as u

from foyer.atomtyper import find_atomtypes
from foyer.exceptions import FoyerError
from foyer import smarts
from foyer.validator import Validator
from foyer.utils.io import import_, has_mbuild

def preprocess_forcefield_files(forcefield_files=None):
    if forcefield_files is None:
        return None

    preprocessed_files = []

    for xml_file in forcefield_files:
        if not hasattr(xml_file,'read'):
            f = open(xml_file)
            _,suffix = os.path.split(xml_file)
        else:
            f = xml_file
            suffix=""

        # read and preprocess
        xml_contents = f.read()
        f.close()
        xml_contents = re.sub(r"(def\w*=\w*[\"\'])(.*)([\"\'])", lambda m: m.group(1) + re.sub(r"&(?!amp;)", r"&amp;", m.group(2)) + m.group(3),
                              xml_contents)

        try:
            '''
            Sort topology objects by precedence, defined by the number of
            `type` attributes specified, where a `type` attribute indicates
            increased specificity as opposed to use of `class`
            '''
            root = ET.fromstring(xml_contents)
            for element in root:
                if 'Force' in element.tag:
                    element[:] = sorted(element, key=lambda child: (
                        -1 * len([attr_name for attr_name in child.keys()
                                    if 'type' in attr_name])))
            xml_contents = ET.tostring(root, method='xml').decode()
        except ET.ParseError:
            '''
            Provide the user with a warning if sorting could not be performed.
            This indicates a bad XML file, which will be passed on to the
            Validator to yield a more descriptive error message.
            '''
            warnings.warn('Invalid XML detected. Could not auto-sort topology '
                          'objects by precedence.')

        # write to temp file
        _, temp_file_name = mkstemp(suffix=suffix)
        with open(temp_file_name, 'w') as temp_f:
            temp_f.write(xml_contents)

        # append temp file name to list
        preprocessed_files.append(temp_file_name)

    return preprocessed_files




def _check_independent_residues(topology):
    """Check to see if residues will constitute independent graphs."""
    for res in topology.residues():
        atoms_in_residue = set([atom for atom in res.atoms()])
        bond_partners_in_residue = [item for sublist in [atom.bond_partners for atom in res.atoms()] for item in sublist]
        # Handle the case of a 'residue' with no neighbors
        if not bond_partners_in_residue:
            continue
        if set(atoms_in_residue) != set(bond_partners_in_residue):
            return False
    return True


def _update_atomtypes(unatomtyped_topology, res_name, prototype):
    """Update atomtypes in residues in a topology using a prototype topology.

    Atomtypes are updated when residues in each topology have matching names.

    Parameters
    ----------
    unatomtyped_topology : openmm.app.Topology
        Topology lacking atomtypes defined by `find_atomtypes`.
    prototype : openmm.app.Topology
        Prototype topology with atomtypes defined by `find_atomtypes`.

    """
    for res in unatomtyped_topology.residues():
        if res.name == res_name:
            for old_atom, new_atom_id in zip([atom for atom in res.atoms()], [atom.id for atom in prototype.atoms()]):
                old_atom.id = new_atom_id

def _error_or_warn(error, msg):
    """Raise an error or warning if topology objects are not fully parameterized.
    
    Parameters
    ----------
    error : bool
        If True, raise an error, else raise a warning
    msg : str
        The message to be provided with the error or warning
    """
    if error:
        raise Exception(msg)
    else:
        warnings.warn(msg)


class Forcefield(object):
    """Forcefield allowing SMARTS based atomtyping.

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.
    name : str, optional, default=None
        Name of a forcefield to load that is packaged within foyer.

    atomtypes: list of Topology.AtomType
    bondtypes : list of Topology.BondType
    angletypes : list of Topology.AngleType


    """
    def __init__(self, forcefield_files=None, name=None, validation=True, debug=False):
        self.atomTypeDefinitions = dict()
        self.atomTypeOverrides = dict()
        self.atomTypeDesc = dict()
        self.atomTypeRefs = dict()
        self._atomTypes = dict()
        self._atomClasses = {'':set()}
        self._included_forcefields = dict()
        self.non_element_types = dict()
        self.atomtypes = list()
        self.bondtypes = list()
        self.angletypes = list()

        all_files_to_load = []
        if forcefield_files is not None:
            if isinstance(forcefield_files, (list, tuple, set)):
                for file in forcefield_files:
                    all_files_to_load.append(file)
            else:
                all_files_to_load.append(forcefield_files)

        if name is not None:
            try:
                file = self.included_forcefields[name]
            except KeyError:
                raise IOError('Forcefield {} cannot be found'.format(name))
            else:
                all_files_to_load.append(file)

        preprocessed_files = preprocess_forcefield_files(all_files_to_load)
        if validation:
            for ff_file_name in preprocessed_files:
                Validator(ff_file_name, debug)
        for xml in preprocessed_files:
            self._generate_potential_terms_from_xml(xml)
        self.parser = smarts.SMARTS(self.non_element_types)

    def _generate_potential_terms_from_xml(self, xml):
        """ From a given xml file, create the topology.Potential objects """
        tree = ET.parse(xml)
        root = tree.getroot()
        self.atomtypes.extend(self._generate_atomtypes(root))
        self.bondtypes.extend(self._generate_bondtypes(root))
        self.angletypes.extend(self._generate_angletypes(root))

    def _generate_atomtypes(self, root):
        """ From an ET, generate topology.AtomType terms"""
        all_atypes =[]
        for atype_ele in root.findall('AtomType'):
            expression = atype_ele.get('expression', default=None)
            ivars = atype_ele.get('independent_variables', default=None)
            mass = atype_ele.get('mass', default="0")
            charge = atype_ele.get('charge', default=None)
            name = atype_ele.get('name', default='')
            aclass = atype_ele.get('class', default='') #This doesn't get used in constructing the AtomType, and seems the same as the name attribute, but we'll parse it just in case
            definition = atype_ele.get('def', default='')
            overrides = atype_ele.get('overrides', default='')
            doi = atype_ele.get('doi', default='')
            desc = atype_ele.get('desc', default='')

            if aclass in self._atomClasses:
                type_set = self._atomClasses[aclass]
            else:
                type_set = set()
                self._atomClasses[aclass] = type_set
            type_set.add(name)
            self._atomClasses[''].add(name)

            if definition:
                self.atomTypeDefinitions[name] = definition
            if overrides:
                overrides = set(atype.strip() for atype
                                in overrides.split(","))
                if overrides:
                    self.atomTypeOverrides[name] = overrides
            if desc:
                self.atomTypeDesc[name] = desc
            if doi:
                dois = set(doi_item.strip() for doi_item in doi.split(','))
                self.atomTypeRefs[name] = dois

            non_parameter_strings = ['expression', 'independent_variables',
                    'mass', 'charge', 'name', 'class', 
                    'def', 'overrides', 'doi', 'desc']
            parameters = {key:_parse_unyt(val) for key,val in atype_ele.items() 
                          if key not in non_parameter_strings}
            new_atype = topology.AtomType(name=name,
                                        expression=expression, 
                                        parameters=parameters,
                                        independent_variables=ivars,
                                        mass=_parse_unyt(mass),
                                        charge=float(charge),
                                        atomclass=aclass)
            all_atypes.append(new_atype)
        return all_atypes


    def _generate_bondtypes(self, root):
        """ From an ET, generate topology.BondType terms"""
        all_btypes =[]
        for btype_ele in root.findall('BondType'):
            class1 = btype_ele.get('class1', default=None)
            class2 = btype_ele.get('class2', default=None)
            type1 = btype_ele.get('type1', default=class1)
            type2 = btype_ele.get('type2', default=class2)

            expression = btype_ele.get('expression', default=None)
            ivars = btype_ele.get('independent_variables', default=None)
            non_parameter_strings = ['expression', 'independent_variables',
                    'class1', 'class2', 'type1', 'type2']
            parameters = {key:_parse_unyt(val) for key,val in btype_ele.items() 
                          if key not in non_parameter_strings}
            
            if len(btype_ele.getchildren()) > 0:
                for subele in btype_ele.getchildren():
                    class1 = subele.get('class1', default=type1)
                    class2 = subele.get('class2', default=type2)
                    type1 = subele.get('type1', default=class1)
                    type2 = subele.get('type2', default=class2)

                    expression = subele.get('expression', default=expression)
                    ivars = subele.get('independent_variables', default=ivars)
                    parameters.update({key:_parse_unyt(val) 
                                        for key,val in subele.items() 
                                        if key not in non_parameter_strings
                                      })

                    new_btype = topology.BondType(expression=expression, 
                                            parameters=copy.copy(parameters), 
                                            independent_variables=ivars,
                                            member_types=[type1, type2])
                    all_btypes.append(new_btype)
            else:
                new_btype = topology.BondType(expression=expression, 
                                        parameters=parameters,
                                        independent_variables=ivars,
                                        member_types=[type1, type2])
                all_btypes.append(new_btype)
        return all_btypes

    def _generate_angletypes(self, root):
        """ From an ET, generate topology.AngleType terms"""
        all_angtypes =[]
        for angtype_ele in root.findall('AngleType'):
            class1 = angtype_ele.get('class1', default=None)
            class2 = angtype_ele.get('class2', default=None)
            class3 = angtype_ele.get('class3', default=None)
            type1 = angtype_ele.get('type1', default=class1)
            type2 = angtype_ele.get('type2', default=class2)
            type3 = angtype_ele.get('type3', default=class3)

            expression = angtype_ele.get('expression', default=None)
            ivars = angtype_ele.get('independent_variables', default=None)
            non_parameter_strings = ['expression', 'independent_variables',
                    'class1', 'class2', 'class3', 'type1', 'type2', 'type3']
            parameters = {key:_parse_unyt(val) for key,val in angtype_ele.items() 
                          if key not in non_parameter_strings}

            if len(angtype_ele.getchildren()) > 0:
                for subele in angtype_ele.getchildren():
                    class1 = subele.get('class1', default=type1)
                    class2 = subele.get('class2', default=type2)
                    class3 = subele.get('class3', default=type3)
                    type1 = subele.get('type1', default=class1)
                    type2 = subele.get('type2', default=class2)
                    type3 = subele.get('type3', default=class3)
                    
                    expression = subele.get('expression', default=expression)
                    ivars = subele.get('independent_variables', default=ivars)
                    parameters.update({key:_parse_unyt(val) 
                                        for key,val in subele.items() 
                                        if key not in non_parameter_strings 
                                      })

                    new_angtype = topology.AngleType(expression=expression, 
                                            parameters=copy.copy(parameters), 
                                            independent_variables=ivars,
                                            member_types=[type1, type2, type3])
                    all_angtypes.append(new_angtype)
            else:
                new_angtype = topology.AngleType(expression=expression, 
                                        parameters=parameters,
                                        independent_variables=ivars,
                                        member_types=[type1, type2, type3])
                all_angtypes.append(new_angtype)
        return all_angtypes

    @property
    def included_forcefields(self):
        if any(self._included_forcefields):
            return self._included_forcefields

        ff_dir = resource_filename('foyer', 'forcefields')
        ff_filepaths = set(glob.glob(os.path.join(ff_dir, '*.xml')))

        for ff_filepath in ff_filepaths:
            _, ff_file = os.path.split(ff_filepath)
            basename, _ = os.path.splitext(ff_file)
            self._included_forcefields[basename] = ff_filepath
        return self._included_forcefields

    def _create_element(self, element, mass):
        if not isinstance(element, elem.Element):
            try:
                element = elem.get_by_symbol(element)
            except KeyError:
                # Enables support for non-atomistic "element types"
                if element not in self.non_element_types:
                    warnings.warn('Non-atomistic element type detected. '
                                  'Creating custom element for {}'.format(element))
                element = custom_elem.Element(number=0,
                                       mass=mass,
                                       name=element,
                                       symbol=element)
            else:
                return element, False

        return element, True

    def apply(self, top, references_file=None, use_residue_map=True,
              assert_bond_params=True, assert_angle_params=True,
              assert_dihedral_params=True, assert_improper_params=False,
              *args, **kwargs):
        """Apply the force field to a molecular structure

        Parameters
        ----------
        top : topoology.Topology or mbuild.Compound
            Molecular structure to apply the force field to
        references_file : str, optional, default=None
            Name of file where force field references will be written (in Bibtex
            format)
        use_residue_map : boolean, optional, default=True
            Store atomtyped topologies of residues to a dictionary that maps
            them to residue names.  Each topology, including atomtypes, will be
            copied to other residues with the same name. This avoids repeatedly
            calling the subgraph isomorphism on idential residues and should
            result in better performance for systems with many identical
            residues, i.e. a box of water. Note that for this to be applied to
            independent molecules, they must each be saved as different
            residues in the topology.
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            proper dihedrals.
        assert_improper_params : bool, optional, default=False
            If True, Foyer will exit if parameters are not found for all system
            improper dihedrals.
        """
        if self.atomTypeDefinitions == {}:
            raise FoyerError('Attempting to atom-type using a force field '
                    'with no atom type defitions.')

        if has_mbuild:
            mb = import_('mbuild')
            if ((not isinstance(top, topology.Topology)) and 
                isinstance(top, mb.Compound)):
                from topology.external.convert_mbuild import from_mbuild
                top = from_mbuild(top)
                #residues = kwargs.get('residues')
                #topology, positions = generate_topology(topology,
                #        self.non_element_types, residues=residues)

        top = self.run_atomtyping(top, use_residue_map=use_residue_map)
        # After we determine atomtypes, we need to apply parameters to ever
        # Connection in our topology. However, this will only
        # look at the Connections the topology is aware of 
        # (just 2-connection Bonds)
        # We would need an extra call to enumerate the angles/dihedrals
        # unless the topology has already enuemrated them
        top = self.parametrize_topology(top)

        '''
        Check that all topology objects (angles, dihedrals, and impropers)
        have parameters assigned. OpenMM will generate an error if bond parameters
        are not assigned.
        '''
        #data = self._SystemData

        #if data.bonds:
        #    missing = [b for b in structure.bonds
        #               if b.type is None]
        #    if missing:
        #        nmissing = len(structure.bonds) - len(missing)
        #        msg = ("Parameters have not been assigned to all bonds. "
        #               "Total system bonds: {}, Parametrized bonds: {}"
        #               "".format(len(structure.bonds), nmissing))
        #        _error_or_warn(assert_bond_params, msg)

        #if data.angles and (len(data.angles) != len(structure.angles)):
        #    msg = ("Parameters have not been assigned to all angles. Total "
        #           "system angles: {}, Parameterized angles: {}"
        #           "".format(len(data.angles), len(structure.angles)))
        #    _error_or_warn(assert_angle_params, msg)

        #proper_dihedrals = [dihedral for dihedral in structure.dihedrals
        #                    if not dihedral.improper]
        #if data.propers and len(data.propers) != \
        #        len(proper_dihedrals) + len(structure.rb_torsions):
        #    msg = ("Parameters have not been assigned to all proper dihedrals. "
        #           "Total system dihedrals: {}, Parameterized dihedrals: {}. "
        #           "Note that if your system contains torsions of Ryckaert-"
        #           "Bellemans functional form, all of these torsions are "
        #           "processed as propers.".format(len(data.propers),
        #           len(proper_dihedrals) + len(structure.rb_torsions)))
        #    _error_or_warn(assert_dihedral_params, msg)

        #improper_dihedrals = [dihedral for dihedral in structure.dihedrals
        #                      if dihedral.improper]
        #if data.impropers and len(data.impropers) != \
        #        len(improper_dihedrals) + len(structure.impropers):
        #    msg = ("Parameters have not been assigned to all impropers. Total "
        #           "system impropers: {}, Parameterized impropers: {}. "
        #           "Note that if your system contains torsions of Ryckaert-"
        #           "Bellemans functional form, all of these torsions are "
        #           "processed as propers".format(len(data.impropers),
        #           len(improper_dihedrals) + len(structure.impropers)))
        #    _error_or_warn(assert_improper_params, msg)

        
        if references_file:
            atom_types = set(site.atom_type for site in top.sites)
            self._write_references_to_file(atom_types, references_file)

        return top

    def parametrize_topology(self, top):
        """ Specify Potential terms for all sites, bonds, angles, etc in top"""
        #for site in top.sites:
            #site.atom_type = _parametrize_atomtype(site)

        for bond in top.bonds:
            bond.connection_type = self._parametrize_bondtype(bond)

        for angle in top.angles:
            angle.connection_type = self._parametrize_angletype(angle)

        #Todo: dihedrals/impropers/etc 

        return top

    def _parametrize_angletype(self, angle):
        """ Specify AngleType term for Angle"""
        found_angletype = [angtype for angtype in self.angletypes 
                if _matching_constituents(angle, angtype)]
        if len(found_angletype) > 1:
            raise FoyerError("Multiple AngleType parameters "
                            "found for {}".format(angle))
        elif len(found_angletype) == 0:
            raise FoyerError("No AngleType parameters "
                            "found for {}".format(bond))
        else:
            return found_angletype[0]

    def _parametrize_bondtype(self, bond):
        """ Specify BondType term for Bond"""
        found_bondtype = [btype for btype in self.bondtypes 
                if _matching_constituents(bond, btype)]
        if len(found_bondtype) > 1:
            raise FoyerError("Multiple BondType parameters "
                            "found for {}".format(bond))
        elif len(found_bondtype) == 0:
            raise FoyerError("No BondType parameters "
                            "found for {}".format(bond))
        else:
            return found_bondtype[0]

    def _parametrize_atomtype(self, site):
        """ Specify AtomType term for Site """
        found_atomtype = [atype for atype in self.atomtypes
                if atype.name == site.id]
        if len(found_atomtype) > 1:
            raise FoyerError("Multiple AtomType parameters "
                            "found for {}".format(site.id))
        else:
            return found_atomtype[0]

    def run_atomtyping(self, top, use_residue_map=True):
        """Atomtype the topology

        Parameters
        ----------
        top : topology.Topology
            Molecular structure to find atom types of
        use_residue_map : boolean, optional, default=True
            Store atomtyped topologies of residues to a dictionary that maps
            them to residue names.  Each topology, including atomtypes, will be
            copied to other residues with the same name. This avoids repeatedly
            calling the subgraph isomorphism on idential residues and should
            result in better performance for systems with many identical
            residues, i.e. a box of water. Note that for this to be applied to
            independent molecules, they must each be saved as different
            residues in the topology.
        """

        find_atomtypes(top, forcefield=self)

        if not all([a.atom_type for a in top.sites]):
            raise ValueError('Not all atoms in topology have atom types')

        return top

    def _write_references_to_file(self, atom_types, references_file):
        atomtype_references = {}
        for atype in atom_types:
            try:
                atomtype_references[atype.name] = self.atomTypeRefs[atype.name]
            except KeyError:
                warnings.warn("Reference not found for atom type '{}'."
                              "".format(atype))
        unique_references = collections.defaultdict(list)
        for atomtype, dois in atomtype_references.items():
            for doi in dois:
                unique_references[doi].append(atomtype)
        unique_references = collections.OrderedDict(sorted(unique_references.items()))
        with open(references_file, 'w') as f:
            for doi, atomtypes in unique_references.items():
                url = "http://dx.doi.org/" + doi
                headers = {"accept": "application/x-bibtex"}
                bibtex_ref = requests.get(url, headers=headers).text
                note = (',\n\tnote = {Parameters for atom types: ' +
                        ', '.join(sorted(atomtypes)) + '}')
                bibtex_ref = bibtex_ref[:-2] + note + bibtex_ref[-2:]
                f.write('{}\n'.format(bibtex_ref))

def _parse_unyt(unyt_string):
    if '*' in unyt_string:
        split_index = unyt_string.find('*')
        val = float(unyt_string[:split_index])
        unit = u.Unit(unyt_string[split_index+1:])
        return val * unit
    else:
        return unyt_string

def _matching_constituents(connection, connectiontype):
    """ Assert that the Connection's members' types are the same as those
    presecribed in the ConnectionType """
    # This is a list of strings based on connection member atomtype names
    c_member_types = [a.atom_type.name for a in connection.connection_members]
    c_member_classes = [a.atom_type.atomclass 
            for a in connection.connection_members]

    # This is a list of strings based on the connectiontype types
    ctype_member_types = [a for a in connectiontype.member_types]

    # Check forward and reverse
    # I think we can just do __eq__ on the list because
    # each item's __eq__ is compared, and each item should just be a string
    if c_member_types == ctype_member_types:
        return True
    if c_member_types[::-1] == ctype_member_types:
        return True
    if c_member_classes == ctype_member_types:
        return True
    if c_member_classes[::-1] == ctype_member_types:
        return True

    return False
