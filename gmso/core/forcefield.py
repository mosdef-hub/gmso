"""Module for working with GMSO forcefields."""

import copy
import itertools
import warnings
from collections import ChainMap
from pathlib import Path
from typing import Iterable

from lxml import etree
from pydantic import ValidationError

from gmso.core.element import element_by_symbol
from gmso.exceptions import (
    ForceFieldParseError,
    GMSOError,
    MissingPotentialError,
)
from gmso.utils._constants import FF_TOKENS_SEPARATOR
from gmso.utils.decorators import deprecate_function, deprecate_kwargs
from gmso.utils.ff_utils import (
    parse_ff_atomtypes,
    parse_ff_connection_types,
    parse_ff_metadata,
    parse_ff_pairpotential_types,
    validate,
)
from gmso.utils.misc import mask_with, validate_type


def _group_by_expression(potential_types):
    """Group a dictionary of potentials by their expression."""
    expr_group = {}

    for potential in potential_types:
        potential_type = potential_types[potential]
        atom_types_list = expr_group.get(str(potential_type.expression), [])
        atom_types_list.append(potential_type)
        expr_group[str(potential_type.expression)] = atom_types_list

    return expr_group


class ForceField(object):
    """A generic implementation of the forcefield class.

    The ForceField class is one of the core data structures in gmso, which is
    used to hold a collection of gmso.core.Potential subclass objects along with some
    metadata to represent a forcefield. The forcefield object can be applied
    to any gmso.Topology which has effects on its Sites, Bonds, Angles and Dihedrals.

    Parameters
    ----------
    xml_loc : str
        Path to the forcefield xml. The forcefield xml can be either in Foyer or GMSO style.
    strict: bool, default=True
        If true, perform a strict validation of the forcefield XML file
    greedy: bool, default=True
        If True, when using strict mode, fail on the first error/mismatch
    backend: str, default="forcefield-utilities"
        Can be "gmso" or "forcefield-utilities". This will define the methods to
        load the forcefield.

    Attributes
    ----------
    name : str
        Name of the forcefield
    version : str
        Version of the forcefield
    atom_types : dict
        A collection of atom types in the forcefield
    bond_types : dict
        A collection of bond types in the forcefield
    angle_types : dict
        A collection of angle types in the forcefield
    dihedral_types : dict
        A collection of dihedral types in the forcefield
    units : dict
        A collection of unyt.Unit objects used in the forcefield
    scaling_factors : dict
        A collection of scaling factors used in the forcefield

    See Also
    --------
    gmso.ForceField.from_xml
        A class method to create forcefield object from XML files
    gmso.utils.ff_utils.validate
        Function to validate the gmso XML file

    """

    @deprecate_kwargs([("backend", "gmso"), ("backend", "GMSO")])
    def __init__(
        self,
        xml_loc=None,
        strict=True,
        greedy=True,
        backend="forcefield-utilities",
    ):
        if xml_loc is not None:
            if backend in ["gmso", "GMSO"]:
                ff = ForceField.from_xml(xml_loc, strict, greedy)
            elif backend in [
                "forcefield-utilities",
                "forcefield_utilities",
                "ff-utils",
                "ff_utils",
                "ffutils",
            ]:
                ff = ForceField.xml_from_forcefield_utilities(xml_loc)

            else:
                raise (
                    GMSOError(
                        f"Backend provided does not exist. Please provide one of `'gmso'` or \
                `'forcefield-utilities'`"
                    )
                )
            self.name = ff.name
            self.version = ff.version
            self.atom_types = ff.atom_types
            self.bond_types = ff.bond_types
            self.angle_types = ff.angle_types
            self.dihedral_types = ff.dihedral_types
            self.improper_types = ff.improper_types
            self.pairpotential_types = ff.pairpotential_types
            self.potential_groups = ff.potential_groups
            self.scaling_factors = ff.scaling_factors
            self.combining_rule = ff.combining_rule
            self.units = ff.units
        else:
            self.name = "ForceField"
            self.version = "1.0.0"
            self.atom_types = {}
            self.bond_types = {}
            self.angle_types = {}
            self.dihedral_types = {}
            self.improper_types = {}
            self.pairpotential_types = {}
            self.potential_groups = {}
            self.scaling_factors = {}
            self.combining_rule = "geometric"
            self.units = {}

    @property
    def non_element_types(self):
        """Get the non-element types in the ForceField."""
        non_element_types = set()

        for name, atom_type in self.atom_types.items():
            element_symbol = atom_type.get_tag(
                "element"
            )  # FixMe: Should we make this a first class citizen?
            if element_symbol:
                element = element_by_symbol(element_symbol)
                non_element_types.add(element_symbol) if not element else None

        return non_element_types

    @property
    def atom_class_groups(self):
        """Return a dictionary of atomClasses in the Forcefield."""
        atom_types = self.atom_types.values()
        atomclass_dict = {}
        for atom_type in atom_types:
            if atom_type.atomclass is not None:
                atomclass_group = atomclass_dict.get(atom_type.atomclass, [])
                atomclass_group.append(atom_type)
                atomclass_dict[atom_type.atomclass] = atomclass_group
        return atomclass_dict

    def group_atom_types_by_expression(self):
        """Return all AtomTypes in this ForceField with grouped by expression.

        See Also
        --------
        _group_by_expression
            Groups a dictionary of gmso.ParametricPotentials by their expression

        Returns
        -------
        dict
            A dictionary where the key, value -> expression, list of atom_types with that expression
        """
        return _group_by_expression(self.atom_types)

    def group_bond_types_by_expression(self):
        """Return all BondTypes in this ForceField with grouped by expression.

        See Also
        --------
        _group_by_expression
            Groups a dictionary of gmso.ParametricPotentials by their expression

        Returns
        -------
        dict
            A dictionary where the key, value -> expression, list of BondTypes with that expression
        """
        return _group_by_expression(self.bond_types)

    def group_angle_types_by_expression(self):
        """Return all AngleTypes in this ForceField with grouped by expression.

        See Also
        --------
        _group_by_expression
            Groups a dictionary of gmso.ParametricPotentials by their expression

        Returns
        -------
        dict
            A dictionary where the key, value -> expression, list of AngleTypes with that expression
        """
        return _group_by_expression(self.angle_types)

    def group_dihedral_types_by_expression(self):
        """Return all DihedralTypes in this ForceField with grouped by expression.

        See Also
        --------
        _group_by_expression
            Groups a dictionary of gmso.ParametricPotentials by their expression

        Returns
        -------
        dict
            A dictionary where the key, value -> expression, list of DihedralTypes with that expression
        """
        return _group_by_expression(self.dihedral_types)

    def group_improper_types_by_expression(self):
        """Return all ImproperTypes in this ForceField with grouped by expression.

        See Also
        --------
        _group_by_expression
            Groups a dictionary of gmso.ParametricPotentials by their expression

        Returns
        -------
        dict
            A dictionary where the key, value -> expression, list of ImproperTypes with that expression
        """
        return _group_by_expression(self.improper_types)

    def group_pairpotential_types_by_expression(self):
        """Return all PairPotentialTypes in this ForceField with grouped by expression

        See Also
        --------
        _group_by_expression
            Groups a dictionary of gmso.ParametricPotentials by their expression

        Returns
        -------
        dict
            A dictionary where the key, value -> expression, list of PairPotentialTypes with that expression
        """
        return _group_by_expression(self.pairpotential_types)

    def get_potential(self, group, key, return_match_order=False, warn=False):
        """Return a specific potential by key in this ForceField.

        Parameters
        ----------
        group:  {'atom_type', 'bond_type', 'angle_type', 'dihedral_type', 'improper_type'}
            The potential group to perform this search on
        key: str (for atom type) or list of str (for connection types)
            The key to lookup for this potential group
        return_match_order : bool, default=False
            If true, return the order of connection member types/classes that got matched
        warn: bool, default=False
            If true, raise a warning instead of Error if no match found

        Returns
        -------
        gmso.ParametricPotential
            The parametric potential requested

        Raises
        ------
        MissingPotentialError
            When the potential specified by `key` is not found in the ForceField
            potential group `group`
        """
        group = group.lower()

        potential_extractors = {
            "atom_type": self._get_atom_type,
            "bond_type": self._get_bond_type,
            "angle_type": self._get_angle_type,
            "dihedral_type": self._get_dihedral_type,
            "improper_type": self._get_improper_type,
        }

        if group not in potential_extractors:
            raise ValueError(f"Cannot get potential for {group}")

        validate_type(
            (
                [key]
                if isinstance(key, str) or not isinstance(key, Iterable)
                else key
            ),
            str,
        )

        return potential_extractors[group](
            key, return_match_order=return_match_order, warn=warn
        )

    def get_parameters(self, group, key, warn=False, copy=False):
        """Return parameters for a specific potential by key in this ForceField.

        This function uses the `get_potential` function to get Parameters

        See Also
        --------
        gmso.ForceField.get_potential
            Get specific potential/parameters from a forcefield potential group by key
        """
        potential = self.get_potential(group, key, warn=warn)
        return potential.get_parameters(copy=copy)

    def _get_atom_type(self, atom_type, return_match_order=False, warn=False):
        """Get a particular atom_type with given `atom_type` from this ForceField."""
        if isinstance(atom_type, list):
            atom_type = atom_type[0]

        if not self.atom_types.get(atom_type):
            msg = f"AtomType {atom_type} is not present in the ForceField"
            if warn:
                warnings.warn(msg)
            else:
                raise MissingPotentialError(msg)

        return self.atom_types.get(atom_type)

    def _get_bond_type(self, atom_types, return_match_order=False, warn=False):
        """Get a particular bond_type between `atom_types` from this ForceField."""
        if len(atom_types) != 2:
            raise ValueError(
                f"BondType potential can only "
                f"be extracted for two atoms. Provided {len(atom_types)}"
            )

        forward = FF_TOKENS_SEPARATOR.join(atom_types)
        reverse = FF_TOKENS_SEPARATOR.join(reversed(atom_types))
        match = None
        if forward in self.bond_types:
            match = self.bond_types[forward], (0, 1)
        if reverse in self.bond_types:
            match = self.bond_types[reverse], (1, 0)

        msg = (
            f"BondType between atoms {atom_types[0]} and {atom_types[1]} "
            f"is missing from the ForceField"
        )
        if match:
            if return_match_order:
                return match
            else:
                return match[0]
        elif warn:
            warnings.warn(msg)
            return None
        else:
            raise MissingPotentialError(msg)

    def _get_angle_type(self, atom_types, return_match_order=False, warn=False):
        """Get a particular angle_type between `atom_types` from this ForceField."""
        if len(atom_types) != 3:
            raise ValueError(
                f"AngleType potential can only "
                f"be extracted for three atoms. Provided {len(atom_types)}"
            )

        forward = FF_TOKENS_SEPARATOR.join(atom_types)
        reverse = FF_TOKENS_SEPARATOR.join(reversed(atom_types))
        match = None
        if forward in self.angle_types:
            match = self.angle_types[forward], (0, 1, 2)
        if reverse in self.angle_types:
            match = self.angle_types[reverse], (2, 1, 0)

        msg = (
            f"AngleType between atoms {atom_types[0]}, {atom_types[1]} "
            f"and {atom_types[2]} is missing from the ForceField"
        )
        if match:
            if return_match_order:
                return match
            else:
                return match[0]
        elif warn:
            warnings.warn(msg)
            return None
        else:
            raise MissingPotentialError(msg)

    def _get_dihedral_type(
        self, atom_types, return_match_order=False, warn=False
    ):
        """Get a particular dihedral_type between `atom_types` from this ForceField."""
        if len(atom_types) != 4:
            raise ValueError(
                f"DihedralType potential can only "
                f"be extracted for four atoms. Provided {len(atom_types)}"
            )

        forward = FF_TOKENS_SEPARATOR.join(atom_types)
        reverse = FF_TOKENS_SEPARATOR.join(reversed(atom_types))

        match = None
        if forward in self.dihedral_types:
            match = self.dihedral_types[forward], (0, 1, 2, 3)
        if reverse in self.dihedral_types:
            match = self.dihedral_types[reverse], (3, 2, 1, 0)

        if match:
            if return_match_order:
                return match
            else:
                return match[0]

        for i in range(1, 5):
            forward_patterns = mask_with(atom_types, i)
            reverse_patterns = mask_with(reversed(atom_types), i)

            for forward_pattern, reverse_pattern in zip(
                forward_patterns, reverse_patterns
            ):
                forward_match_key = FF_TOKENS_SEPARATOR.join(forward_pattern)
                reverse_match_key = FF_TOKENS_SEPARATOR.join(reverse_pattern)

                if forward_match_key in self.dihedral_types:
                    match = self.dihedral_types[forward_match_key], (0, 1, 2, 3)
                    break

                if reverse_match_key in self.dihedral_types:
                    match = self.dihedral_types[reverse_match_key], (3, 2, 1, 0)
                    break

            if match:
                break

        msg = (
            f"DihedralType between atoms {atom_types[0]}, {atom_types[1]}, "
            f"{atom_types[2]} and {atom_types[3]} is missing from the ForceField."
        )
        if match:
            if return_match_order:
                return match
            else:
                return match[0]
        elif warn:
            warnings.warn(msg)
            return None
        else:
            raise MissingPotentialError(msg)

    def _get_improper_type(
        self, atom_types, return_match_order=False, warn=False
    ):
        """Get a particular improper_type between `atom_types` from this ForceField."""
        if len(atom_types) != 4:
            raise ValueError(
                f"ImproperType potential can only "
                f"be extracted for four atoms. Provided {len(atom_types)}"
            )

        forward = FF_TOKENS_SEPARATOR.join(atom_types)

        if forward in self.improper_types:
            if return_match_order:
                return self.improper_types[forward], (0, 1, 2, 3)
            else:
                return self.improper_types[forward]

        equiv_idx = [
            (0, i, j, k) for (i, j, k) in itertools.permutations((1, 2, 3), 3)
        ]
        equivalent = [
            [atom_types[m], atom_types[n], atom_types[o], atom_types[p]]
            for (m, n, o, p) in equiv_idx
        ]
        for eq, order in zip(equivalent, equiv_idx):
            eq_key = FF_TOKENS_SEPARATOR.join(eq)
            if eq_key in self.improper_types:
                if return_match_order:
                    return self.improper_types[eq_key], order
                else:
                    return self.improper_types[eq_key]

        match = None
        for i in range(1, 5):
            forward_patterns = mask_with(atom_types, i)
            for forward_pattern in forward_patterns:
                forward_match_key = FF_TOKENS_SEPARATOR.join(forward_pattern)
                if forward_match_key in self.improper_types:
                    match = self.improper_types[forward_match_key], (0, 1, 2, 3)
                    break
            if match:
                break
        if not match:
            for i in range(1, 5):
                for eq, order in zip(equivalent, equiv_idx):
                    equiv_patterns = mask_with(eq, i)
                    for equiv_pattern in equiv_patterns:
                        equiv_pattern_key = FF_TOKENS_SEPARATOR.join(
                            equiv_pattern
                        )
                        if equiv_pattern_key in self.improper_types:
                            match = (
                                self.improper_types[equiv_pattern_key],
                                order,
                            )
                            break
                    if match:
                        break
                if match:
                    break

        msg = (
            f"ImproperType between atoms {atom_types[0]}, {atom_types[1]}, "
            f"{atom_types[2]} and {atom_types[3]} is missing from the ForceField."
        )
        if match:
            return match
        elif warn:
            warnings.warn(msg)
            return None
        else:
            raise MissingPotentialError(msg)

    def __repr__(self):
        """Return a formatted representation of the Forcefield."""
        return (
            f"<ForceField {self.name},\n "
            f"{len(self.atom_types)} AtomTypes,\n "
            f"{len(self.bond_types)} BondTypes,\n "
            f"{len(self.angle_types)} AngleTypes,\n "
            f"{len(self.dihedral_types)} DihedralTypes,\n "
            f"{len(self.improper_types)} ImproperType,\n "
            f"{len(self.pairpotential_types)} PairPotentialType,\n "
            f"id: {id(self)}>"
        )

    def __str__(self):
        """Return a string representation of the ForceField."""
        return f"<ForceField {self.name}, id: {id(self)}>"

    def __eq__(self, other):
        # TODO: units don't match between gmso and ffutils loading
        return all(
            [
                self.name == other.name,
                self.version == other.version,
                self.atom_types == other.atom_types,
                self.bond_types == other.bond_types,
                self.angle_types == other.angle_types,
                self.dihedral_types == other.dihedral_types,
                self.improper_types == other.improper_types,
                self.pairpotential_types == other.pairpotential_types,
                self.potential_groups == other.potential_groups,
                self.scaling_factors == other.scaling_factors,
                self.combining_rule == other.combining_rule,
                # self.units == other.units,
            ]
        )

    @classmethod
    def xml_from_forcefield_utilities(cls, filename):
        from forcefield_utilities.xml_loader import FoyerFFs, GMSOFFs

        try:
            loader = GMSOFFs()
            ff = loader.load(filename).to_gmso_ff()
        # Temporarily opt out, pending new forcefield-utilities release
        # except (ForceFieldParseError, FileNotFoundError, ValidationError):
        except:
            loader = FoyerFFs()
            ff = loader.load(filename).to_gmso_ff()
            ff.units = {
                "energy": "kJ",
                "distance": "nm",
                "mass": "amu",
                "charge": "elementary_charge",
            }
        return ff

    def to_xml(self, filename, overwrite=False, backend="gmso"):
        """Get an lxml ElementTree representation of this ForceField

        Parameters
        ----------
        filename: Union[str, pathlib.Path], default=None
            The filename to write the XML file to

        overwrite: bool, default=False
            If True, overwrite an existing file if it exists

        backend: str, default="gmso"
            Can be "gmso" or "forcefield-utilities". This will define the methods to
            write the xml.
        """
        if backend.lower() == "gmso":
            self._xml_from_gmso(filename, overwrite)
        elif backend in [
            "forcefield_utilities",
            "forcefield-utilities",
            "ffutils",
        ]:
            raise NotImplementedError(
                "The forcefield utilities module does not have an xml writer as of yet."
            )
        else:
            raise (
                GMSOError(
                    f"Backend provided does not exist. Please provide one of `'gmso'` or \
            `'forcefield-utilities'`"
                )
            )

    def _xml_from_gmso(self, filename, overwrite=False):
        """Write out an xml file with GMSO as the backend."""
        ff_el = etree.Element(
            "ForceField",
            attrib={"name": str(self.name), "version": str(self.version)},
        )

        metadata = etree.SubElement(ff_el, "FFMetaData")
        if self.scaling_factors.get("electrostatics14Scale") is not None:
            metadata.attrib["electrostatics14Scale"] = str(
                self.scaling_factors.get("electrostatics14Scale")
            )
        if self.scaling_factors.get("nonBonded14Scale") is not None:
            metadata.attrib["nonBonded14Scale"] = str(
                self.scaling_factors.get("nonBonded14Scale")
            )
        if self.combining_rule is not None:
            metadata.attrib["combiningRule"] = str(self.combining_rule)

        # ToDo: ParameterUnitsDefintions and DefaultUnits
        if self.units:
            str_unytDict = copy.copy(self.units)
            for key, value in str_unytDict.items():
                str_unytDict[key] = str(value)
            etree.SubElement(metadata, "Units", attrib=str_unytDict)
        else:
            etree.SubElement(
                metadata,
                "Units",
                attrib={
                    "energy": "kJ",
                    "distance": "nm",
                    "mass": "amu",
                    "charge": "elementary_charge",
                },
            )

        at_groups = self.group_atom_types_by_expression()
        for expr, atom_types in at_groups.items():
            atypes = etree.SubElement(
                ff_el, "AtomTypes", attrib={"expression": expr}
            )
            params_units_def = None
            for atom_type in atom_types:
                if params_units_def is None:
                    params_units_def = {}
                    for param, value in atom_type.parameters.items():
                        params_units_def[param] = value.units
                        etree.SubElement(
                            atypes,
                            "ParametersUnitDef",
                            attrib={
                                "parameter": param,
                                "unit": str(value.units),
                            },
                        )

                atypes.append(atom_type.etree(units=params_units_def))

        bond_types_groups = self.group_bond_types_by_expression()
        angle_types_groups = self.group_angle_types_by_expression()
        dihedral_types_groups = self.group_dihedral_types_by_expression()
        improper_types_groups = self.group_improper_types_by_expression()

        for tag, potential_group in [
            ("BondTypes", bond_types_groups),
            ("AngleTypes", angle_types_groups),
            ("DihedralTypes", dihedral_types_groups),
            ("DihedralTypes", improper_types_groups),
        ]:
            for expr, potentials in potential_group.items():
                potential_group = etree.SubElement(
                    ff_el, tag, attrib={"expression": expr}
                )
                params_units_def = None
                for potential in potentials:
                    if params_units_def is None:
                        params_units_def = {}
                        for param, value in potential.parameters.items():
                            params_units_def[param] = value.units
                            etree.SubElement(
                                potential_group,
                                "ParametersUnitDef",
                                attrib={
                                    "parameter": param,
                                    "unit": str(value.units),
                                },
                            )

                    potential_group.append(potential.etree(params_units_def))

        ff_etree = etree.ElementTree(element=ff_el)

        if not isinstance(filename, Path):
            filename = Path(filename)

        if filename.suffix != ".xml":
            from gmso.exceptions import ForceFieldError

            raise ForceFieldError(
                f"The filename {str(filename)} is not an XML file. "
                f"Please provide filename with .xml extension"
            )

        if not overwrite and filename.exists():
            raise FileExistsError(
                f"File {filename} already exists. Consider "
                f"using overwrite=True if you want to overwrite "
                f"the existing file."
            )

        ff_etree.write(
            str(filename),
            pretty_print=True,
            xml_declaration=True,
            encoding="utf-8",
        )

    @classmethod
    @deprecate_function(
        "The internal `from_xml` will be deprecated soon. Please load the XML with the `xml_from_forcefield_utilities`."
    )
    def from_xml(cls, xmls_or_etrees, strict=True, greedy=True):
        """Create a gmso.Forcefield object from XML File(s).

        This class method creates a ForceField object from the reference
        XML file. This method takes in a single or collection of XML files
        with information about gmso.AtomTypes, gmso.BondTypes, gmso.AngleTypes
        , gmso.PairPotentialTypes and gmso.DihedralTypes to create the ForceField object.

        Parameters
        ----------
        xmls_or_etrees : Union[str, Iterable[str], etree._ElementTree, Iterable[etree._ElementTree]]
            The forcefield XML locations or XML Element Trees
        strict: bool, default=True
            If true, perform a strict validation of the forcefield XML file
        greedy: bool, default=True
            If True, when using strict mode, fail on the first error/mismatch

        Returns
        -------
        forcefield : gmso.ForceField
            A gmso.Forcefield object with a collection of Potential objects
            created using the information in the XML file
        """

        if not isinstance(xmls_or_etrees, Iterable) or isinstance(
            xmls_or_etrees, str
        ):
            xmls_or_etrees = [xmls_or_etrees]

        should_parse_xml = False
        if not (
            all(map(lambda x: isinstance(x, str), xmls_or_etrees))
            or all(
                map(lambda x: isinstance(x, etree._ElementTree), xmls_or_etrees)
            )
        ):
            raise TypeError(
                "Please provide an iterable of strings "
                "as locations of the XML files "
                "or equivalent element Trees"
            )

        if all(map(lambda x: isinstance(x, str), xmls_or_etrees)):
            should_parse_xml = True

        versions = []
        names = []
        ff_atomtypes_list = []
        ff_bondtypes_list = []
        ff_angletypes_list = []
        ff_dihedraltypes_list = []
        ff_pairpotentialtypes_list = []

        atom_types_dict = ChainMap()
        bond_types_dict = {}
        angle_types_dict = {}
        dihedral_types_dict = {}
        improper_types_dict = {}
        pairpotential_types_dict = {}
        potential_groups = {}

        for loc_or_etree in set(xmls_or_etrees):
            validate(loc_or_etree, strict=strict, greedy=greedy)
            ff_tree = loc_or_etree

            if should_parse_xml:
                ff_tree = etree.parse(loc_or_etree)

            ff_el = ff_tree.getroot()
            versions.append(ff_el.attrib["version"])
            names.append(ff_el.attrib["name"])
            ff_meta_tree = ff_tree.find("FFMetaData")

            if ff_meta_tree is not None:
                ff_meta_map = parse_ff_metadata(ff_meta_tree)

            ff_atomtypes_list.extend(ff_tree.findall("AtomTypes"))
            ff_bondtypes_list.extend(ff_tree.findall("BondTypes"))
            ff_angletypes_list.extend(ff_tree.findall("AngleTypes"))
            ff_dihedraltypes_list.extend(ff_tree.findall("DihedralTypes"))
            ff_pairpotentialtypes_list.extend(
                ff_tree.findall("PairPotentialTypes")
            )

        # Consolidate AtomTypes
        for atom_types in ff_atomtypes_list:
            this_atom_types_group = parse_ff_atomtypes(atom_types, ff_meta_map)
            this_atom_group_name = atom_types.attrib.get("name", None)
            if this_atom_group_name:
                potential_groups[this_atom_group_name] = this_atom_types_group
            atom_types_dict.update(this_atom_types_group)

        # Consolidate BondTypes
        for bond_types in ff_bondtypes_list:
            this_bond_types_group = parse_ff_connection_types(
                bond_types, child_tag="BondType"
            )
            this_bond_types_group_name = bond_types.attrib.get("name", None)

            if this_bond_types_group_name:
                potential_groups[this_bond_types_group_name] = (
                    this_bond_types_group
                )

            bond_types_dict.update(this_bond_types_group)

        # Consolidate AngleTypes
        for angle_types in ff_angletypes_list:
            this_angle_types_group = parse_ff_connection_types(
                angle_types, child_tag="AngleType"
            )
            this_angle_types_group_name = angle_types.attrib.get("name", None)

            if this_angle_types_group_name:
                potential_groups[this_angle_types_group_name] = (
                    this_angle_types_group
                )

            angle_types_dict.update(this_angle_types_group)

        # Consolidate DihedralTypes
        for dihedral_types in ff_dihedraltypes_list:
            this_dihedral_types_group = parse_ff_connection_types(
                dihedral_types, child_tag="DihedralType"
            )
            this_improper_types_group = parse_ff_connection_types(
                dihedral_types, child_tag="ImproperType"
            )
            this_group_name = dihedral_types.attrib.get("name", None)

            dihedral_types_dict.update(this_dihedral_types_group)
            improper_types_dict.update(this_improper_types_group)

            if this_group_name:
                this_dihedral_types_group.update(this_improper_types_group)
                potential_groups[this_group_name] = this_dihedral_types_group

        # Consolidate PairPotentialType
        for pairpotential_types in ff_pairpotentialtypes_list:
            this_pairpotential_types_group = parse_ff_pairpotential_types(
                pairpotential_types
            )
            this_pairpotential_types_group_name = (
                pairpotential_types.attrib.get("name", None)
            )

            if this_pairpotential_types_group_name:
                potential_groups[this_pairpotential_types_group_name] = (
                    this_pairpotential_types_group
                )

            pairpotential_types_dict.update(this_pairpotential_types_group)

        ff = cls()
        ff.name = names[0]
        ff.version = versions[0]
        ff.scaling_factors = ff_meta_map["scaling_factors"]
        ff.combining_rule = ff_meta_map["combining_rule"]
        ff.units = ff_meta_map["Units"]
        ff.atom_types = atom_types_dict.maps[0]
        ff.bond_types = bond_types_dict
        ff.angle_types = angle_types_dict
        ff.dihedral_types = dihedral_types_dict
        ff.improper_types = improper_types_dict
        ff.pairpotential_types = pairpotential_types_dict
        ff.potential_groups = potential_groups
        return ff
