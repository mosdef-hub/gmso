from warnings import warn

import numpy as np
import unyt as u
import os
import sympy

from lxml import etree
from gmso.utils.io import has_foyer  # Need for foyer
from gmso.exceptions import ForceFieldParseError

if has_foyer:
    import foyer


def from_foyer(foyer_xml, gmso_xml=None):
    """Convert a foyer XML to a gmso XML

    Parameters
    ----------
    foyer_xml : XML
        An XML file in the foyer format 
   
    Returns
    -------
    gmso_xml : XML, default=None
        An XML file in the gmso format
    """

    if gmso_xml is None:
        gmso_xml = os.path.splitext(foyer_xml)[0] + "_gmso.xml"

    foyer_xml_tree = etree.parse(foyer_xml)
    f_kwargs = {
        "atom_types": [],
        "non_bonded_forces": [],
        "harmonic_bond_types": [],
        "harmonic_angle_types": [],
        "urey_bradley_angle_types": [],
        "rb_torsion_dihedral_types": [],
        "periodic_torsion_dihedral_types": [],
        "periodic_improper_types": [],
    }

    # Try to load in AtomType section
    # Load in AtomTypes section otherwise
    atom_types_el = foyer_xml_tree.findall("AtomType")
    if len(atom_types_el) == 0:
        atom_types_el = foyer_xml_tree.findall("AtomTypes")
        if len(atom_types_el) == 0:
            raise ForceFieldParseError

    for atom_types in atom_types_el:
        for atom_type in atom_types.getiterator("Type"):
            f_kwargs["atom_types"].append(atom_type)
    nonbonded_force_el = foyer_xml_tree.findall("NonbondedForce")
    for atom_types in nonbonded_force_el:
        f_kwargs["coulomb14scale"] = atom_types.attrib["coulomb14scale"]
        f_kwargs["lj14scale"] = atom_types.attrib["lj14scale"]
        for atom_type in atom_types.getiterator("Atom"):
            f_kwargs["non_bonded_forces"].append(atom_type)

    harmonic_bond_force_el = foyer_xml_tree.findall("HarmonicBondForce")
    for hbf in harmonic_bond_force_el:
        for bond_type in hbf.getiterator("Bond"):
            f_kwargs["harmonic_bond_types"].append(bond_type)

    harmonic_angle_force_el = foyer_xml_tree.findall("HarmonicAngleForce")
    for haf in harmonic_angle_force_el:
        for angle_type in haf.getiterator("Angle"):
            f_kwargs["harmonic_angle_types"].append(angle_type)

    urey_bradley_angle_el = foyer_xml_tree.findall("AmoebaUreyBradleyForce")
    for ubf in urey_bradley_angle_el:
        for angle_type in ubf.getiterator("UreyBradley"):
            f_kwargs["urey_bradley_angle_types"].append(angle_type)

    periodic_torsion_force_el = foyer_xml_tree.findall("PeriodicTorsionForce")
    for ptf in periodic_torsion_force_el:
        for dihedral_type in ptf.getiterator("Proper"):
            f_kwargs["periodic_torsion_dihedral_types"].append(dihedral_type)
        for dihedral_type in ptf.getiterator("Improper"):
            f_kwargs["periodic_improper_types"].append(dihedral_type)

    rb_torsion_force_el = foyer_xml_tree.findall("RBTorsionForce")
    for rbf in rb_torsion_force_el:
        for dihedral_type in rbf.getiterator("Proper"):
            f_kwargs["rb_torsion_dihedral_types"].append(dihedral_type)

    _write_gmso_xml(gmso_xml, **f_kwargs)


def _write_gmso_xml(gmso_xml, **kwargs):
    """Given the set of keyword arguments, write a gmso topology's forcefield xml file"""
    ff_kwargs = {
        "atom_types": kwargs.get("atom_types", []),
        "non_bonded_forces": kwargs.get("non_bonded_forces", []),
        "harmonic_bond_types": kwargs.get("harmonic_bond_types", []),
        "harmonic_angle_types": kwargs.get("harmonic_angle_types", []),
        "urey_bradley_angle_types": kwargs.get("urey_bradley_angle_types", []),
        "periodic_torsion_dihedral_types": kwargs.get(
            "periodic_torsion_dihedral_types", []
        ),
        "periodic_improper_types": kwargs.get("periodic_improper_types", []),
        "rb_torsion_dihedral_types": kwargs.get("rb_torsion_dihedral_types", []),
    }
    forceField = etree.Element("ForceField")
    forceField.attrib["name"] = os.path.splitext(gmso_xml)[0]
    forceField.attrib["version"] = "0.0.1"

    ffMeta = _create_subelement(forceField, "FFMetaData")
    units = _create_subelement(
            ffMeta,
            name="Units",
            attrib_dict={
                "energy": "kJ/mol",
                "mass": "amu",
                "charge": "coulomb",
                "distance": "nm",
                }
            )
    # AtomTypes
    nonBondedAtomTypes = _create_subelement(forceField,
            "AtomTypes",
            attrib_dict={
                "expression": "ep * ((sigma/r)**12 - (sigma/r)**6) + q / (e0 * r)",
                }
            )


    # NonBondedForces
    nonBondedAtomTypesParamsUnitsDef_ep = _create_subelement(
            nonBondedAtomTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "ep",
                "unit": "kJ/mol",
                }
            )
    nonBondedAtomTypesParamsUnitsDef_sigma = _create_subelement(
            nonBondedAtomTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "sigma",
                "unit": "nm",
                }
            )
    nonBondedAtomTypesParamsUnitsDef_e0 = _create_subelement(
            nonBondedAtomTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "e0",
                "unit": "A**2*s**4/(kg*m**3)",
                }
            )
    nonBondedAtomTypesParamsUnitsDef_q = _create_subelement(
            nonBondedAtomTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "q",
                "unit": "coulomb",
                }
            )

    # HarmonicBondTypes
    harmonicBondTypes = _create_subelement(
            forceField,
            "BondTypes",
            attrib_dict={
                "expression": "k * (r-r_eq)**2",
                }
            )

    harmonicBondTypesParamsUnitsDef_k = _create_subelement(
            harmonicBondTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "k",
                "unit": "kJ/nm**2",
                }
            )

    harmonicBondTypesParamsUnitsDef_r_eq = _create_subelement(
            harmonicBondTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "r_eq",
                "unit": "nm"
                }
            )

    # HarmonicAngleTypes
    harmonicAngleTypes = _create_subelement(
            forceField,
            "AngleTypes",
            attrib_dict={
                "expression": "k * (theta - theta_eq)**2",
                }
            )

    harmonicAngleTypesParamsUnitsDef_k = _create_subelement(
            harmonicAngleTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "k",
                "unit": "kJ/radian**2",
                }
    )

    harmonicBondTypesParamsUnitsDef_theta_eq = _create_subelement(
            harmonicAngleTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "theta_eq",
                "unit": "radian",
                }
            )

    # UreyBradleyAngleTypes
    ureybradleyAngleTypes = _create_subelement(
            forceField,
            "AngleTypes",
            attrib_dict={
                "expression": "k * (w - w_0) ** 2",
                }
            )

    ureybradleyAngleTypesParamsUnitsDef_k = _create_subelement(
            ureybradleyAngleTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "k",
                "unit": "kJ/radian**2"
                }
            )

    ureybradleyAngleTypesParamsUnitsDef_w_0 = _create_subelement(
            ureybradleyAngleTypes,
            "ParametersUnitDef",
            attrib_dict={
                "parameter": "w_0",
                "unit": "nm"
                }
            )

    # PeriodicTorsionDihedralTypes
    periodicTorsionDihedralTypes = _create_subelement(
            forceField,
            "DihedralTypes",
            attrib_dict={
                "expression": "k * (1 + cos(n * phi - delta))",
                }
            )

    # PeriodicImproperDihedralTypes
    periodicImproperTypes = _create_subelement(
            forceField,
            "DihedralTypes",
            attrib_dict={
                "expression": "k * (1 + cos(n * phi - delta))",
                }
            )

    # RBTorsionDihedralTypes
    rbTorsionDihedralTypes = _create_subelement(
            forceField,
            "DihedralTypes",
            attrib_dict={
                "expression": "c0 * cos(phi)**0 + c1 * cos(phi)**1 + c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + c5 * cos(phi)**5",
                }
            )

    for atom_type in ff_kwargs["atom_types"]:
        thisAtomType = etree.SubElement(nonBondedAtomTypes, "AtomType")
        thisAtomType.attrib["name"] = atom_type.get("name", "AtomType")
        thisAtomType.attrib["atomclass"] = atom_type.get("class", "")
        thisAtomType.attrib["element"] = atom_type.get("element", "")
        thisAtomType.attrib["charge"] = atom_type.get("charge", "0.0")
        thisAtomType.attrib["mass"] = atom_type.get("mass", "0.0")
        thisAtomType.attrib["definition"] = atom_type.get("def", "")
        thisAtomType.attrib["description"] = atom_type.get("desc", "")

    for i, atom_type in enumerate(ff_kwargs["non_bonded_forces"]):
        thisAtomType = nonBondedAtomTypes.find(
            './/AtomType[@name="{}"]'.format(atom_type.get("type"))
        )
        thisAtomType.attrib["name"] = atom_type.get("type", "AtomType")
        parameters = etree.SubElement(thisAtomType, "Parameters")
        parameter_ep = etree.SubElement(parameters, "Parameter")
        parameter_ep.attrib["name"] = "ep"
        parameter_ep.attrib["value"] = atom_type.get("epsilon")
        parameter_sigma = etree.SubElement(parameters, "Parameter")
        parameter_sigma.attrib["name"] = "sigma"
        parameter_sigma.attrib["value"] = atom_type.get("sigma")

        parameter_e0 = etree.SubElement(parameters, "Parameter")
        parameter_e0.attrib["name"] = "e0"
        parameter_e0.attrib["value"] = "8.8542e-12"
        parameter_q = etree.SubElement(parameters, "Parameter")
        parameter_q.attrib["name"] = "q"
        parameter_q.attrib["value"] = atom_type.get("charge")

    for i, bond_type in enumerate(ff_kwargs["harmonic_bond_types"]):
        thisBondType = etree.SubElement(harmonicBondTypes, "BondType")
        thisBondType.attrib["name"] = bond_type.get(
            "name", "BondType-Harmonic-{}".format(i + 1)
        )
        for j, item in enumerate(bond_type.items()):
            if "type" in item[0]:
                thisBondType.attrib["type{}".format(j + 1)] = bond_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisBondType.attrib["type{}".format(j + 1)] = bond_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisBondType, "Parameters")
        parameter_k = etree.SubElement(parameters, "Parameter")
        parameter_k.attrib["name"] = "k"
        parameter_k.attrib["value"] = bond_type.get("k", "1.0")
        parameter_r_eq = etree.SubElement(parameters, "Parameter")
        parameter_r_eq.attrib["name"] = "r_eq"
        parameter_r_eq.attrib["value"] = bond_type.get("length", "1.0")

    for i, angle_type in enumerate(ff_kwargs["harmonic_angle_types"]):
        thisAngleType = etree.SubElement(harmonicAngleTypes, "AngleType")
        thisAngleType.attrib["name"] = angle_type.get(
            "name", "AngleType-Harmonic-{}".format(i + 1)
        )
        for j, item in enumerate(angle_type.items()):
            if "type" in item[0]:
                thisAngleType.attrib["type{}".format(j + 1)] = angle_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisAngleType.attrib["type{}".format(j + 1)] = angle_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisAngleType, "Parameters")
        parameter_k = etree.SubElement(parameters, "Parameter")
        parameter_k.attrib["name"] = "k"
        parameter_k.attrib["value"] = angle_type.get("k", "1.0")

        parameter_r_eq = etree.SubElement(parameters, "Parameter")
        parameter_r_eq.attrib["name"] = "theta_eq"
        parameter_r_eq.attrib["value"] = angle_type.get("angle", "1.0")

    for i, angle_type in enumerate(ff_kwargs["urey_bradley_angle_types"]):
        thisAngleType = etree.SubElement(ureybradleyAngleTypes, "AngleType")
        thisAngleType.attrib["name"] = angle_type.get(
            "name", "AngleType-UreyBradley-{}".format(i + 1)
        )
        thisAngleType.attrib["type1"] = angle_type.get("type1", "t1")
        thisAngleType.attrib["type2"] = angle_type.get("type2", "t2")
        thisAngleType.attrib["type3"] = angle_type.get("type3", "t2")

        parameters = etree.SubElement(thisAngleType, "Parameters")
        parameter_k = etree.SubElement(parameters, "Parameter")
        parameter_k.attrib["name"] = "k"
        parameter_k.attrib["value"] = angle_type.get("k", "1.0")

        parameter_r_eq = etree.SubElement(parameters, "Parameter")
        parameter_r_eq.attrib["name"] = "w_0"
        parameter_r_eq.attrib["value"] = angle_type.get("d", "1.0")

    max_j = 0

    for i, dihedral_type in enumerate(ff_kwargs["periodic_torsion_dihedral_types"]):
        thisDihedralType = etree.SubElement(
            periodicTorsionDihedralTypes, "DihedralType"
        )
        thisDihedralType.attrib["name"] = dihedral_type.get(
            "name", "DihedralType-Periodic-Proper-{}".format(i + 1)
        )

        for j, item in enumerate(dihedral_type.items()):
            if "type" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisDihedralType, "Parameters")

        j = 1
        while dihedral_type.get("k{}".format(j)):
            parameter_k_name = "k{}".format(j)
            parameter_k = etree.SubElement(parameters, "Parameter")
            parameter_k.attrib["name"] = parameter_k_name
            parameter_k.attrib["value"] = dihedral_type.get(parameter_k_name)
            parameter_n = etree.SubElement(parameters, "Parameter")
            parameter_n.attrib["name"] = "n{}".format(j)
            parameter_n.attrib["value"] = dihedral_type.get("periodicity{}".format(j))

            parameter_delta = etree.SubElement(parameters, "Parameter")
            parameter_delta.attrib["name"] = "delta{}".format(j)
            parameter_delta.attrib["value"] = dihedral_type.get("phase{}".format(j))
            j += 1
        if j > max_j:
            max_j = j
    for k in range(0, max_j):
        periodicTorsionDihedralTypesParamsUnitsDef_k = etree.Element(
            "ParametersUnitDef"
        )
        periodicTorsionDihedralTypesParamsUnitsDef_k.attrib["parameter"] = "k{}".format(
            k
        )
        periodicTorsionDihedralTypesParamsUnitsDef_k.attrib["unit"] = "kJ"
        periodicTorsionDihedralTypes.insert(
            0, periodicTorsionDihedralTypesParamsUnitsDef_k
        )

        periodicTorsionDihedralTypesParamsUnitsDef_n = etree.Element(
            "ParametersUnitDef"
        )
        periodicTorsionDihedralTypesParamsUnitsDef_n.attrib["parameter"] = "n{}".format(
            k
        )
        periodicTorsionDihedralTypesParamsUnitsDef_n.attrib["unit"] = "dimensionless"
        periodicTorsionDihedralTypes.insert(
            0, periodicTorsionDihedralTypesParamsUnitsDef_n
        )

        periodicTorsionDihedralTypesParamsUnitsDef_del = etree.Element(
            "ParametersUnitDef"
        )
        periodicTorsionDihedralTypesParamsUnitsDef_del.attrib[
            "parameter"
        ] = "delta{}".format(k)
        periodicTorsionDihedralTypesParamsUnitsDef_del.attrib["unit"] = "radian"
        periodicTorsionDihedralTypes.insert(
            0, periodicTorsionDihedralTypesParamsUnitsDef_del
        )

    for i, dihedral_type in enumerate(ff_kwargs["rb_torsion_dihedral_types"]):
        thisDihedralType = etree.SubElement(rbTorsionDihedralTypes, "DihedralType")
        thisDihedralType.attrib["name"] = dihedral_type.get(
            "name", "DihedralType-RB-Proper-{}".format(i + 1)
        )

        for j, item in enumerate(dihedral_type.items()):
            if "type" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisDihedralType, "Parameters")

        j = 0
        while dihedral_type.get("c{}".format(j)):
            parameter_c_name = "c{}".format(j)
            parameter_c = etree.SubElement(parameters, "Parameter")
            parameter_c.attrib["name"] = parameter_c_name
            parameter_c.attrib["value"] = dihedral_type.get(parameter_c_name)
            j += 1
        if j > max_j:
            max_j = j
    for k in range(0, max_j):
        rbTorsionDihedralTypesParamsUnitsDef_c = etree.Element("ParametersUnitDef")
        rbTorsionDihedralTypesParamsUnitsDef_c.attrib["parameter"] = "c{}".format(k)
        rbTorsionDihedralTypesParamsUnitsDef_c.attrib["unit"] = "kJ/mol"
        rbTorsionDihedralTypes.insert(0, rbTorsionDihedralTypesParamsUnitsDef_c)

    for i, dihedral_type in enumerate(ff_kwargs["periodic_improper_types"]):
        thisDihedralType = etree.SubElement(periodicImproperTypes, "DihedralType")
        thisDihedralType.attrib["name"] = dihedral_type.get(
            "name", "DihedralType-Periodic-Improper-{}".format(i + 1)
        )

        for j, item in enumerate(dihedral_type.items()):
            if "type" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisDihedralType, "Parameters")

        j = 1
        while dihedral_type.get("k{}".format(j)):
            parameter_k_name = "k{}".format(j)
            parameter_k = etree.SubElement(parameters, "Parameter")
            parameter_k.attrib["name"] = parameter_k_name
            parameter_k.attrib["value"] = dihedral_type.get(parameter_k_name)
            parameter_n = etree.SubElement(parameters, "Parameter")
            parameter_n.attrib["name"] = "n{}".format(j)
            parameter_n.attrib["value"] = dihedral_type.get("periodicity{}".format(j))

            parameter_delta = etree.SubElement(parameters, "Parameter")
            parameter_delta.attrib["name"] = "delta{}".format(j)
            parameter_delta.attrib["value"] = dihedral_type.get("phase{}".format(j))
            j += 1
        if j > max_j:
            max_j = j
    for k in range(0, max_j):
        periodicImproperTypesParamsUnitsDef_k = etree.Element("ParametersUnitDef")
        periodicImproperTypesParamsUnitsDef_k.attrib["parameter"] = "k{}".format(k)
        periodicImproperTypesParamsUnitsDef_k.attrib["unit"] = "kJ"
        periodicImproperTypes.insert(0, periodicImproperTypesParamsUnitsDef_k)

        periodicImproperTypesParamsUnitsDef_n = etree.Element("ParametersUnitDef")
        periodicImproperTypesParamsUnitsDef_n.attrib["parameter"] = "n{}".format(k)
        periodicImproperTypesParamsUnitsDef_n.attrib["unit"] = "dimensionless"
        periodicImproperTypes.insert(0, periodicImproperTypesParamsUnitsDef_n)

        periodicImproperTypesParamsUnitsDef_del = etree.Element("ParametersUnitDef")
        periodicImproperTypesParamsUnitsDef_del.attrib["parameter"] = "delta{}".format(
            k
        )
        periodicImproperTypesParamsUnitsDef_del.attrib["unit"] = "degree"
        periodicImproperTypes.insert(0, periodicImproperTypesParamsUnitsDef_del)

    ff_tree = etree.ElementTree(forceField)
    ff_tree.write(gmso_xml, pretty_print=True, xml_declaration=True, encoding="utf-8")

def _create_subelement(root_el, name, attrib_dict=None):
    sub_el = etree.SubElement(root_el, name, attrib_dict)
    return sub_el
