import pathlib

from lxml import etree
from gmso.exceptions import ForceFieldParseError


def from_foyer(foyer_xml, gmso_xml=None, overwrite=False):
    """Convert a foyer XML to a gmso XML

    Parameters
    ----------
    foyer_xml : XML file
        An XML file in the foyer format
    gmso_xml: (str, pathlib.Path) default=None
        The output GMSO xml filename. If None, the output XML
        file's name is set by using foyer_xml file suffixed with
        `_gmso`.
    overwrite: bool, default=False
        If true, overwrite the output XML file if it exists

    Raises
    ------
    FileExistsError
        If overwrite is False and the output file (`gmso_xml`) already
        exists
    """
    if gmso_xml is None:
        file_name = pathlib.Path(str(foyer_xml))
        stem = file_name.stem
        gmso_xml = file_name.with_name(stem + "_gmso.xml")
    else:
        gmso_xml = pathlib.Path(gmso_xml)

    if not overwrite and gmso_xml.resolve().exists():
        raise FileExistsError(
            f'The file {gmso_xml.name} already exists. '
            f'Please use a different file name or set '
            f'overwrite=True to overwrite it'
        )

    foyer_xml_tree = etree.parse(foyer_xml)
    f_kwargs = {
        "coulomb14scale": [],
        "lj14scale": [],
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

    _write_gmso_xml(str(gmso_xml), **f_kwargs)


def _write_gmso_xml(gmso_xml, **kwargs):
    """Given the set of keyword arguments, write a gmso topology's forcefield xml file"""
    ff_kwargs = {
        "coulomb14scale": kwargs.get("coulomb14scale"),
        "lj14scale": kwargs.get("lj14scale"),
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
    forcefield = etree.Element("ForceField")
    forcefield.attrib["name"] = pathlib.Path(str(gmso_xml)).stem
    forcefield.attrib["version"] = "0.0.1"

    ffMeta = _create_sub_element(forcefield, "FFMetaData")
    if ff_kwargs["coulomb14scale"]:
        ffMeta.attrib["electrostatics14Scale"] = ff_kwargs["coulomb14scale"]

    if ff_kwargs["lj14scale"]:
        ffMeta.attrib["nonBonded14Scale"] = ff_kwargs["lj14scale"]

    units = _create_sub_element(
        ffMeta,
        name="Units",
        attrib_dict={
            "energy": "kJ/mol",
            "mass": "amu",
            "charge": "coulomb",
            "distance": "nm",
        },
    )

    # AtomTypes and NonBonded Forces
    _write_nbforces(forcefield, ff_kwargs)

    # HarmonicBondTypes
    if len(ff_kwargs["harmonic_bond_types"]) > 0:
        _write_harmonic_bonds(forcefield, ff_kwargs)

    # HarmonicAngleTypes
    if len(ff_kwargs["harmonic_angle_types"]) > 0:
        _write_harmonic_angles(forcefield, ff_kwargs)

    # UreyBradleyAngleTypes
    if len(ff_kwargs["urey_bradley_angle_types"]) > 0:
        _write_ub_angles(forcefield, ff_kwargs)

    # PeriodicTorsionDihedralTypes and PeriodicImproperDihedralTypes
    if len(ff_kwargs["periodic_torsion_dihedral_types"]) > 0:
        _write_periodic_dihedrals(forcefield, ff_kwargs)

    if len(ff_kwargs["periodic_improper_types"]) > 0:
        _write_periodic_impropers(forcefield, ff_kwargs)

    # RBTorsionDihedralTypes
    if len(ff_kwargs["rb_torsion_dihedral_types"]) > 0:
        _write_rb_torsions(forcefield, ff_kwargs)

    ff_tree = etree.ElementTree(forcefield)
    ff_tree.write(
        str(gmso_xml), pretty_print=True, xml_declaration=True, encoding="utf-8"
    )


def _insert_parameters_units_def(root, name, unit):
    params_units_def = _create_sub_element(
        root,
        "ParametersUnitDef",
        attrib_dict={"parameter": name, "unit": unit},
    )
    root.insert(0, params_units_def)


def _add_parameters(root, params_dict):
    parameters = _create_sub_element(root, "Parameters")
    for param_name, param_value in params_dict.items():
        _create_sub_element(
            parameters,
            "Parameter",
            attrib_dict={
                'name': param_name,
                'value': param_value
            }
        )


def _get_dihedral_or_improper_parameters(dihedral_type):
    parameters = {}
    j = 1
    while dihedral_type.get("k{}".format(j)):
        param_k_name = "k{}".format(j)
        param_k_value = dihedral_type.get(param_k_name)
        param_n_name = "n{}".format(j)
        param_n_value = dihedral_type.get("periodicity{}".format(j))
        param_delta_name = "delta{}".format(j)
        param_delta_value = dihedral_type.get("phase{}".format(j))
        parameters[param_k_name] = param_k_value
        parameters[param_n_name] = param_n_value
        parameters[param_delta_name] = param_delta_value
        j += 1
    return parameters, j


def _populate_class_or_type_attrib(root, type_):
    for j, item in enumerate(type_.items()):
        if "type" in item[0]:
            root.attrib["type{}".format(j + 1)] = type_.get(
                "type{}".format(j + 1), "c{}".format(j + 1)
            )
        elif "class" in item[0]:
            root.attrib["type{}".format(j + 1)] = type_.get(
                "class{}".format(j + 1), "c{}".format(j + 1)
            )


def _write_nbforces(forcefield, ff_kwargs):
    # AtomTypes
    nonBondedAtomTypes = _create_sub_element(
        forcefield,
        "AtomTypes",
        attrib_dict={
            "expression": "ep * ((sigma/r)**12 - (sigma/r)**6) + q / (e0 * r)",
        },
    )
    parameters_units = {
        "ep": "kJ/mol",
        "sigma": "nm",
        "e0": "A**2*s**4/(kg*m**3)",
        "q": "coulomb"
    }

    # NonBondedForces
    for name, unit in parameters_units.items():
        _insert_parameters_units_def(
            nonBondedAtomTypes,
            name,
            unit
        )

    for atom_type in ff_kwargs["atom_types"]:
        thisAtomType = _create_sub_element(
            nonBondedAtomTypes,
            "AtomType",
            attrib_dict={
                "name": atom_type.get("name", "AtomType"),
                "atomclass": atom_type.get("class", ""),
                "element": atom_type.get("element", ""),
                "charge": atom_type.get("charge", "0.0"),
                "mass": atom_type.get("mass", "0.0"),
                "definition": atom_type.get("def", ""),
                "description": atom_type.get("desc", ""),
            },
        )

    for i, atom_type in enumerate(ff_kwargs["non_bonded_forces"]):
        thisAtomType = nonBondedAtomTypes.find(
            './/AtomType[@name="{}"]'.format(atom_type.get("type"))
        )
        thisAtomType.attrib["name"] = atom_type.get("type", "AtomType")
        parameters = {
            "ep": atom_type.get("epsilon"),
            "sigma": atom_type.get("sigma"),
            "e0": "8.8542e-12",
            "q": atom_type.get("charge")
        }
        _add_parameters(thisAtomType, parameters)


def _write_harmonic_bonds(forcefield, ff_kwargs):
    harmonicBondTypes = _create_sub_element(
        forcefield, "BondTypes", attrib_dict={"expression": "k * (r-r_eq)**2", }
    )

    parameters_units = {
        "k": "kJ/nm**2",
        "r_eq": "nm"
    }

    for name, unit in parameters_units.items():
        _insert_parameters_units_def(
            harmonicBondTypes,
            name,
            unit
        )

    for i, bond_type in enumerate(ff_kwargs["harmonic_bond_types"]):
        thisBondType = _create_sub_element(
            harmonicBondTypes,
            "BondType",
            attrib_dict={
                "name": bond_type.get("name", "BondType-Harmonic-{}".format(i + 1)),
            },
        )
        _populate_class_or_type_attrib(thisBondType, bond_type)

        parameters = {
            "k": bond_type.get("k", "1.0"),
            "r_eq": bond_type.get("length", "1.0")
        }
        _add_parameters(thisBondType, parameters)


def _write_harmonic_angles(forcefield, ff_kwargs):
    harmonicAngleTypes = _create_sub_element(
        forcefield,
        "AngleTypes",
        attrib_dict={"expression": "k * (theta - theta_eq)**2", },
    )

    parameters_units = {
        "k": "kJ/radian**2",
        "theta_eq": "radian"
    }

    for name, unit in parameters_units.items():
        _insert_parameters_units_def(
            harmonicAngleTypes,
            name,
            unit
        )

    for i, angle_type in enumerate(ff_kwargs["harmonic_angle_types"]):
        thisAngleType = _create_sub_element(
            harmonicAngleTypes,
            "AngleType",
            attrib_dict={
                "name": angle_type.get("name", "AngleType-Harmonic-{}".format(i + 1)),
            },
        )

        parameters = {
            "k": angle_type.get("k", "1.0"),
            "theta_eq": angle_type.get("angle", "1.0")
        }
        _add_parameters(thisAngleType, parameters)

        _populate_class_or_type_attrib(thisAngleType, angle_type)


def _write_ub_angles(forcefield, ff_kwargs):
    ureybradleyAngleTypes = _create_sub_element(
        forcefield, "AngleTypes", attrib_dict={"expression": "k * (w - w_0) ** 2", }
    )

    parameters_units = {
        "k": "kJ/radian**2",
        "w_0": "nm"
    }

    for name, unit in parameters_units.items():
        _insert_parameters_units_def(
            ureybradleyAngleTypes,
            name,
            unit
        )

    for i, angle_type in enumerate(ff_kwargs["urey_bradley_angle_types"]):
        thisAngleType = _create_sub_element(
            ureybradleyAngleTypes,
            "AngleType",
            attrib_dict={
                "name": angle_type.get(
                    "name", "AngleType-UreyBradley-{}".format(i + 1)
                ),
                "type1": angle_type.get("type1", "t1"),
                "type2": angle_type.get("type2", "t2"),
                "type3": angle_type.get("type3", "t3"),
            },
        )

        parameters = {
            "k": angle_type.get("k", "1.0"),
            "w_0": angle_type.get("d", "1.0")
        }
        _add_parameters(thisAngleType, parameters)


def _write_periodic_dihedrals(forcefield, ff_kwargs):
    periodicTorsionDihedralTypes = _create_sub_element(
        forcefield,
        "DihedralTypes",
        attrib_dict={"expression": "k * (1 + cos(n * phi - delta))", },
    )
    max_j = 0
    for i, dihedral_type in enumerate(ff_kwargs["periodic_torsion_dihedral_types"]):
        thisDihedralType = _create_sub_element(
            periodicTorsionDihedralTypes,
            "DihedralType",
            attrib_dict={
                "name": dihedral_type.get(
                    "name", "DihedralType-Periodic-Proper-{}".format(i + 1)
                ),
            },
        )

        _populate_class_or_type_attrib(thisDihedralType, dihedral_type)

        parameters, max_index = _get_dihedral_or_improper_parameters(dihedral_type)
        if max_index > max_j:
            max_j = max_index

        _add_parameters(thisDihedralType, parameters)

    for k in range(0, max_j):
        _insert_parameters_units_def(
            periodicTorsionDihedralTypes,
            'k{}'.format(k),
            'kJ'
        )
        _insert_parameters_units_def(
            periodicTorsionDihedralTypes,
            'n{}'.format(k),
            'dimensionless'
        )
        _insert_parameters_units_def(
            periodicTorsionDihedralTypes,
            'delta{}'.format(k),
            'radian'
        )


def _write_periodic_impropers(forcefield, ff_kwargs):
    max_j = 0
    periodicImproperTypes = _create_sub_element(
        forcefield,
        "DihedralTypes",
        attrib_dict={"expression": "k * (1 + cos(n * phi - delta))", },
    )
    for i, dihedral_type in enumerate(ff_kwargs["periodic_improper_types"]):
        thisDihedralType = _create_sub_element(
            periodicImproperTypes,
            "DihedralType",
            attrib_dict={
                "name": dihedral_type.get(
                    "name", "DihedralType-Periodic-Improper-{}".format(i + 1)
                ),
            },
        )

        _populate_class_or_type_attrib(thisDihedralType, dihedral_type)

        parameters, max_index = _get_dihedral_or_improper_parameters(dihedral_type)
        if max_index > max_j:
            max_j = max_index

        _add_parameters(thisDihedralType, parameters)

    for k in range(0, max_j):
        _insert_parameters_units_def(
            periodicImproperTypes,
            'k{}'.format(k),
            'kJ'
        )
        _insert_parameters_units_def(
            periodicImproperTypes,
            'n{}'.format(k),
            'dimensionless'
        )
        _insert_parameters_units_def(
            periodicImproperTypes,
            'delta{}'.format(k),
            'degree'
        )


def _write_rb_torsions(forcefield, ff_kwargs):
    rbTorsionDihedralTypes = _create_sub_element(
        forcefield,
        "DihedralTypes",
        attrib_dict={
            "expression": "c0 * cos(phi)**0 + c1 * cos(phi)**1 + "
                          "c2 * cos(phi)**2 + c3 * cos(phi)**3 + "
                          "c4 * cos(phi)**4 + c5 * cos(phi)**5",
        },
    )

    max_j = 0
    for i, dihedral_type in enumerate(ff_kwargs["rb_torsion_dihedral_types"]):
        thisDihedralType = _create_sub_element(
            rbTorsionDihedralTypes,
            "DihedralType",
            attrib_dict={
                "name": dihedral_type.get(
                    "name", "DihedralType-RB-Proper-{}".format(i + 1)
                ),
            },
        )

        _populate_class_or_type_attrib(thisDihedralType, dihedral_type)

        parameters = {}

        j = 0
        while dihedral_type.get("c{}".format(j)):
            param_c_name = "c{}".format(j)
            param_c_value = dihedral_type.get(param_c_name)
            parameters[param_c_name] = param_c_value
            j += 1

        _add_parameters(thisDihedralType, parameters)

        if j > max_j:
            max_j = j
    for k in range(0, max_j):
        _insert_parameters_units_def(
            rbTorsionDihedralTypes,
            "c{}".format(k),
            "kJ/mol"
        )


def _create_sub_element(root_el, name, attrib_dict=None):
    sub_el = etree.SubElement(root_el, name, attrib_dict)
    return sub_el
