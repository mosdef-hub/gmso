"""Serialization to json."""
from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType


def to_json(top, types=False, update=False):
    """Return a json serializable from a topology.

    This is used for json serializing the topology

    Parameters
    ----------
    top: gmso.Topology, required
        The topology
    types: bool, default=False
        If true, include type info (i.e. Potentials)
    update: bool, default=False
        If true, update the topology before iterating through the files

    Returns
    -------
    dict
        A json serializable dictionary representing members of this Topology
    """
    if types and not top.is_typed():
        raise ValueError(
            "Cannot incorporate types because the topology is not typed."
        )
    if update:
        top.update_topology()

    json_dict = {
        "name": top._name,
        "scaling_factors": top.scaling_factors,
        "subtopologies": [],
        "atoms": [],
        "bonds": [],
        "angles": [],
        "dihedrals": [],
        "impropers": [],
        "atom_types": [],
        "bond_types": [],
        "angle_types": [],
        "dihedral_types": [],
        "improper_types": [],
    }

    for atom in top._sites:
        atom_dict = atom.json_dict(exclude={"atom_type"})
        if types and atom.atom_type:
            atom_dict["atom_type"] = id(atom.atom_type)

        json_dict["atoms"].append(atom_dict)

    targets = {
        Bond: json_dict["bonds"],
        Angle: json_dict["angles"],
        Dihedral: json_dict["dihedrals"],
        Improper: json_dict["impropers"],
        AtomType: json_dict["atom_types"],
        BondType: json_dict["bond_types"],
        AngleType: json_dict["angle_types"],
        DihedralType: json_dict["dihedral_types"],
        ImproperType: json_dict["improper_types"],
    }

    for connections, exclude_attr in [
        (top._bonds, "bond_type"),
        (top._angles, "angle_type"),
        (top._dihedrals, "dihedral_type"),
        (top._impropers, "improper_type"),
    ]:
        for connection in connections:
            connection_dict = connection.json_dict(
                exclude={exclude_attr, "connection_members"}
            )
            target = targets[type(connection)]
            connection_dict["connection_members"] = [
                top.get_index(member)
                for member in connection.connection_members
            ]
            target.append(connection_dict)
            connection_type = getattr(connection, exclude_attr)
            if types and connection_type:
                connection_dict[exclude_attr] = id(connection_type)
    if types:
        for potentials in [
            top._atom_types.values(),
            top._bond_types.values(),
            top._angle_types.values(),
            top._dihedral_types.values(),
            top._improper_types.values(),
        ]:
            for potential in potentials:
                potential_dict = potential.json_dict(
                    exclude={"topology", "set_ref"}
                )
                target = targets[type(potential)]
                potential_dict["id"] = id(potential)
                target.append(potential_dict)

    for subtop in top.subtops:
        subtop_dict = subtop.json_dict()
        json_dict["subtopologies"].append(subtop_dict)

    return json_dict


def from_json(json_dict):
    """Convert a json_dict into a topology.

    Parameters
    ----------
    json_dict: dict
        The json (dictionary) representation of a Topology

    Returns
    -------
    gmso.Topology
        the equivalent Topology representation from the dictionary
    """
    from gmso.core.topology import Topology

    top = Topology(
        name=json_dict["name"],
    )
    top.scaling_factors = json_dict["scaling_factors"]
    id_to_type_map = {}
    for atom_dict in json_dict["atoms"]:
        atom_type_id = atom_dict.pop("atom_type", None)
        atom = Atom.parse_obj(atom_dict)
        top.add_site(atom)
        if atom_type_id:
            if not id_to_type_map[atom_type_id]:
                id_to_type_map[atom_type_id] = []
            id_to_type_map[atom_type_id].append(atom)

    for bond_dict in json_dict["bonds"]:
        bond_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in bond_dict["connection_members"]
        ]
        bond = Bond.parse_obj(bond_dict)
        top.add_connection(bond)

    for angle_dict in json_dict["angles"]:
        angle_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in angle_dict["connection_members"]
        ]
        angle = Angle.parse_obj(angle_dict)
        top.add_connection(angle)

    for dihedral_dict in json_dict["dihedrals"]:
        dihedral_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in dihedral_dict["connection_members"]
        ]
        dihedral = Dihedral.parse_obj(dihedral_dict)
        top.add_connection(dihedral)

    for improper_dict in json_dict["impropers"]:
        improper_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in improper_dict["connection_members"]
        ]
        improper = Improper.parse_obj(improper_dict)
        top.add_connection(improper)

    for atom_type_dict in json_dict["atom_types"]:
        atom_type_id = atom_type_dict.pop("id", None)
        atom_type = AtomType.parse_obj(atom_type_dict)
        if atom_type_id in id_to_type_map:
            for associated_atom in id_to_type_map[atom_type_id]:
                associated_atom.atom_type = atom_type

    top.update_topology()
    return top
