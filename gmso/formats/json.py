"""Serialization to json."""

import json
import warnings
from copy import deepcopy
from pathlib import Path

import unyt as u

from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.box import Box
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.pairpotential_type import PairPotentialType
from gmso.formats.formats_registry import loads_as, saves_as


def _to_json(top, types=True, update=False):
    """Return a json serializable dictionary from a topology.

    This method is used for json serializing the topology

    Parameters
    ----------
    top: gmso.Topology, required
        The topology
    types: bool, default=True
        If true, include type info (i.e. Potentials)
    update: bool, default=False
        If true, update the topology before iterating through the files

    Returns
    -------
    dict
        A json serializable dictionary representing members of this Topology
    """
    if types and not top.is_typed():
        warnings.warn(
            "Cannot incorporate types because the topology is not typed."
        )
        types = False

    if not types and top.is_typed():
        warnings.warn(
            "The provided topology is typed and `types` is set to False. "
            "The types(potentials) info will be lost in the serialized representation. "
            "Please consider using `types=True` if this behavior is not intended. "
        )

    if types and not top.is_fully_typed():
        warnings.warn(
            "The provided topology is not full typed and `types` is set to True. "
            "Please consider using `types=False` if this behavior is not intended. "
        )

    if update:
        top.update_topology()

    json_dict = {
        "name": top._name,
        "scaling_factors": {
            "scaling_factors": top.scaling_factors.tolist(),
            "molecule_scaling_factors": {
                k: v.tolist() for k, v in top.molecule_scaling_factors.items()
            },
        },
        "box": top.box.json_dict() if top.box else None,
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
        "pair_potentialtypes": [],
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
            top.atom_types,
            top.bond_types,
            top.angle_types,
            top.dihedral_types,
            top.improper_types,
        ]:
            for potential in potentials:
                potential_dict = potential.json_dict(
                    exclude={"topology", "set_ref"}
                )
                target = targets[type(potential)]
                potential_dict["id"] = id(potential)
                target.append(potential_dict)

        for pairpotential_type in top._pairpotential_types:
            json_dict["pair_potentialtypes"].append(
                pairpotential_type.json_dict()
            )

    return json_dict


def _set_scaling_factors(top, scaling_factors):
    """Set the global/permolecule scaling factors."""
    global_scaling_factor = scaling_factors["scaling_factors"]
    top.set_scaling_factors(global_scaling_factor[0], global_scaling_factor[1])
    for k, v in scaling_factors["molecule_scaling_factors"].items():
        top.set_scaling_factors(v[0], v[1], molecule_id=k)


def _from_json(json_dict):
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

    # FixMe: DeepCopying a dictionary might not be the most efficient
    # DeepCopying the following structure is done because of the subsequent
    # updates to the dictionary will modify the original one passed as function's
    # argument
    json_dict = deepcopy(json_dict)

    top = Topology(
        name=json_dict["name"],
    )
    _set_scaling_factors(top, json_dict["scaling_factors"])
    id_to_type_map = {}
    for atom_dict in json_dict["atoms"]:
        atom_type_id = atom_dict.pop("atom_type", None)
        atom = Atom.model_validate(atom_dict)
        top.add_site(atom)
        if atom_type_id:
            if not id_to_type_map.get(atom_type_id):
                id_to_type_map[atom_type_id] = []
            id_to_type_map[atom_type_id].append(atom)

    for bond_dict in json_dict["bonds"]:
        bond_type_id = bond_dict.pop("bond_type", None)
        bond_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in bond_dict["connection_members"]
        ]
        bond = Bond.model_validate(bond_dict)
        top.add_connection(bond)
        if bond_type_id:
            if not id_to_type_map.get(bond_type_id):
                id_to_type_map[bond_type_id] = []
            id_to_type_map[bond_type_id].append(bond)

    for angle_dict in json_dict["angles"]:
        angle_type_id = angle_dict.pop("angle_type", None)
        angle_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in angle_dict["connection_members"]
        ]
        angle = Angle.model_validate(angle_dict)
        top.add_connection(angle)
        if angle_type_id:
            if not id_to_type_map.get(angle_type_id):
                id_to_type_map[angle_type_id] = []
            id_to_type_map[angle_type_id].append(angle)

    for dihedral_dict in json_dict["dihedrals"]:
        dihedral_type_id = dihedral_dict.pop("dihedral_type", None)
        dihedral_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in dihedral_dict["connection_members"]
        ]
        dihedral = Dihedral.model_validate(dihedral_dict)
        top.add_connection(dihedral)
        if dihedral_type_id:
            if not id_to_type_map.get(dihedral_type_id):
                id_to_type_map[dihedral_type_id] = []
            id_to_type_map[dihedral_type_id].append(dihedral)

    for improper_dict in json_dict["impropers"]:
        improper_type_id = improper_dict.pop("improper_type", None)
        improper_dict["connection_members"] = [
            top._sites[member_idx]
            for member_idx in improper_dict["connection_members"]
        ]
        improper = Improper.model_validate(improper_dict)
        if improper_type_id:
            if not id_to_type_map.get(improper_type_id):
                id_to_type_map[improper_type_id] = []
            id_to_type_map[improper_type_id].append(improper)

    for atom_type_dict in json_dict["atom_types"]:
        atom_type_id = atom_type_dict.pop("id", None)
        atom_type = AtomType.model_validate(atom_type_dict)
        if atom_type_id in id_to_type_map:
            for associated_atom in id_to_type_map[atom_type_id]:
                associated_atom.atom_type = atom_type

    for connection_types, Creator, attr in [
        (json_dict["bond_types"], BondType, "bond_type"),
        (json_dict["angle_types"], AngleType, "angle_type"),
        (json_dict["dihedral_types"], DihedralType, "dihedral_type"),
        (json_dict["improper_types"], ImproperType, "improper_type"),
    ]:
        for connection_type_dict in connection_types:
            connection_type_id = connection_type_dict.pop("id")
            connection_type = Creator.model_validate(connection_type_dict)
            if connection_type_id in id_to_type_map:
                for associated_connection in id_to_type_map[connection_type_id]:
                    setattr(associated_connection, attr, connection_type)

    if json_dict.get("box") is not None:
        box_dict = json_dict["box"]
        lengths = u.unyt_array(
            box_dict["lengths"]["array"], box_dict["lengths"]["unit"]
        )
        angles = u.unyt_array(
            box_dict["angles"]["array"], box_dict["angles"]["unit"]
        )
        top.box = Box(lengths=lengths, angles=angles)

    top.update_topology()

    # AtomTypes need to be updated for pairpotentialtype addition
    for pair_potentialtype_dict in json_dict["pair_potentialtypes"]:
        pair_potentialtype = PairPotentialType.model_validate(
            pair_potentialtype_dict
        )
        top.add_pairpotentialtype(pair_potentialtype, update=False)

    return top


@saves_as(".json")
def write_json(top, filename, types=True, update=False, **kwargs):
    """Save the topology as a JSON file.

    Parameters
    ----------
    top: gmso.Topology
        The topology to save
    filename: str, pathlib.Path
        The file to save to topology to, must be suffixed with .json
    types: bool, default=True
        If true, include type info (i.e. Potentials)
    update: bool, default=False
        If true, update the topology before iterating through the files
    **kwargs: dict
        The keyword arguments to _to_json and json.dump methods
    """
    json_dict = _to_json(
        top,
        update=update,
        types=types,
    )
    if not isinstance(filename, Path):
        filename = Path(filename).resolve()

    with filename.open("w") as json_file:
        json.dump(json_dict, json_file, **kwargs)


@loads_as(".json")
def load_json(filename):
    """Load a topology from a json file.

    Parameters
    ----------
    filename: str, pathlib.Path
        The file to load the topology from, must be suffixed with .json

    Returns
    -------
    gmso.Topology
        The Topology object obtained by loading the json file
    """
    if not isinstance(filename, Path):
        filename = Path(filename).resolve()

    with filename.open("r") as json_file:
        json_dict = json.load(json_file)
        top = _from_json(json_dict)
        return top
