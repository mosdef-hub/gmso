"""GMSO and foyer use specific residues to apply force fields and mapping molecule number to atom numbers."""

import os
from warnings import warn
from xml.dom import minidom

import mbuild as mb
from forcefield_utilities.xml_loader import FoyerFFs, GMSOFFs
from mbuild.compound import Compound
from mbuild.utils.io import has_foyer

import gmso
from gmso.core.views import PotentialFilters
from gmso.external.convert_mbuild import from_mbuild as mb_convert
from gmso.parameterization import apply as gmso_apply


def specific_ff_to_residue(
    structure,
    forcefield_selection=None,
    gmso_match_ff_by="molecule",
    residues=None,
    boxes_for_simulation=1,
):
    """
    Take the mbuild Compound or mbuild Box and applies the selected FF to the corresponding residue via foyer and GMSO.

    Note: a residue is defined as a molecule in this case, so it is not
    designed for applying a force field to a protein.

    Parameters
    ----------
    structure: mbuild Compound object or mbuild Box object;
        The mBuild Compound object with box lengths and angles or mbuild Box object, which contains the molecules
        (or empty box) that will have the force field applied to them.
    forcefield_selection: str or dictionary, default=None
        Apply a force field to the output file by selecting a force field xml file with
        its path or by using the standard force field name provided the `foyer` package.
        Example dict for FF file: {'ETH': 'oplsaa.xml', 'OCT': 'path_to file/trappe-ua.xml'}
        Example str for FF file: 'path_to_file/trappe-ua.xml'
        Example dict for standard FF names: {'ETH': 'oplsaa', 'OCT': 'trappe-ua'}
        Example str for standard FF names: 'trappe-ua'
        Example of a mixed dict with both: {'ETH': 'oplsaa', 'OCT': 'path_to_file/'trappe-ua.xml'}
    gmso_match_ff_by: str ("group" or "molecule"), default = "molecule"
        How the GMSO force field is applied, using the molecules name/residue name (mbuild.Compound.name)
        for GOMC and NAMD.  This is regardless number of levels in the mbuild.Compound.

        * "molecule" applies the force field using the molecule's name or atom's name for a
        single atom molecule (1 atom/bead  = molecule).
            -  Molecule > 1 atom ----> uses the "mbuild.Compound.name" (1 level above the atoms/beads)
            as the molecule's name.  This "mb.Compound.name" (1 level above the atoms/beads) needs
            to be used in the Charmm object's residue_list and forcefield_selection (if >1 force field), and
            will be the residue name in the PSF, PDB, and FF files.
            -  Molecule = 1 atom/bead  ----> uses the "atom/bead's name" as the molecule's name.
            This "atom/bead's name" needs to be used in the Charmm object's residue_list and
            forcefield_selection (if >1 force field), and will be the residue name in the PSF, PDB, and FF files.

            NOTE: Non-bonded zeolites or other fixed structures without bonds will use the
            "atom/bead's name" as the molecule's name, if they are single non-bonded atoms.
            However, the user may want to use the "group" option instead for this type of system,
            if applicable (see the "group" option).

            - Example (Charmm_writer selects the user changeable "ETH", in the Charmm object residue_list and
            forcefield_selection (if >1 force field), which sets the residue "ETH" in the PSF, PDB, and FF files):
                ethane = mbuild.load("CC", smiles=True)
                ethane.name = "ETH"

                ethane_box = mbuild.fill_box(
                compound=[ethane],
                n_compounds=[100],
                box=[4, 4, 4]
                )

            - Example (Charmm_writer must to select the non-user changeable "_CH4" (per the foyer TraPPE force field),
            in the Charmm object residue_list and forcefield_selection (if >1 force field),
            which sets the residue "_CH4" in the PSF, PDB, and FF files):

            methane_ua_bead_name = "_CH4"
            methane_child_bead = mbuild.Compound(name=methane_ua_bead_name)
            methane_box = mbuild.fill_box(
                compound=methane_child_bead, n_compounds=4, box=[1, 2, 3]
            )
            methane_box.name = "MET"

            - Example (Charmm_writer must to select the non-user changeable "Na" and "Cl"
            in the Charmm object residue_list and forcefield_selection (if >1 force field),
            which sets the residues "Na" and "Cl" in the PSF, PDB, and FF files):

            sodium_atom_name = "Na"
            sodium_child_atom = mbuild.Compound(name=sodium_atom_name)
            sodium = mb.Compound(name="SOD")
            sodium.add(sodium_child_atom, inherit_periodicity=False)

            chloride_atom_name = "Cl"
            chloride_child_bead = mbuild.Compound(name=chloride_atom_name)
            chloride = mb.Compound(name="CHL")
            chloride.add(chloride_child_atom, inherit_periodicity=False)

            sodium_chloride_box = mbuild.fill_box(
                compound=[sodium, chloride],
                n_compounds=[4, 4],
                box=[1, 2, 3]
            )

            - Example zeolite (Charmm_writer must to select the non-user changeable "Si" and "O"
            in the Charmm object residue_list and forcefield_selection (if >1 force field),
            which sets the residues "Si" and "O" in the PSF, PDB, and FF files):

            lattice_cif_ETV_triclinic = load_cif(file_or_path=get_mosdef_gomc_fn("ETV_triclinic.cif"))
            ETV_triclinic = lattice_cif_ETV_triclinic.populate(x=1, y=1, z=1)
            ETV_triclinic.name = "ETV"

        * "group" applies the force field 1 level under the top level mbuild.Compound, if only 2 levels exist,
        taking the top levels name mbuild.Compound. Or in other words:
            - For only 2 level (mbuild container-particles) group will grab the name of the mbuild container
            - For > 2 levels (e.g., mbuild container-etc-molecule-residue-particle),
            group will grab 1 level down from top

        This is ideal to use when you are building simulation box(es) using mbuild.fill_box(),
        with molecules, and it allows you to add another level to single atom molecules
        (1 atom/bead  = molecule) to rename the mbuild.Compound().name, changing the residue's
        name and allowing keeping the atom/bead name so the force field is applied properly.

        WARNING: This "group" option will take all the molecule below it, regardless if they
        are selected to be separte from the group via the residue_list and forcefield_selection
        (if >1 force field).

        NOTE: This "group" option may be best for non-bonded zeolites or other fixed structures
        without bonds, if they are single non-bonded atoms. Using this "group" option, the user
        can select the residue name for the Charmm_writer's residue_list and forcefield_selection
        (if >1 force field) to force field all the atoms with a single residue name, and output
        this residue name in the PSF, PDB, and FF files.

            - Example (Charmm_writer select the user changeable "MET" in the Charmm object residue_list
            and forcefield_selection (if >1 force field), which sets the residue "MET" in the
            PSF, PDB, and FF files):

            methane_ua_bead_name = "_CH4"
            methane_child_bead = mbuild.Compound(name=methane_ua_bead_name)
            methane_box = mbuild.fill_box(
                compound=methane_child_bead, n_compounds=4, box=[1, 2, 3]
            )
            methane_box.name = "MET"

            - Example (Charmm_writer select the user changeable "MET" in the Charmm object residue_list
            and forcefield_selection (if >1 force field), which sets the residue "MET" in the
            PSF, PDB, and FF files):

            methane_ua_bead_name = "_CH4"
            methane_molecule_name = "MET"
            methane = mb.Compound(name=methane_molecule_name)
            methane_child_bead = mb.Compound(name=methane_ua_bead_name)
            methane.add(methane_child_bead, inherit_periodicity=False)

            methane_box = mb.fill_box(
                compound=methane, n_compounds=10, box=[1, 2, 3]
            )

            - Example (Charmm_writer select the user changeable "MET" in the Charmm object residue_list
            and forcefield_selection (if >1 force field), which sets the residue "MET" in the
            PSF, PDB, and FF files):

            methane_child_bead = mb.Compound(name="_CH4")
            methane = mb.Compound(name="MET")
            methane.add(methane_child_bead, inherit_periodicity=False)

            box_liq = mb.fill_box(
                compound=methane,
                n_compounds=1230,
                box=[4.5, 4.5, 4.5]
            )

            - Example zeolite (Charmm_writer select the user changeable "ETV" in the Charmm object residue_list
            and forcefield_selection (if >1 force field), which sets the residue
            "ETV" in the PSF, PDB, and FF files):

            lattice_cif_ETV_triclinic = load_cif(file_or_path=get_mosdef_gomc_fn("ETV_triclinic.cif"))
            ETV_triclinic = lattice_cif_ETV_triclinic.populate(x=1, y=1, z=1)
            ETV_triclinic.name = "ETV"

    residues: list, [str, ..., str], default=None
        Labels of unique residues in the Compound. Residues are assigned by
        checking against Compound.name.  Only supply residue names as 4 characters
        strings, as the residue names are truncated to 4 characters to fit in the
        psf and pdb file.
    boxes_for_simulation: either 1 or 2, default=1
        Gibbs (GEMC) or grand canonical (GCMC) ensembles are examples of where the boxes_for_simulation would be 2.
        Canonical (NVT) or isothermal–isobaric (NPT) ensembles are example with the boxes_for_simulation equal to 1.
        Note: the only valid options are 1 or 2.

    Returns
    -------
    list, [
        topology,
        unique_topology_groups_lists,
        residues_applied_list,
        electrostatics14Scale_dict,
        nonBonded14Scale_dict,
        atom_types_dict,
        bond_types_dict,
        angle_types_dict,
        dihedral_types_dict,
        improper_types_dict,
        combining_rule,
        ]

    topology: gmso.Topology
        gmso Topology with applied force field
    unique_topology_groups_list: list
        list of residues (i.e., list of strings).
        These are all the residues in which the force field actually applied.
    electrostatics14Scale_dict: dict
        A dictionary with the 1,4-electrostatic/Coulombic scalars for each residue,
        as the forcefields are specified by residue {'residue_name': '1-4_electrostatic_scaler'}.
    nonBonded14Scale_dict: dict
        A dictionary with the 1,4-non-bonded scalars for each residue,
        as the forcefields are specified by residue {'residue_name': '1-4_nonBonded_scaler'}.
    atom_types_dict: dict
        A dict with the all the residues as the keys. The unique values are a list containing,
        {'expression': confirmed singular atom types expression or equation,
        'atom_types': gmso Topology.atom_types}.
    bond_types_dict: dict
        A dict with the all the residues as the keys. The unique  values are a list containing,
        {'expression': confirmed singular bond types expression or equation,
        'bond_types': gmso Topology.bond_types}.
    angle_types_dict: dict
        A dict with the all the residues as the keys. The unique  values are a list containing,
        {'expression': confirmed singular angle types expression or equation,
        'angle_types': gmso Topology.angle_types}.
    dihedral_types_dict: dict
        A dict with the all the residues as the keys. The unique  values are a list containing,
        {'expression': confirmed singular dihedral types expression or equation,
        'dihedral_types': gmso Topology.dihedral_types}.
    improper_types_dict: dict
        A dict with the all the residues as the keys. The unique  values are a list containing,
        {'expression': confirmed singular improper types expression or equation,
        'improper_types': gmso Topology.improper_types}.
    combining_rule: str
        The possible mixing/combining  rules are 'geometric' or 'lorentz',
        which provide the  geometric and arithmetic mixing rule, respectively.
        NOTE: Arithmetic means the 'lorentz' combining or mixing rule.
        NOTE: GMSO default to the 'lorentz' mixing rule if none is provided,
        and this writers default is the GMSO default.

    Notes
    -----
    To write the NAMD/GOMC force field, pdb, psf, and force field
    (.inp) files, the residues and forcefields must be provided in
    a str or dictionary. If a dictionary is provided all residues must
    be specified to a force field if the boxes_for_simulation is equal to 1.

    Generating an empty box (i.e., pdb and psf files):
    Enter residues = [], but the accompanying structure must be an empty mb.Box.
    However, when doing this, the forcefield_selection must be supplied,
    or it will provide an error (i.e., forcefield_selection can not be equal to None).

    In this current FF/psf/pdb writer, a residue type is essentially a molecule type.
    Therefore, it can only correctly write systems where every bead/atom in the molecule
    has the same residue name, and the residue name is specific to that molecule type.
    For example: a protein molecule with many residue names is not currently supported,
    but is planned to be supported in the future.
    """
    if has_foyer:
        from foyer import Forcefield
        from foyer.forcefields import forcefields
    else:
        error_msg = (
            "Package foyer is not installed. "
            "Please install it using conda install -c conda-forge foyer"
        )
        raise ImportError(error_msg)

    # Validate inputs
    _validate_boxes_for_simulation(boxes_for_simulation)
    _validate_forcefield_selection_and_residues(forcefield_selection, residues)
    new_gmso_topology, initial_no_atoms = _validate_structure(
        structure, residues
    )
    gmso_compatible_forcefield_selection = _validate_forcefields(
        forcefield_selection, residues
    )

    # can use  match_ff_by="group" or "molecule", group was only chosen so everything is using the
    # user selected mb.Compound.name...
    gmso_apply(
        new_gmso_topology,
        gmso_compatible_forcefield_selection,
        speedup_by_molgraph=True,
        identify_connections=True,
        match_ff_by=gmso_match_ff_by,
        speedup_by_moltag=True,
        remove_untyped=True,
    )
    new_gmso_topology.update_topology()

    # find mixing rule.  If an empty.box mixing rule is set to None
    combining_rule = (
        new_gmso_topology._combining_rule
        if isinstance(structure, Compound)
        else None
    )

    # identify the bonded atoms and hence the molecule, label the GMSO objects
    # and create the function outputs.
    molecule_number = 0  # 0 sets the 1st molecule_number at 1
    molecules_atom_number_dict = {}
    unique_topology_groups_list = []
    unique_topologies_groups_dict = {}
    atom_types_dict = {}
    bond_types_dict = {}
    angle_types_dict = {}
    dihedral_types_dict = {}
    improper_types_dict = {}
    nonBonded14Scale_dict = {}
    electrostatics14Scale_dict = {}

    for unique_group in new_gmso_topology.unique_site_labels(
        gmso_match_ff_by, name_only=True
    ):
        if unique_group is not None:
            unique_topology_groups_list.append(unique_group)

    for unique_group in unique_topology_groups_list:
        unique_subtop_group = new_gmso_topology.create_subtop(
            label_type=gmso_match_ff_by, label=unique_group
        )
        unique_topologies_groups_dict[unique_group] = unique_subtop_group

        nb_scalers_list = new_gmso_topology.get_lj_scale(
            molecule_id=unique_group
        )
        electro_scalers_list = new_gmso_topology.get_electrostatics_scale(
            molecule_id=unique_group
        )
        nonBonded14Scale_dict[unique_group] = (
            None if nb_scalers_list is None else nb_scalers_list[2]
        )
        electrostatics14Scale_dict[unique_group] = (
            None if electro_scalers_list is None else electro_scalers_list[2]
        )

    _cross_check_residues_and_unique_site_labels(
        structure, residues, unique_topology_groups_list, boxes_for_simulation
    )

    # get all the bonded atoms, which is used for the bonded map to identify molecules
    bonded_atom_number_set = set()
    all_bonded_atoms_list = set()
    for bond in new_gmso_topology.bonds:
        bonded_atom_0_iter = new_gmso_topology.get_index(
            bond.connection_members[0]
        )
        bonded_atom_1_iter = new_gmso_topology.get_index(
            bond.connection_members[1]
        )
        bonded_atom_tuple_iter = sorted(
            [bonded_atom_0_iter, bonded_atom_1_iter]
        )
        bonded_atom_number_set.add(tuple(bonded_atom_tuple_iter))
        all_bonded_atoms_list.update(bonded_atom_tuple_iter)

    # TODO: Refactor this section, might be able to use unique_site_by_labels(name_only=False)?
    # map all bonded atoms as molecules
    molecules_atom_number_list = []
    for site_idx, site in enumerate(new_gmso_topology.sites):
        if site_idx in all_bonded_atoms_list:
            for bonded_atoms_n in bonded_atom_number_set:
                if site_idx in bonded_atoms_n:
                    if len(molecules_atom_number_list) != 0:
                        for initiated_molecule in molecules_atom_number_list:
                            if site_idx in initiated_molecule:
                                initiated_molecule.update(bonded_atoms_n)
                                break
                            elif (
                                initiated_molecule
                                == molecules_atom_number_list[-1]
                            ):
                                molecules_atom_number_list.append(
                                    {bonded_atoms_n[0], bonded_atoms_n[1]}
                                )
                                break
                    else:
                        molecules_atom_number_list.append(
                            {bonded_atoms_n[0], bonded_atoms_n[1]}
                        )
        else:
            molecules_atom_number_list.append({site_idx})

    # create a molecule number to atom number dict
    # Example:  {molecule_number_x: {atom_number_1, ..., atom_number_y}, ...}
    for molecule in molecules_atom_number_list:
        molecules_atom_number_dict.update({molecule_number: molecule})
        molecule_number += 1

    for site in new_gmso_topology.sites:
        site_atom_number_iter = new_gmso_topology.get_index(site)
        # get molecule number
        for mol_n, atom_set_n in molecules_atom_number_dict.items():
            if site_atom_number_iter in atom_set_n:
                molecule_p_number = mol_n

        if gmso_match_ff_by == "group":
            site.__dict__["residue_name_"] = site.__dict__["group_"]
        elif gmso_match_ff_by == "molecule":
            site.__dict__["residue_name_"] = site.__dict__["molecule_"].name

        site.__dict__["residue_number_"] = molecule_p_number + 1

    # create a topolgy only with the bonded parameters, including their residue/molecule type
    # which permit force fielding in GOMC easier in the charmm_writer
    # iterate thru the unique topologies and get the unique atom, bond, angle, dihedral, improper types
    for (
        unique_top_group_name_iter,
        unique_top_iter,
    ) in unique_topologies_groups_dict.items():
        # get the unique non-bonded data, equations, and other info
        atom_type_expression_set = set()
        for atom_type in unique_top_iter.atom_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        ):
            atom_type.__dict__["tags_"] = {
                "resname": unique_top_group_name_iter
            }
            atom_type_expression_set.add(atom_type.expression)
        if len(atom_type_expression_set) == 1:
            atom_types_dict.update(
                {
                    unique_top_group_name_iter: {
                        "expression": list(atom_type_expression_set)[0],
                        "atom_types": unique_top_iter.atom_types(
                            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
                        ),
                    }
                }
            )
        elif len(atom_type_expression_set) == 0:
            atom_types_dict.update({unique_top_group_name_iter: None})
        else:
            raise ValueError(
                "There is more than 1 nonbonded equation types per residue or molecules "
            )

        # get the unique bond data, equations, and other info
        bond_type_expression_set = set()
        for bond_type in unique_top_iter.bond_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        ):
            bond_type.__dict__["tags_"] = {
                "resname": unique_top_group_name_iter
            }
            bond_type_expression_set.add(bond_type.expression)
        if len(bond_type_expression_set) == 1:
            bond_types_dict.update(
                {
                    unique_top_group_name_iter: {
                        "expression": list(bond_type_expression_set)[0],
                        "bond_types": unique_top_iter.bond_types(
                            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
                        ),
                    }
                }
            )
        elif len(bond_type_expression_set) == 0:
            bond_types_dict.update({unique_top_group_name_iter: None})
        else:
            raise ValueError(
                "There is more than 1 bond equation types per residue or molecules "
            )

        # get the unique angle data, equations, and other info
        angle_type_expression_set = set()
        for angle_type in unique_top_iter.angle_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        ):
            angle_type.__dict__["tags_"] = {
                "resname": unique_top_group_name_iter
            }
            angle_type_expression_set.add(angle_type.expression)
        if len(angle_type_expression_set) == 1:
            angle_types_dict.update(
                {
                    unique_top_group_name_iter: {
                        "expression": list(angle_type_expression_set)[0],
                        "angle_types": unique_top_iter.angle_types(
                            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
                        ),
                    }
                }
            )
        elif len(angle_type_expression_set) == 0:
            angle_types_dict.update({unique_top_group_name_iter: None})
        else:
            raise ValueError(
                "There is more than 1 angle equation types per residue or molecules "
            )

        # get the unique dihedral data, equations, and other info
        dihedral_type_expression_set = set()
        for dihedral_type in unique_top_iter.dihedral_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        ):
            dihedral_type.__dict__["tags_"] = {
                "resname": unique_top_group_name_iter
            }
            dihedral_type_expression_set.add(dihedral_type.expression)
        if len(dihedral_type_expression_set) == 1:
            dihedral_types_dict.update(
                {
                    unique_top_group_name_iter: {
                        "expression": list(dihedral_type_expression_set)[0],
                        "dihedral_types": unique_top_iter.dihedral_types(
                            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
                        ),
                    }
                }
            )
        elif len(dihedral_type_expression_set) == 0:
            dihedral_types_dict.update({unique_top_group_name_iter: None})
        else:
            raise ValueError(
                "There is more than 1 dihedral equation types per residue or molecules "
            )

        # get the unique improper data, equations, and other info
        improper_type_expression_set = set()
        for improper_type in unique_top_iter.improper_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
        ):
            improper_type.__dict__["tags_"] = {
                "resname": unique_top_group_name_iter
            }
            improper_type_expression_set.add(improper_type.expression)
        if len(improper_type_expression_set) == 1:
            improper_types_dict.update(
                {
                    unique_top_group_name_iter: {
                        "expression": list(improper_type_expression_set)[0],
                        "improper_types": unique_top_iter.improper_types(
                            filter_by=PotentialFilters.UNIQUE_NAME_CLASS
                        ),
                    }
                }
            )
        elif len(improper_type_expression_set) == 0:
            improper_types_dict.update({unique_top_group_name_iter: None})
        else:
            raise ValueError(
                "There is more than 1 improper equation types per residue or molecules "
            )

        # check to see if the non-bonded and electrostatic 1-4 interactions are in each group/molecule/residue
        if unique_top_group_name_iter not in list(
            nonBonded14Scale_dict.keys()
        ) or unique_top_group_name_iter not in list(
            electrostatics14Scale_dict.keys()
        ):
            raise ValueError(
                f"The {unique_top_group_name_iter} residue is not provided for the "
                f'{"nonBonded14Scale"} and {"electrostatics14Scale"} values'
            )

    topology = new_gmso_topology
    # calculate the final number of atoms
    final_no_atoms = topology.n_sites

    if final_no_atoms != initial_no_atoms:
        error_msg = (
            "The initial number of atoms sent to the force field analysis is "
            "not the same as the final number of atoms analyzed. "
            f"The initial number of atoms was {initial_no_atoms}"
            f"and the final number of atoms was {final_no_atoms}. "
            "Please ensure that all the residues names that are in the initial "
            "Compound are listed in the residues list "
            "(i.e., the residues variable)."
        )
        raise ValueError(error_msg)

    return [
        topology,
        unique_topology_groups_list,
        electrostatics14Scale_dict,
        nonBonded14Scale_dict,
        atom_types_dict,
        bond_types_dict,
        angle_types_dict,
        dihedral_types_dict,
        improper_types_dict,
        combining_rule,
    ]


def _validate_structure(structure, residues):
    """Validate if input is an mb.Compound with initialized box or mb.Box."""
    if isinstance(structure, (Compound, mb.Box)):
        error_msg = f"The structure, {mb.Compound} or {mb.Box}, needs to have have box lengths and angles."
        if isinstance(structure, Compound):
            if structure.box is None:
                raise TypeError(error_msg)

        elif isinstance(structure, mb.Box):
            if structure.lengths is None or structure.angles is None:
                raise TypeError(error_msg)
    else:
        error_msg = (
            "The structure expected to be of type: "
            f"{mb.Compound} or {mb.Box}, received: {type(structure)}"
        )
        raise TypeError(error_msg)

    # Check to see if it is an empty mbuild.Compound and set intial atoms to 0
    # note empty mbuild.Compound will read 1 atoms but there is really noting there
    # calculate the initial number of atoms for later comparison
    # flatten the mbuild.compound, which is needed to get the mbuild.compound.names correctly from
    # unflattend mbuild.compound.
    # if mbuild.box do not flatten as it must be an empty box.
    if isinstance(structure, Compound):
        if len(structure.children) == 0:
            # there are no real atoms in the Compound so the test fails. User should use mbuild.Box
            error_msg = (
                "If you are not providing an empty box, "
                "you need to specify the atoms/beads as children in the mb.Compound. "
                "If you are providing and empty box, please do so by specifying and "
                f"mbuild Box ({mb.Box})"
            )
            raise TypeError(error_msg)

        else:
            initial_no_atoms = len(structure.to_parmed().atoms)
            new_gmso_topology = mb_convert(structure, custom_groups=residues)
    else:
        initial_no_atoms = 0
        lengths = structure.lengths
        angles = structure.angles

        # create a new empty gmso topology.  This is needed because an empty mbuild.Compound
        # can not be converted to gmso.Topology without counting at 1 atom
        new_gmso_topology = gmso.Topology()
        new_gmso_topology.box = gmso.Box(lengths=lengths, angles=angles)

    return new_gmso_topology, initial_no_atoms


def _validate_forcefield_selection_and_residues(forcefield_selection, residues):
    """Validate if input forcefield_selection and reisudes is of the correct form."""
    if forcefield_selection is None:
        error_msg = (
            "Please the force field selection (forcefield_selection) as a dictionary "
            "with all the residues specified to a force field "
            '-> Ex: {"Water": "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package."
        )
        raise TypeError(error_msg)

    elif not isinstance(forcefield_selection, dict):
        error_msg = (
            "The force field selection (forcefield_selection) "
            "is not a dictionary. Please enter a dictionary "
            "with all the residues specified to a force field "
            '-> Ex: {"Water": "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package."
        )
        raise TypeError(error_msg)

    if not isinstance(residues, (list, tuple)):
        error_msg = (
            "Please enter the residues list in the specific_ff_to_residue."
        )
        raise TypeError(error_msg)

    forcefield_keys_list = list(forcefield_selection.keys())
    if forcefield_keys_list == [] and len(residues) != 0:
        print_error_message = "The forcefield_selection variable are not provided, but there are residues provided."
        raise ValueError(print_error_message)

    elif forcefield_keys_list != [] and len(residues) == 0:
        print_error_message = (
            "The residues variable is an empty list but there are "
            "forcefield_selection variables provided."
        )
        raise ValueError(print_error_message)


def _validate_boxes_for_simulation(boxes_for_simulation):
    """Validate if input boxes_for_simulation is of the correct form."""
    if boxes_for_simulation not in [1, 2]:
        boxes_for_simulation_error_msg = (
            "boxes_for_simulation must be either 1 or 2."
        )
        raise ValueError(boxes_for_simulation_error_msg)


def _validate_forcefields(forcefield_selection, residues):
    """Validate and create GMSO ForceField object from the forcefield_selection."""
    if has_foyer:
        from foyer import Forcefield
        from foyer.forcefields import forcefields

    forcefield_keys_list = list(forcefield_selection.keys())
    user_entered_ff_with_path_dict = {}
    # True means user entered the path, False is a standard foyer FF with no path
    for residue in residues:
        if residue in forcefield_keys_list:
            ff_extension = os.path.splitext(forcefield_selection[residue])[1]
            if ff_extension == ".xml":
                user_entered_ff_with_path_dict[residue] = True
            elif ff_extension == "":
                user_entered_ff_with_path_dict[residue] = False
            else:
                error_msg = "Please make sure you are enterning the correct FF name or path with xml extension"
                # "Please make sure you are entering the correct foyer FF name or a path to a FF file (with .xml extension)."
                raise ValueError(error_msg)

    # check if FF files exist and create a forcefield selection with directory paths
    # forcefield_selection_with_paths
    forcefield_selection_with_paths = {}
    for residue in forcefield_keys_list:
        ff_for_residue = forcefield_selection[residue]
        if user_entered_ff_with_path_dict[residue]:
            ff_names_path_iteration = forcefield_selection[residue]
            try:
                read_xlm_iteration = minidom.parse(ff_names_path_iteration)
                forcefield_selection_with_paths[residue] = (
                    ff_names_path_iteration
                )

            except:
                error_msg = (
                    "Please make sure you are entering the correct foyer FF path, "
                    "including the FF file name.xml. "
                    "If you are using the pre-build FF files in foyer, "
                    "only use the string name without any extension. "
                    "The selected FF file could also could not formated properly, or "
                    "there may be errors in the FF file itself."
                )
                raise ValueError(error_msg)
        elif not user_entered_ff_with_path_dict[residue]:
            ff_for_residue = forcefield_selection[residue]
            ff_names_path_iteration = (
                f"{forcefields.get_ff_path()[0]}/xml/{ff_for_residue}.xml"
            )
            try:
                read_xlm_iteration = minidom.parse(ff_names_path_iteration)
                forcefield_selection_with_paths[residue] = (
                    ff_names_path_iteration
                )
            except:
                error_msg = (
                    "Please make sure you are entering the correct foyer FF name, or the "
                    "correct file extension (i.e., .xml, if required)."
                )
                raise ValueError(error_msg)

    # push the FF paths and/or name to the GMSO format and create the new GMSO topology format
    gmso_compatable_forcefield_selection = {}
    for ff_key_iter, ff_value_iter in forcefield_selection_with_paths.items():
        # try to load the Foyer and GMSO FFs, if Foyer convert to GMSO; otherwise, it is an error.
        try:
            try:
                ff_new_gmso_value_iter = FoyerFFs.get_ff(
                    ff_value_iter
                ).to_gmso_ff()
            except:
                ff_new_gmso_value_iter = GMSOFFs.get_ff(
                    ff_value_iter
                ).to_gmso_ff()

        except:
            error_msg = (
                f"The supplied force field xml for the "
                f"{ff_key_iter} residue is not a foyer or gmso xml, "
                f"or the xml has errors and it not able to load properly."
            )
            raise TypeError(error_msg)

        gmso_compatable_forcefield_selection.update(
            {ff_key_iter: ff_new_gmso_value_iter}
        )
    return gmso_compatable_forcefield_selection


def _cross_check_residues_and_unique_site_labels(
    structure, residues, unique_topology_groups_list, boxes_for_simulation
):
    """Cross checking the residues list and unique_site_labels list."""
    # Check that all residues in box are in the residue names
    if isinstance(structure, mb.Compound):
        for applied_res_i in unique_topology_groups_list:
            if applied_res_i not in residues:
                error_msg_all_res_not_specified = (
                    f"All the residues are not specified in the residue list, or "
                    f"the {applied_res_i} residue does not match the residues that "
                    f"were found in the foyer and GMSO force field application. "
                )
                raise ValueError(error_msg_all_res_not_specified)

    # check if all the molecules/residues were found in in the mb.Compound/allowable input
    msg2 = (
        "All the residues were not used from the forcefield_selection "
        "string or dictionary. There may be residues below other "
        "specified residues in the mbuild.Compound hierarchy. "
        "If so, all the highest listed residues pass down the force "
        "fields through the hierarchy. Alternatively, residues that "
        "are not in the structure may have been specified. "
    )
    msg3 = (
        f"NOTE: This warning will appear if you are using the CHARMM pdb and psf writers "
        f"2 boxes, and the boxes do not contain all the residues in each box."
    )
    for res_i in residues:
        if res_i not in unique_topology_groups_list:
            msg1 = f"The {res_i} residues were not used from the forcefield_selection string or dictionary. "
            if boxes_for_simulation == 1:
                raise ValueError(f"{msg1}{msg2}")
            if boxes_for_simulation == 2:
                warn(f"{msg1}{msg2}{msg3}")
