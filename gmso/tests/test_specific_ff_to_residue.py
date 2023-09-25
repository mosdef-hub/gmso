import mbuild as mb
import pytest
from foyer.forcefields import forcefields
from mbuild import Box, Compound
from mbuild.utils.io import has_foyer

from gmso.exceptions import GMSOError
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn
from gmso.utils.specific_ff_to_residue import specific_ff_to_residue


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestSpecificFFToResidue(BaseTest):
    # Tests for the mbuild.utils.specific_FF_to_residue.Specific_FF_to_residue() function
    def test_specific_ff_ff_is_none(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"Please the force field selection \(forcefield_selection\) as a "
            r"dictionary with all the residues specified to a force field "
            '-> Ex: {"Water": "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package.",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection=None,
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_ff_wrong_ff_extension(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match="Please make sure you are enterning the correct FF name or path with xml extension",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection={ethane_gomc.name: "oplsaa.pdb"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_all_residue_not_input(self, ethane_gomc, ethanol_gomc):
        with pytest.raises(
            GMSOError,
            match=f"A particle named C cannot be associated with the\n        "
            f"custom_groups \['ETH'\]. "
            f"Be sure to specify a list of group names that will cover\n        "
            f"all particles in the compound. This particle is one level below ETO.",
        ):
            box = mb.fill_box(
                compound=[ethane_gomc, ethanol_gomc],
                box=[1, 1, 1],
                n_compounds=[1, 1],
            )

            specific_ff_to_residue(
                box,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=2,
            )

    def test_specific_ff_to_residue_ff_selection_not_dict(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"The force field selection \(forcefield_selection\) "
            "is not a dictionary. Please enter a dictionary "
            "with all the residues specified to a force field "
            '-> Ex: {"Water": "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package.",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection="oplsaa",
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_is_none(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"Please enter the residues list in the specific_ff_to_residue.",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=None,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_simulation_boxes_not_1_or_2(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"boxes_for_simulation must be either 1 or 2",
        ):
            test_box_ethane_gomc = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[2, 3, 4]
            )

            specific_ff_to_residue(
                test_box_ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=3,
            )

    def test_specific_ff_to_residue_ffselection_wrong_path(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            # match=r"FileNotFoundError: \[Errno 2\] No such file or directory: 'oplsaa.xml'"
            match=r"Please make sure you are entering the correct foyer FF path, "
            r"including the FF file name.xml. "
            r"If you are using the pre-build FF files in foyer, "
            r"only use the string name without any extension. "
            r"The selected FF file could also could not formated properly, "
            r"or there may be errors in the FF file itself.",
        ):
            test_box_ethane_gomc = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 5, 6]
            )

            specific_ff_to_residue(
                test_box_ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa.xml"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_ff_wrong_path(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"Please make sure you are entering the correct foyer FF path, "
            r"including the FF file name.xml. "
            r"If you are using the pre-build FF files in foyer, "
            r"only use the string name without any extension. "
            r"The selected FF file could also could not formated properly, "
            r"or there may be errors in the FF file itself.",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection={ethane_gomc.name: "oplsaa.xml"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_input_string_as_compound(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"The structure expected to be of type: "
            r"<class 'mbuild.compound.Compound'> or <class 'mbuild.box.Box'>, "
            r"received: <class 'str'>",
        ):
            specific_ff_to_residue(
                "ethane_gomc",
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_boxes_for_simulation_not_int(
        self, ethane_gomc
    ):
        with pytest.raises(
            ValueError, match=r"boxes_for_simulation must be either 1 or 2."
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=1.1,
            )

    def test_specific_ff_to_residues_no_ff(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match="The forcefield_selection variable are not provided, but there are residues provided.",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection={},
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_no_residues(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"The residues variable is an empty list but there are "
            "forcefield_selection variables provided.",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[],
                boxes_for_simulation=1,
            )

    def test_specific_ff_wrong_foyer_name(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"Please make sure you are entering the correct foyer FF name, "
            r"or the correct file extension \(i.e., .xml, if required\).",
        ):
            box_0 = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
            )

            specific_ff_to_residue(
                box_0,
                forcefield_selection={ethane_gomc.name: "xxx"},
                residues=[ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_ff_selection_run(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[4, 5, 6]
        )

        [
            test_topology,
            test_residues_applied_list,
            test_electrostatics14Scale_dict,
            test_nonBonded14Scale_dict,
            test_atom_types_dict,
            test_bond_types_dict,
            test_angle_types_dict,
            test_dihedral_types_dict,
            test_improper_types_dict,
            test_combining_rule_dict,
        ] = specific_ff_to_residue(
            test_box_ethane_gomc,
            forcefield_selection={
                ethane_gomc.name: f"{forcefields.get_ff_path()[0]}/xml/oplsaa.xml"
            },
            residues=[ethane_gomc.name],
            boxes_for_simulation=1,
        )
        assert test_electrostatics14Scale_dict == {"ETH": 0.5}
        assert test_nonBonded14Scale_dict == {"ETH": 0.5}
        assert test_residues_applied_list == ["ETH"]

    def test_specific_ff_to_no_atoms_no_box_dims_in_residue(self):
        with pytest.raises(
            TypeError,
            match=f"The structure, {mb.Compound} or {mb.Box}, needs to have have box lengths and angles.",
        ):
            empty_compound = mb.Compound()

            specific_ff_to_residue(
                empty_compound,
                forcefield_selection={"empty_compound": "oplsaa"},
                residues=["empty_compound"],
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_no_atoms_in_residue(self):
        with pytest.raises(
            ValueError,
            match=r"The residues variable is an empty list but there "
            r"are forcefield_selection variables provided.",
        ):
            empty_compound = mb.Compound()
            empty_compound.box = Box([4, 4, 4])

            specific_ff_to_residue(
                empty_compound,
                forcefield_selection={"empty_compound": "oplsaa"},
                residues=[],
                boxes_for_simulation=1,
            )

    # this test is not usable until the individual mb.Compounds can be force fielded in MosDeF-GOMC

    def test_charmm_empty_compound_test_no_children(self, methane_ua_gomc):
        empty_box = mb.Compound()
        empty_box.box = mb.Box(lengths=[4, 4, 4])

        with pytest.raises(
            TypeError,
            match=r"If you are not providing an empty box, "
            r"you need to specify the atoms/beads as children in the mb.Compound. "
            r"If you are providing and empty box, please do so by specifying and "
            r"mbuild Box \({}\)".format(type(Box(lengths=[1, 1, 1]))),
        ):
            specific_ff_to_residue(
                empty_box,
                forcefield_selection={"AAA": "trappe-ua"},
                residues=["AAA"],
                boxes_for_simulation=1,
            )

    def test_charmm_a_few_mbuild_layers(self, ethane_gomc, ethanol_gomc):
        box_reservior_1 = mb.fill_box(
            compound=[ethane_gomc], box=[1, 1, 1], n_compounds=[1]
        )
        box_reservior_1.name = ethane_gomc.name
        box_reservior_1.periodicity = (True, True, True)
        box_reservior_2 = mb.fill_box(
            compound=[ethanol_gomc], box=[1, 1, 1], n_compounds=[1]
        )
        box_reservior_2.name = ethanol_gomc.name
        box_reservior_2.translate([0, 0, 1])

        box_reservior_3 = mb.Compound()
        box_reservior_3.name = "A"
        box_reservior_3.box = Box(lengths=[3, 3, 3])
        box_reservior_3.add(box_reservior_1, inherit_periodicity=False)
        box_reservior_3.add(box_reservior_2, inherit_periodicity=False)

        [
            test_topology,
            test_residues_applied_list,
            test_electrostatics14Scale_dict,
            test_nonBonded14Scale_dict,
            test_atom_types_dict,
            test_bond_types_dict,
            test_angle_types_dict,
            test_dihedral_types_dict,
            test_improper_types_dict,
            test_combining_rule_dict,
        ] = specific_ff_to_residue(
            box_reservior_3,
            forcefield_selection={
                ethanol_gomc.name: "oplsaa",
                ethane_gomc.name: "oplsaa",
                # box_reservior_3.name: "oplsaa",
            },
            residues=[ethanol_gomc.name, ethane_gomc.name],
            # residues=[box_reservior_3.name],
            boxes_for_simulation=1,
        )

        assert test_topology.n_sites == 17
        assert test_electrostatics14Scale_dict == {"ETO": 0.5, "ETH": 0.5}
        assert test_nonBonded14Scale_dict == {"ETO": 0.5, "ETH": 0.5}
        assert test_residues_applied_list.sort() == ["ETO", "ETH"].sort()

    def test_charmm_all_residues_not_in_dict_boxes_for_simulation_1(
        self, ethane_gomc, ethanol_gomc
    ):
        with pytest.raises(
            ValueError,
            match=f"The {'ETO'} residues were not used from the forcefield_selection string or dictionary. "
            "All the residues were not used from the forcefield_selection "
            "string or dictionary. There may be residues below other "
            "specified residues in the mbuild.Compound hierarchy. "
            "If so, all the highest listed residues pass down the force "
            "fields through the hierarchy. Alternatively, residues that "
            "are not in the structure may have been specified. ",
        ):
            box_reservior_0 = mb.fill_box(
                compound=[ethane_gomc], box=[1, 1, 1], n_compounds=[1]
            )
            specific_ff_to_residue(
                box_reservior_0,
                forcefield_selection={
                    ethanol_gomc.name: "oplsaa",
                    ethane_gomc.name: "oplsaa",
                },
                residues=[ethanol_gomc.name, ethane_gomc.name],
                boxes_for_simulation=1,
            )

    def test_charmm_all_residues_not_in_dict_boxes_for_simulation_2(
        self, ethane_gomc, ethanol_gomc
    ):
        with pytest.warns(
            UserWarning,
            match=f"The {'ETO'} residues were not used from the forcefield_selection string or dictionary. "
            "All the residues were not used from the forcefield_selection "
            "string or dictionary. There may be residues below other "
            "specified residues in the mbuild.Compound hierarchy. "
            "If so, all the highest listed residues pass down the force "
            "fields through the hierarchy. Alternatively, residues that "
            "are not in the structure may have been specified. "
            f"NOTE: This warning will appear if you are using the CHARMM pdb and psf writers "
            f"2 boxes, and the boxes do not contain all the residues in each box.",
        ):
            box_reservior_0 = mb.fill_box(
                compound=[ethane_gomc], box=[1, 1, 1], n_compounds=[1]
            )
            specific_ff_to_residue(
                box_reservior_0,
                forcefield_selection={
                    ethanol_gomc.name: "oplsaa",
                    ethane_gomc.name: "oplsaa",
                },
                residues=[ethanol_gomc.name, ethane_gomc.name],
                boxes_for_simulation=2,
            )

    def test_specific_ff_params_benzene_water_aa(self):
        benzene_aa = mb.load(get_fn("benzene_aa.mol2"))
        benzene_aa.name = "BEN"
        water_aa = mb.load("O", smiles=True)
        water_aa.name = "WAT"
        test_box = mb.fill_box(
            compound=[benzene_aa, water_aa], n_compounds=[1, 1], box=[4, 4, 4]
        )

        [
            test_topology,
            test_residues_applied_list,
            test_electrostatics14Scale_dict,
            test_nonBonded14Scale_dict,
            test_atom_types_dict,
            test_bond_types_dict,
            test_angle_types_dict,
            test_dihedral_types_dict,
            test_improper_types_dict,
            test_combining_rule_dict,
        ] = specific_ff_to_residue(
            test_box,
            forcefield_selection={
                benzene_aa.name: f"{get_fn('gmso_xmls/test_ffstyles/benzene_GAFF.xml')}",
                water_aa.name: f"{get_fn('gmso_xmls/test_ffstyles/spce_water__lorentz_combining.xml')}",
            },
            residues=[benzene_aa.name, water_aa.name],
            boxes_for_simulation=1,
        )
        assert test_topology.n_sites == 15
        assert test_topology.n_bonds == 14
        assert test_topology.n_angles == 19
        assert test_topology.n_dihedrals == 24
        assert test_topology.n_impropers == 6

        assert test_electrostatics14Scale_dict == {"BEN": 0.833333333, "WAT": 0}
        assert test_nonBonded14Scale_dict == {"BEN": 0.5, "WAT": 0}
        assert test_residues_applied_list == ["BEN", "WAT"]
        assert test_combining_rule_dict == "lorentz"

        # atom tests
        assert (
            str(test_atom_types_dict["BEN"]["expression"])
            == "4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        )
        assert (
            len(list(test_atom_types_dict["BEN"]["atom_types"].yield_view()))
            == 2
        )

        assert (
            str(test_atom_types_dict["WAT"]["expression"])
            == "4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        )
        assert (
            len(list(test_atom_types_dict["WAT"]["atom_types"].yield_view()))
            == 2
        )

        # bond tests
        assert (
            str(test_bond_types_dict["BEN"]["expression"])
            == "k*(r - r_eq)**2/2"
        )
        assert (
            len(list(test_bond_types_dict["BEN"]["bond_types"].yield_view()))
            == 2
        )

        assert (
            str(test_bond_types_dict["WAT"]["expression"]) == "k*(r - r_eq)**2"
        )
        assert (
            len(list(test_bond_types_dict["WAT"]["bond_types"].yield_view()))
            == 1
        )

        # angle tests
        assert (
            str(test_angle_types_dict["BEN"]["expression"])
            == "k*(theta - theta_eq)**2/2"
        )
        assert (
            len(list(test_angle_types_dict["BEN"]["angle_types"].yield_view()))
            == 2
        )

        assert (
            str(test_angle_types_dict["WAT"]["expression"])
            == "k*(theta - theta_eq)**2"
        )
        assert (
            len(list(test_angle_types_dict["WAT"]["angle_types"].yield_view()))
            == 1
        )

        # dihedral tests
        assert (
            str(test_dihedral_types_dict["BEN"]["expression"])
            == "k*(cos(n*phi - phi_eq) + 1)"
        )
        assert (
            len(
                list(
                    test_dihedral_types_dict["BEN"][
                        "dihedral_types"
                    ].yield_view()
                )
            )
            == 3
        )

        # improper tests
        assert (
            str(test_improper_types_dict["BEN"]["expression"])
            == "k*(cos(n*phi - phi_eq) + 1)"
        )
        assert (
            len(
                list(
                    test_improper_types_dict["BEN"][
                        "improper_types"
                    ].yield_view()
                )
            )
            == 1
        )

    def test_specific_ff_params_benzene_aa_grouped(self):
        methane_ua_bead_name = "_CH4"
        methane_child_bead = mb.Compound(name=methane_ua_bead_name)
        methane_box = mb.fill_box(
            compound=methane_child_bead, n_compounds=4, box=[1, 2, 3]
        )
        methane_box.name = "MET"

        [
            test_topology,
            test_residues_applied_list,
            test_electrostatics14Scale_dict,
            test_nonBonded14Scale_dict,
            test_atom_types_dict,
            test_bond_types_dict,
            test_angle_types_dict,
            test_dihedral_types_dict,
            test_improper_types_dict,
            test_combining_rule_dict,
        ] = specific_ff_to_residue(
            methane_box,
            forcefield_selection={
                methane_box.name: "trappe-ua",
            },
            residues=[methane_box.name],
            gmso_match_ff_by="group",
            boxes_for_simulation=1,
        )
        assert test_topology.n_sites == 4
        assert test_topology.n_bonds == 0
        assert test_topology.n_angles == 0
        assert test_topology.n_dihedrals == 0
        assert test_topology.n_impropers == 0

        assert test_electrostatics14Scale_dict == {"MET": 0}
        assert test_nonBonded14Scale_dict == {"MET": 0}
        assert test_residues_applied_list == ["MET"]
        assert test_combining_rule_dict == "lorentz"

        # atom tests
        assert (
            str(test_atom_types_dict["MET"]["expression"])
            == "4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        )
        assert (
            len(list(test_atom_types_dict["MET"]["atom_types"].yield_view()))
            == 1
        )
