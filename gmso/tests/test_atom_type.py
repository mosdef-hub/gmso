from copy import deepcopy

import pytest
import sympy
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.topology import Topology
from gmso.tests.base_test import BaseTest


class TestAtomType(BaseTest):
    @pytest.fixture(scope="session")
    def atomtype_metadata(self):
        return AtomType()

    def test_atom_type_tag_kwarg(self):
        at = AtomType(
            tags={"element": "Li", "comesFrom": "ForceFieldExperiments"}
        )
        assert at.tag_names == ["element", "comesFrom"]

    def test_new_atom_type(self, charge, mass):
        new_type = AtomType(
            name="mytype",
            charge=charge,
            mass=mass,
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            parameters={
                "sigma": 1 * u.nm,
                "epsilon": 10 * u.Unit("kcal / mol"),
            },
            independent_variables={"r"},
        )
        assert new_type.name == "mytype"
        assert_allclose_units(new_type.charge, charge, rtol=1e-5, atol=1e-8)
        assert_allclose_units(
            new_type.parameters["sigma"], 1 * u.nm, rtol=1e-5, atol=1e-8
        )
        assert_allclose_units(
            new_type.parameters["epsilon"],
            10 * u.Unit("kcal / mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(new_type.mass, mass, rtol=1e-5, atol=1e-8)

    def test_setters(self, charge, mass):
        new_type = AtomType()
        new_type.name = "SettingName"
        new_type.charge = -1.0 * charge
        new_type.mass = 1 * mass
        new_type.independent_variables = "r"
        new_type.parameters = {
            "sigma": 1 * u.nm,
            "epsilon": 10 * u.Unit("kcal / mol"),
        }
        new_type.expression = "r * sigma * epsilon"
        assert new_type.name == "SettingName"
        assert_allclose_units(
            new_type.charge, -1.0 * charge, rtol=1e-5, atol=1e-8
        )
        assert_allclose_units(new_type.mass, 1 * mass, rtol=1e-5, atol=1e-8)
        assert new_type.independent_variables == {sympy.symbols("r")}
        assert new_type.parameters == {
            "sigma": 1 * u.nm,
            "epsilon": 10 * u.Unit("kcal / mol"),
        }
        assert new_type.expression == sympy.sympify("r * sigma * epsilon")

    def test_incorrect_indep_vars(self):
        with pytest.raises(ValueError):
            AtomType(expression="x*y", independent_variables="z")

    def test_incorrect_expression(self, charge):
        with pytest.raises(ValueError):
            AtomType(name="mytype", charge=charge, expression=4.2)

    def test_expression_consistency(self, charge):
        # Test nb-func symbol consistency with parameter consistency in init
        new_type = AtomType(
            name="mytype",
            charge=charge,
            parameters={"x": 1 * u.m, "y": 10 * u.m},
            expression="x + y * z",
            independent_variables="z",
        )

        symbol_x, symbol_y, symbol_z = sympy.symbols("x y z")
        correct_expr = sympy.sympify("x+y*z")
        assert new_type.expression.free_symbols == set(
            [symbol_x, symbol_y, symbol_z]
        )
        assert correct_expr == new_type.expression

    def test_equivalance(self, charge):
        first_type = AtomType(
            name="mytype",
            charge=charge,
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.m},
            independent_variables={"r"},
        )
        same_type = AtomType(
            name="mytype",
            charge=charge,
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.m},
            independent_variables={"r"},
        )
        different_name = AtomType(
            name="difftype",
            charge=charge,
        )
        different_charge = AtomType(
            name="mytype",
            charge=4.0 * charge,
        )
        different_function = AtomType(
            name="mytype",
            charge=charge,
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.m},
            expression="r * sigma * epsilon",
            independent_variables={"r"},
        )
        different_params = AtomType(
            name="mytype",
            charge=charge,
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            parameters={"sigma": 42 * u.m, "epsilon": 100000 * u.m},
            independent_variables={"r"},
        )
        different_mass = AtomType(
            name="mytype",
            charge=charge,
            mass=5 * u.kg / u.mol,
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.m},
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            independent_variables={"r"},
        )

        assert first_type == same_type
        assert first_type != different_name
        assert first_type != different_charge
        assert first_type != different_function
        assert first_type != different_params
        assert first_type != different_mass

    def test_set_nb_func(self, charge):
        # Try changing the nonbonded function, but keep the parameters
        first_type = AtomType(
            name="mytype",
            charge=charge,
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.m},
            independent_variables="r",
        )
        first_type.set_expression(expression="r * (sigma + epsilon)")
        correct_expr = sympy.sympify("r * (sigma + epsilon)")
        assert first_type.expression == correct_expr
        assert first_type.parameters == {"sigma": 1 * u.m, "epsilon": 10 * u.m}

    def test_set_nb_func_bad(self):
        # Try changing the nonbonded function, keep the parameters,
        # but the nb function uses different symbols
        first_type = AtomType(
            expression="r*sigma*epsilon",
            parameters={"sigma": 1 * u.m, "epsilon": 100 * u.m},
            independent_variables="r",
        )
        with pytest.raises(ValueError):
            first_type.set_expression(expression="a + b")

    def test_set_nb_func_params(self, charge):
        # Try changing the nb parameters, but keeping the function
        first_type = AtomType(
            name="mytype",
            charge=charge,
            expression="r * sigma * epsilon",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.m},
            independent_variables="r",
        )
        first_type.set_expression(
            parameters={"sigma": 42 * u.m, "epsilon": 24 * u.m}
        )
        correct_expr = sympy.sympify("r * sigma * epsilon")
        assert first_type.expression == correct_expr
        assert first_type.parameters == {"sigma": 42, "epsilon": 24}

    def test_set_nb_func_params_bad(self):
        # Try changing the parameters, keep the function,
        # but the new parameters use different symbols
        first_type = AtomType(
            expression="r*sigma*epsilon",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.J},
            independent_variables={"r"},
        )

        with pytest.raises(ValueError):
            first_type.set_expression(parameters={"a": 1 * u.g, "b": 10 * u.m})

    def test_set_nb_func_params_partial(self):
        # Try changing the parameters, keep the function,
        # but change only one symbol
        first_type = AtomType(
            expression="r*sigma*epsilon",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.J},
            independent_variables={"r"},
        )
        first_type.set_expression(parameters={"sigma": 42 * u.m})
        correct_expr = sympy.sympify("r*sigma*epsilon")
        assert first_type.parameters == {"sigma": 42 * u.m, "epsilon": 10 * u.J}
        assert first_type.expression == correct_expr

    def test_set_nb_func_params_both_correct(self):
        # Try correctly changing both the nb function and parameters
        first_type = AtomType(
            expression="r*sigma*epsilon",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.J},
            independent_variables="r",
        )
        first_type.set_expression(
            expression="a+b*x",
            parameters={"a": 100 * u.m, "b": 42 * u.J},
            independent_variables="x",
        )
        correct_expr = sympy.sympify("a+b*x")
        correct_params = {"a": 100 * u.m, "b": 42 * u.J}
        assert first_type.expression == correct_expr
        assert first_type.parameters == correct_params

    def test_set_nb_func_params_both_incorrect(self):
        # Try incorrectly changing both the nb function and the parameters
        first_type = AtomType(
            expression="r*sigma*epsilon",
            parameters={"sigma": 1 * u.m, "epsilon": 10 * u.J},
            independent_variables="r",
        )
        with pytest.raises(ValueError):
            first_type.set_expression(
                expression="a*x+b",
                parameters={"c": 100 * u.year, "d": 42 * u.newton},
                independent_variables="x",
            )

    def test_metadata(self):
        valid_type = AtomType(
            doi="123",
            definition="[c]",
            overrides={"122"},
            description="Some type solely for testing",
        )
        assert valid_type.doi
        assert valid_type.overrides

        with pytest.raises(ValueError):
            bad_doi = AtomType(
                doi=123,
                definition="[c]",
                overrides={"122"},
                description="Some type solely for testing",
            )
            bad_defn = AtomType(
                doi="123",
                definition=123,
                overrides={"122"},
                description="Some type solely for testing",
            )
            bad_over = AtomType(
                doi="123",
                definition="[c]",
                overrides="122",
                description="Some type solely for testing",
            )
            bad_desc = AtomType(
                doi="123", definition="[c]", overrides="122", description=123
            )
            valid_type.doi = 123
            valid_type.definition = 123
            valid_type.overrides = ["123"]
            valid_type.description = 123

    def test_atom_type_with_topology_and_site(self):
        site1 = Atom()
        site2 = Atom()
        top = Topology()
        atom_type1 = AtomType()
        atom_type2 = AtomType()
        site1.atom_type = atom_type1
        site2.atom_type = atom_type2
        top.add_site(site1)
        top.add_site(site2)
        assert id(site1.atom_type) != id(site2.atom_type)
        assert site1.atom_type is not None
        assert len(top.atom_types) == 2

    def test_atom_type_with_topology_and_site_change_properties(self):
        site1 = Atom()
        site2 = Atom()
        top = Topology()
        atom_type1 = AtomType()
        atom_type2 = AtomType()
        site1.atom_type = atom_type1
        site2.atom_type = atom_type2
        top.add_site(site1)
        top.add_site(site2)
        site1.atom_type.mass = 250
        assert site1.atom_type.mass == 250
        assert next(iter(top.atom_types)).mass == 250

    def test_with_1000_atom_types(self):
        top = Topology()
        for i in range(1000):
            site = Atom()
            atom_type = AtomType()
            site.atom_type = atom_type
            top.add_site(site, update_types=False)
        top.update_topology()
        assert len(top.atom_types) == 1000
        assert top.n_sites == 1000

    def test_atom_type_copy(self, typed_ethane):
        for atom_type in typed_ethane.atom_types:
            assert atom_type.copy(deep=True) == atom_type
            assert deepcopy(atom_type) == atom_type

    def test_metadata_empty_tags(self, atomtype_metadata):
        assert atomtype_metadata.tag_names == []
        assert list(atomtype_metadata.tag_names_iter) == []

    def test_metadata_add_tags(self, atomtype_metadata):
        atomtype_metadata.add_tag("tag1", dict([("tag_name_1", "value_1")]))
        atomtype_metadata.add_tag("tag2", dict([("tag_name_2", "value_2")]))
        atomtype_metadata.add_tag("int_tag", 1)
        assert len(atomtype_metadata.tag_names) == 3

    def test_metadata_add_tags_overwrite(self, atomtype_metadata):
        with pytest.raises(ValueError):
            atomtype_metadata.add_tag("tag2", "new_value", overwrite=False)
        atomtype_metadata.add_tag("tag2", "new_value", overwrite=True)
        assert atomtype_metadata.get_tag("tag2") == "new_value"
        assert len(atomtype_metadata.tag_names) == 3

    def test_metadata_get_tags(self, atomtype_metadata):
        assert atomtype_metadata.get_tag("tag1").get("tag_name_1") == "value_1"
        assert atomtype_metadata.get_tag("int_tag") == 1
        assert atomtype_metadata.get_tag("non_existent_tag") is None
        with pytest.raises(KeyError):
            atomtype_metadata.get_tag("non_existent_tag", throw=True)

    def test_metadata_all_tags(self, atomtype_metadata):
        assert "int_tag" in atomtype_metadata.tags

    def test_metadata_delete_tags(self, atomtype_metadata):
        with pytest.raises(KeyError):
            atomtype_metadata.delete_tag("non_existent_tag")
        assert atomtype_metadata.pop_tag("non_existent_tag") is None
        atomtype_metadata.delete_tag("int_tag")
        assert len(atomtype_metadata.tag_names) == 2

    def test_atom_type_dict(self):
        atype = AtomType()
        atype_dict = atype.model_dump(exclude={"potential_expression"})
        assert "potential_expression" not in atype_dict
        assert "charge" in atype_dict

    def test_atom_type_clone(self):
        top = Topology()
        atype = AtomType(
            name="ff_255",
            expression="a*x+b+c*y",
            independent_variables={"x", "y"},
            parameters={"a": 200 * u.g, "b": 300 * u.K, "c": 400 * u.J},
            mass=2.0 * u.g / u.mol,
            charge=2.0 * u.elementary_charge,
            atomclass="CX",
            overrides={"ff_234"},
            definition="CC-C",
            description="Dummy Description",
        )
        atype_clone = atype.clone()

        atom1 = Atom(name="1")
        atom2 = Atom(name="2")
        atom1.atom_type = atype
        atom2.atom_type = atype_clone

        top.add_site(atom1)
        top.add_site(atom2)
        top.update_topology()

        assert len(top.atom_types) == 2

        atype_dict = atype.model_dump(exclude={"topology", "set_ref"})
        atype_clone_dict = atype_clone.model_dump(
            exclude={"topology", "set_ref"}
        )

        for key, value in atype_dict.items():
            cloned = atype_clone_dict[key]
            assert value == cloned
            if id(value) == id(cloned):
                assert isinstance(value, str)
                assert isinstance(cloned, str)
