import pytest
import unyt as u
from sympy import sympify

from gmso.core.forcefield import ForceField
from gmso.core.topology import Topology
from gmso.exceptions import NotYetImplementedWarning
from gmso.parameterization.parameterize import apply
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn


class TestTIP4PGMSO(ParameterizationBaseTest):
    def test_tip4p_files(self):
        mol2_file = get_path("tip3p.mol2")
        gmso_top = Topology.load(mol2_file)
        ff = ForceField(get_fn("gmso_xmls/test_ffstyles/tip4p_2005.xml"))
        gmso_top = apply(
            gmso_top,
            ff,
            speedup_by_molgraph=False,
            identify_connections=True,
        )
        assert len(gmso_top.virtual_sites) == 1
        vtype = gmso_top.virtual_sites[0].virtual_type
        assert ("HW", "OW", "HW") == vtype.member_classes
        assert vtype.virtual_potential.expression == sympify(
            "4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        )
        assert vtype.virtual_potential.parameters["epsilon"] == 0.0 * u.kJ / u.mol
        assert vtype.virtual_potential.parameters["sigma"] == 0.0 * u.nm
        assert vtype.virtual_position.expression == sympify(
            "ri + b*(rj-ri+a*(rk-rj))/norm(rj-ri+a*(rk-rj))"
        )
        assert vtype.virtual_position.parameters["a"] == 0.5 * u.dimensionless
        assert vtype.virtual_position.parameters["b"] == 0.15 * u.dimensionless

    def test_default_virtual_site_apply(
        self, basic_virtual_site_ff, empty_virtual_site_top
    ):
        gmso_top = apply(
            empty_virtual_site_top,
            basic_virtual_site_ff,
            speedup_by_molgraph=False,
            identify_connections=True,
        )
        assert len(gmso_top.virtual_sites) == 1
        vtype = gmso_top.virtual_sites[0].virtual_type
        assert vtype.member_classes is None
        assert ("c1", "c2", "c2", "c2", "c1") == vtype.member_types
        assert vtype.virtual_potential.expression == sympify("a+b")
        assert vtype.virtual_potential.parameters["a"] == 1.0 * u.kJ / u.mol
        assert vtype.virtual_potential.parameters["b"] == 1.0 * u.kJ / u.mol
        assert vtype.virtual_position.expression == sympify("ri+rj+rk+rl+rm+c+d")
        assert u.allclose_units(
            vtype.virtual_position.parameters["c"], [1, 1, 1] * u.nm
        )
        assert u.allclose_units(
            vtype.virtual_position.parameters["d"], [0.1, 0.1, 0.1] * u.nm
        )
        assert gmso_top.sites == gmso_top.virtual_sites[0].parent_sites

    def missing_vsite_position_function(
        self, basic_virtual_site_ff, empty_virtual_site_top
    ):
        gmso_top = apply(
            empty_virtual_site_top,
            basic_virtual_site_ff,
            speedup_by_molgraph=False,
            identify_connections=True,
        )
        assert len(gmso_top.virtual_sites) == 1
        with pytest.raises(NotYetImplementedWarning):
            gmso_top.virtual_sites[0].position()
