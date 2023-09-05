import glob
import os

import pytest

from gmso.core.forcefield import ForceField
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn

# Make source directory for all xmls to grab from
XML_DIR = get_fn("gmso_xmls")
TEST_XMLS = glob.glob(os.path.join(XML_DIR, "*/*.xml"))


def compare_xml_files(fn1, fn2):
    """Hash files to check for lossless conversion."""
    with open(fn1, "r") as f:
        line1 = f.readlines()
    with open(fn2, "r") as f:
        line2 = f.readlines()
    for l1, l2 in zip(line1, line2):
        assert l1.replace(" ", "") == l2.replace(" ", "")
    return True


class TestXMLHandling(BaseTest):
    @pytest.fixture
    def ff(self):
        return ForceField(get_path("ff-example0.xml"))

    @pytest.fixture
    def named_groups_ff(self):
        return ForceField(get_path("ff-example1.xml"))

    @pytest.fixture
    def opls_ethane_foyer(self):
        return ForceField(
            get_path(filename=get_path("oplsaa-ethane_foyer.xml"))
        )

    def test_write_xml(self, opls_ethane_foyer):
        opls_ethane_foyer.to_xml("test_xml_writer.xml")
        reloaded_xml = ForceField("test_xml_writer.xml")
        get_names = lambda ff, param: [
            typed for typed in getattr(ff, param).keys()
        ]
        for param in [
            "atom_types",
            "bond_types",
            "angle_types",
            "dihedral_types",
        ]:
            assert get_names(opls_ethane_foyer, param) == get_names(
                reloaded_xml, param
            )

    def test_foyer_xml_conversion(self):
        """Validate xml converted from Foyer can be written out correctly."""
        pass

    def test_write_xml_from_topology(self):
        """Validate xml from a typed topology matches loaded xmls."""
        pass

    @pytest.mark.parametrize("xml", TEST_XMLS)
    def test_load__direct_from_forcefield_utilities(self, xml):
        """Validate loaded xmls from ff-utils match original file."""
        ff = ForceField(xml)
        assert isinstance(ff, ForceField)

    @pytest.mark.parametrize("xml", TEST_XMLS)
    def test_ffutils_backend(self, xml):
        ff1 = ForceField(xml)
        assert isinstance(ff1, ForceField)
        ff2 = ForceField(xml, backend="gmso", strict=False)
        assert isinstance(ff2, ForceField)
        assert ff1 == ff2

    @pytest.mark.parametrize("xml", TEST_XMLS)
    def test_gmso_backend(self, xml):
        ff = ForceField(xml, backend="gmso", strict=False)
        assert isinstance(ff, ForceField)

    @pytest.mark.parametrize("xml", TEST_XMLS)
    def test_load_write_xmls_gmso_backend(self, xml):
        """Validate loaded xmls written out match original file."""
        ff1 = ForceField(xml, backend="forcefield_utilities")
        ff1.to_xml("tmp.xml", overwrite=True)
        ff2 = ForceField("tmp.xml", strict=False)
        if "test_ffstyles" not in xml:
            assert compare_xml_files(xml, "tmp.xml")
        assert ff1 == ff2

    @pytest.mark.parametrize("xml", TEST_XMLS)
    def test_load_write_xmls_ffutils_backend(self, xml):
        """Validate loaded xmls written out match original file."""
        ff1 = ForceField(xml, backend="forcefield-utilities")
        ff1.to_xml("tmp.xml", overwrite=True)
        ff2 = ForceField("tmp.xml")
        if "test_ffstyles" not in xml:
            assert compare_xml_files("tmp.xml", xml)
        assert ff1 == ff2

    def test_xml_error_handling(self):
        """Validate bad xml formatting in xmls."""
        file_path = "dummy_name.xml"
        with pytest.raises(FileNotFoundError):
            ForceField(file_path)
        with pytest.raises(IndexError):
            ForceField(get_path("empty_foyer.xml"))

    def test_kb_in_ffutils(self):
        xml_path = get_path("ff-example0.xml")
        ff = ForceField(xml_path)
        assert ff
