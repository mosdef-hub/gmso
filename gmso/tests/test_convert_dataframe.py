import numpy as np
import pytest
import unyt as u

from gmso.external.convert_dataframe import (
    _recursive_getattr,
    multi_topology_dataframe,
    to_dataframeDict,
)
from gmso.tests.base_test import BaseTest
from gmso.utils.io import has_pandas


class TestConvertDataFrame(BaseTest):
    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_recursive_sites(self, typed_ethane):
        out = list(_recursive_getattr(typed_ethane, "sites", "atom_type.atomclass"))
        expected = [site.atom_type.atomclass for site in typed_ethane.sites]
        assert out == expected

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_recursive_dihedrals(self, typed_ethane):
        out = list(
            _recursive_getattr(typed_ethane, "dihedrals", "dihedral_type.member_types")
        )
        expected = [
            dihedral.dihedral_type.member_types for dihedral in typed_ethane.dihedrals
        ]
        assert out == expected

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_to_dataframeDict(self, typed_ethane):
        expected_valuesList = [8, 7, 12, 9]
        checkList = ["sites", "bonds", "angles", "dihedrals"]
        for parameter, val in zip(checkList, expected_valuesList):
            assert (
                len(to_dataframeDict(typed_ethane, parameter=parameter)[parameter])
                == val
            )
        allDict = to_dataframeDict(typed_ethane, parameter="all")
        dfList = [allDict.get(key) for key in checkList]
        assert list(map(len, dfList)) == expected_valuesList

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_dataframe_impropers(self, benzeneTopology):
        expected_valuesList = [12, 12, 18, 24, 6]
        checkList = ["sites", "bonds", "angles", "dihedrals", "impropers"]
        for parameter, val in zip(checkList, expected_valuesList):
            assert (
                len(to_dataframeDict(benzeneTopology, parameter=parameter)[parameter])
                == val
            )
        allDict = to_dataframeDict(benzeneTopology, parameter="all")
        dfList = [allDict.get(key) for key in checkList]
        assert list(map(len, dfList)) == expected_valuesList

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_dataframe_default_columns(self, typed_ethane):
        expected_columns = [
            "name",
            "atom_type.name",
            "sigma",
            "epsilon",
            "charge",
            "mass",
        ]
        assert np.all(
            list(to_dataframeDict(typed_ethane, "sites")["sites"].columns)
            == expected_columns
        )

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_dataframe_specified_columns(self, typed_ethane):
        input_columns = ["name", "position", "group"]
        expected_columns = ["name", "x", "y", "z", "group"]
        assert np.all(
            list(
                to_dataframeDict(
                    typed_ethane,
                    "sites",
                    columns=input_columns,
                    format="specific_columns",
                )["sites"].columns
            )
            == expected_columns
        )

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_dataframe_publication(self, benzeneTopology):
        dfDict = to_dataframeDict(benzeneTopology, "all", format="publication")
        df = dfDict["sites"]
        assert len(df.index) == 2
        assert len(df.columns) == 8
        assert "Atom Indices" in df.columns
        assert df["Atom Indices"].loc[0] == ", ".join(str(v) for v in np.arange(6))
        assert df["Atom Indices"].loc[1] == ", ".join(str(v) for v in np.arange(6, 12))

        connectList = ["bonds", "angles", "dihedrals", "impropers"]
        checkList = [2, 2, 3, 1]
        for connect, check in zip(connectList, checkList):
            df = dfDict[connect]
            assert len(df.index) == check

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_dataframe_remove_duplicates(self, benzeneTopology):
        dfDict = to_dataframeDict(benzeneTopology, "sites", format="remove_duplicates")
        df = dfDict["sites"]
        assert len(df.index) == 2
        assert len(df.columns) == 1
        assert "Atom Indices" not in df.columns

        connectList = ["bonds", "angles", "dihedrals", "impropers"]
        checkList = [2, 2, 3, 1]
        for connect, check in zip(connectList, checkList):
            dfDict = to_dataframeDict(
                benzeneTopology, connect, format="remove_duplicates"
            )
            df = dfDict[connect]
            assert len(df.index) == check

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_dataframe_unyts(self, typed_ethane):
        dfDict = to_dataframeDict(
            typed_ethane, "all", format="publication", handle_unyts="with_data"
        )
        df = dfDict["sites"]
        assert isinstance(df["charge"].loc[0], u.unyt_quantity)

        dfDict = to_dataframeDict(
            typed_ethane, "all", format="publication", handle_unyts="all_floats"
        )
        df = dfDict["sites"]
        assert isinstance(df["charge"].loc[0], float)

    @pytest.mark.skipif(not has_pandas, reason="Pandas is not installed")
    def test_multi_topology_dataframe(self, benzeneTopology, spce_water):
        dfDict = multi_topology_dataframe(
            [benzeneTopology, spce_water, benzeneTopology]
        )
        connectList = ["sites", "bonds", "angles", "dihedrals", "impropers"]
        checkList = [4, 3, 3, 3, 1]
        for connect, check in zip(connectList, checkList):
            df = dfDict.get(connect)
            if df is None:
                assert df == check
            else:
                assert len(df.index) == check
