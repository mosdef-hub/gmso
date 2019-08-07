import pytest

from topology.core.connection import Connection
from topology.core.site import Site
from topology.tests.base_test import BaseTest


class TestConnectivity(BaseTest):
    def test_methane_connectivity(self, methane):
        assert methane.n_bonds == 4
        assert methane.n_angles == 0
        methane.enumerate_connectivity()
        assert methane.n_bonds == 4
        assert methane.n_angles == 6

    def test_ethane_connectivity(self, ethane):
        assert ethane.n_bonds == 7
        assert ethane.n_angles == 0
        ethane.enumerate_connectivity()
        assert ethane.n_bonds == 7
        assert ethane.n_angles == 12
