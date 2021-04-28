import pytest

from gmso.core.atom_type import AtomType
from gmso.core.atom import Atom
from gmso.core.topology import Topology
from gmso.tests.base_test import BaseTest
from gmso.formats.lammpsdata import write_lammpsdata
from gmso.formats.gro import write_gro
from gmso.formats.mcf import write_mcf
from gmso.formats.xyz import write_xyz
from gmso.formats.top import write_top
from gmso.external import from_mbuild

class TestPerformance(BaseTest):
    @pytest.mark.timeout(15)
    def test_lammps_performance(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=10000)
        write_lammpsdata(top, "data.ar")

    @pytest.mark.timeout(15)
    def test_mcf_performance(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=10000)
        write_mcf(top, "ar.mcf")

    @pytest.mark.timeout(15)
    def test_top_performance(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=10000)
        write_top(top, "ar.top")

    @pytest.mark.timeout(5)
    def test_xyz_performance(self, n_ar_system):
        cmpd = n_ar_system(n_sites=10000)
        top = from_mbuild(cmpd) 
        write_xyz(top, "ar.xyz")

    @pytest.mark.timeout(5)
    def test_gro_performance(self, n_ar_system):
        cmpd = n_ar_system(n_sites=10000)
        top = from_mbuild(cmpd) 
        write_gro(top, "ar.gro")
