import gmso
from gmso.core.box import Box
from gmso.formats.lammpsdata import write_lammpsdata
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class TestLammpsWriter(BaseTest):
    def test_write_lammps(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')

    def test_write_lammps_triclinic(self, topology_site):
        top = topology_site()
        top.box = Box(lengths=[1,1,1], angles=[60,90,120])
        write_lammpsdata(top, filename='data.triclinic')

    def test_water_lammps(self, water_system):
        top = water_system

        ff = gmso.ForceField(get_path('tip3p.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types[site.name]

        top.update_sites()
        top.update_atom_types()

        for bond in top.bonds:
            bond.bond_type = bond.connection_type = ff.bond_types['opls_111~opls_112']

        top.update_bonds()
        top.update_bond_types()

        for subtop in top.subtops:
            angle = gmso.core.angle.Angle(
                connection_members=[site for site in subtop.sites],
                name="opls_112~opls_111~opls_112",
                connection_type=ff.angle_types["opls_112~opls_111~opls_112"]
            )
            top.add_connection(angle)

        top.update_angles()
        top.update_angle_types()

        write_lammpsdata(top, 'data.water')
