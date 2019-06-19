import json
import numpy as np

import unyt as u

import topology as topo
from topology.formats.hoomd_metadata import write_metadata
from topology.tests.base_test import BaseTest
from topology.utils.testing import allclose

class TestHoomdMetaData(BaseTest):
    def test_lj(self):
        top = topo.Topology()
        atype = topo.AtomType(name='A', 
                parameters={'sigma':2 * u.nm, 'epsilon':2 * u.Unit('kJ/mol')})
        top.add_site(topo.Site(name='a',
            atom_type=atype))
        write_metadata(top, 'out.json')
        meta = json.load(open('out.json', 'r'))
        assert np.isclose(meta['objects'][1]['hoomd.md.pair.lj']['tracked_fields']['parameters']['A,A']['sigma'], 2.0)
        assert np.isclose(meta['objects'][1]['hoomd.md.pair.lj']['tracked_fields']['parameters']['A,A']['epsilon'], 2.0)

    def test_lj_mix(self):
        top = topo.Topology()
        top.combining_rule = 'lorentz'
        atype = topo.AtomType(name='A', 
                parameters={'sigma':2 * u.nm, 'epsilon':2 * u.Unit('kJ/mol')})

        btype = topo.AtomType(name='B', 
                parameters={'sigma':4 * u.nm, 'epsilon':4 * u.Unit('kJ/mol')})

        top.add_site(topo.Site(name='a',
            atom_type=atype))
        top.add_site(topo.Site(name='b',
            atom_type=btype))

        write_metadata(top, 'out.json')
        meta = json.load(open('out.json', 'r'))
        assert np.isclose(meta['objects'][1]['hoomd.md.pair.lj']['tracked_fields']['parameters']['A,B']['sigma'], 3)
        assert np.isclose(meta['objects'][1]['hoomd.md.pair.lj']['tracked_fields']['parameters']['A,B']['epsilon'], (2*4)**0.5)

    def test_harmonic_bond(self):
        top = topo.Topology()
        top.combining_rule = 'lorentz'
        atype = topo.AtomType(name='A', 
                parameters={'sigma':2 * u.nm, 'epsilon':2 * u.Unit('kJ/mol')})

        btype = topo.AtomType(name='B', 
                parameters={'sigma':4 * u.nm, 'epsilon':4 * u.Unit('kJ/mol')})

        asite = topo.Site(name='a',
            atom_type=atype)
        bsite = topo.Site(name='b',
            atom_type=btype)


        bondtype = topo.BondType(parameters={'r_eq': 0.14 * u.nm,
                                            'k': 3.0 * u.Unit('kJ/mol/nm**2')},
                                member_types=['A', 'B'])
        bond = topo.Bond(connection_members=[asite, bsite], 
                connection_type=bondtype)
        top.add_connection(bond)

        top.update_top()
        write_metadata(top, 'out.json')
        meta = json.load(open('out.json', 'r'))

        assert np.isclose(meta['objects'][1]['hoomd.md.bond.harmonic']['tracked_fields']['parameters']['A-B']['r0'], 0.14)
        assert np.isclose(meta['objects'][1]['hoomd.md.bond.harmonic']['tracked_fields']['parameters']['A-B']['k'], 3.0)

    def test_harmonic_angle(self):
        top = topo.Topology()
        top.combining_rule = 'lorentz'
        atype = topo.AtomType(name='A', 
                parameters={'sigma':2 * u.nm, 'epsilon':2 * u.Unit('kJ/mol')})
        btype = topo.AtomType(name='B', 
                parameters={'sigma':4 * u.nm, 'epsilon':4 * u.Unit('kJ/mol')})
        ctype = topo.AtomType(name='C', 
                parameters={'sigma':4 * u.nm, 'epsilon':4 * u.Unit('kJ/mol')})


        asite = topo.Site(name='a',
            atom_type=atype)
        bsite = topo.Site(name='b',
            atom_type=btype)
        csite = topo.Site(name='c',
            atom_type=ctype)


        angletype = topo.AngleType(parameters={'theta_eq': 3.14 * u.rad,
                                            'k': 3.0 * u.Unit('kJ/mol/rad**2')},
                                member_types=['A', 'B', 'C'])
        angle = topo.Angle(connection_members=[asite, bsite, csite], 
                connection_type=angletype)
        top.add_connection(angle)

        top.update_top()
        write_metadata(top, 'out.json')
        meta = json.load(open('out.json', 'r'))

        assert np.isclose(meta['objects'][1]['hoomd.md.angle.harmonic']['tracked_fields']['parameters']['A-B-C']['t0'], 3.14)
        assert np.isclose(meta['objects'][1]['hoomd.md.angle.harmonic']['tracked_fields']['parameters']['A-B-C']['k'], 3.0)
