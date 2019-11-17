import warnings

import numpy as np
import unyt as u
from boltons.setutils import IndexedSet

from topology.core.bond import Bond
from topology.core.angle import Angle
from topology.core.dihedral import Dihedral
from topology.core.potential import Potential
from topology.core.atom_type import AtomType
from topology.core.bond_type import BondType
from topology.core.angle_type import AngleType
from topology.core.dihedral_type import DihedralType
from topology.exceptions import TopologyError


class Topology(object):
    """A topology.

    Parameters
    ----------
    name : str, optional
        A name for the Topology.
    """
    def __init__(self, name="Topology", box=None):
        if name is not None:
            self._name = name
        else:
            self._name = "Topology"

        self._box = box
        self._sites = IndexedSet()
        self._typed = False
        self._connections = IndexedSet()
        self._bonds = IndexedSet()
        self._angles = IndexedSet()
        self._dihedrals = IndexedSet()
        self._subtops = IndexedSet()
        self._atom_types = IndexedSet()
        self._connection_types = IndexedSet()
        self._bond_types = IndexedSet()
        self._angle_types = IndexedSet()
        self._dihedral_types = IndexedSet()
        self._combining_rule = 'lorentz'
        self._set_refs = {
            'atom_type_set': self._atom_types,
            'bond_type_set': self._bond_types,
            'angle_type_set': self._angle_types,
            'dihedral_type_set': self._dihedral_types,
        }

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = str(name)

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, box):
        self._box = box

    @property
    def typed(self):
        return self._typed

    @typed.setter
    def typed(self, typed):
        self._typed = typed

    @property
    def combining_rule(self):
        return self._combining_rule

    @combining_rule.setter
    def combining_rule(self, rule):
        if rule not in ['lorentz', 'geometric']:
            raise TopologyError('Combining rule must be `lorentz` or `geometric`')
        self._combining_rule = rule

    @property
    def positions(self):
        xyz = np.empty(shape=(self.n_sites, 3)) * u.nm
        for i, site in enumerate(self.sites):
            xyz[i, :] = site.position
        return xyz

    @property
    def n_sites(self):
        return len(self.sites)

    @property
    def n_connections(self):
        return len(self.connections)

    @property
    def n_bonds(self):
        return len(self.bonds)

    @property
    def n_angles(self):
        return len(self.angles)

    @property
    def n_dihedrals(self):
        return len(self.dihedrals)

    @property
    def subtops(self):
        return self._subtops

    @property
    def n_subtops(self):
        return len(self._subtops)

    @property
    def sites(self):
        return self._sites

    @property
    def connections(self):
        return self._connections

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles

    @property
    def dihedrals(self):
        return self._dihedrals

    @property
    def atom_types(self):
        return self._atom_types

    @property
    def connection_types(self):
        return self._connection_types

    @property
    def bond_types(self):
        return self._bond_types

    @property
    def angle_types(self):
        return self._angle_types

    @property
    def dihedral_types(self):
        return self._dihedral_types

    @property
    def atom_type_expressions(self):
        return list(set([atype.expression for atype in self.atom_types]))

    @property
    def connection_type_expressions(self):
        return list(set([contype.expression for contype in self.connection_types]))

    @property
    def bond_type_expressions(self):
        return list(set([btype.expression for btype in self.bond_types]))

    @property
    def angle_type_expressions(self):
        return list(set([atype.expression for atype in self.angle_types]))

    @property
    def dihedral_type_expressions(self):
        return list(set([atype.expression for atype in self.dihedral_types]))

    def add_site(self, site, update_types=False):
        """Add the site to the topology
        Parameters
        -----------
        site: (topology.core.Site), site to be added to this topology

        Returns
        -------
        None
        """
        self.sites.add(site)
        if update_types and site.atom_type:
            self.atom_types.add(site.atom_type)
            site.atom_type = self.atom_types[self.atom_types.index(site.atom_type)]
            site.atom_type.topology = self
            self.update_atom_types()

    def update_sites(self):
        for connection in self.connections:
            for member in connection.connection_members:
                if member not in self.sites:
                    self.add_site(member)

    def add_connection(self, connection, update=True):
        for conn_member in connection.connection_members:
            if conn_member not in self.sites:
                self.add_site(conn_member)
        self._connections.add(connection)
        if isinstance(connection, Bond):
            self._bonds.add(Bond)
        if isinstance(connection, Angle):
            self._angles.add(connection)
        if isinstance(connection, Dihedral):
            self._dihedrals.add(connection)
        if update:
            self.update_connection_types()

    def update_connections(self):
        for site in self.sites:
            for conn in site.connections:
                if conn not in self.connections:
                    self.add_connection(conn, update=False)
        self.update_connection_types()
        self.is_typed()

    def update_bonds(self):
        pass

    def update_angles(self):
        pass

    def update_dihedrals(self):
        pass

    def update_connection_types(self):
        """Update the connection types based on the connection set"""
        for c in self.connections:
            if c.connection_type is None:
                warnings.warn('Non-parametrized Connection {} detected'.format(c))
            elif not isinstance(c.connection_type, Potential):
                raise TopologyError('Non-Potential {} found'
                                    'in Connection {}'.format(c.connection_type, c))
            elif c.connection_type not in self.connection_types:
                c.connection_type.topology = self
                self.connection_types.add(c.connection_type)
                if isinstance(c.connection_type, BondType):
                    self._bond_types.add(c.connection_type)
                if isinstance(c.connection_type, AngleType):
                    self._angle_types.add(c.connection_type)
                if isinstance(c.connection_type, DihedralType):
                    self._dihedral_types.add(c.connection_type)

    def update_atom_types(self):
        """Update atom types in the topology"""
        if not self._typed:
            self._typed = True
        for site in self.sites:
            if site.atom_type is None:
                warnings.warn('Non-parametrized site detected {}'.format(site))
            elif not isinstance(site.atom_type, AtomType):
                raise TopologyError('Non AtomType instance found in site {}'.format(site))
            elif site.atom_type not in self.atom_types:
                site.atom_type.topology = self
                self.atom_types.add(site.atom_type)
                site.atom_type = self.atom_types[self.atom_types.index(site.atom_type)]

    def add_subtopology(self, subtop):
        self._subtops.add(subtop)
        subtop.parent = self
        self._sites.union(subtop.sites)

    def is_typed(self):
        self.update_connection_types()
        self.update_atom_types()

        if len(self.atom_types) > 0 or len(self.connection_types) > 0:
            self._typed = True
        else:
            self._typed = False
        return self.typed

    def update_angle_types(self):
        pass

    def update_bond_types(self):
        pass

    def update_dihedral_types(self):
        pass

    def update_topology(self):
        """Update the entire topology"""
        self.update_sites()
        self.update_connections()
        self.update_atom_types()
        self.update_connection_types()
        self.is_typed()

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ' ')
        descr.append('{:d} sites, '.format(self.n_sites))
        descr.append('{:d} connections, '.format(self.n_connections))
        descr.append('id: {}>'.format(id(self)))

        return ''.join(descr)

    update_angle_types = update_connection_types
    update_bond_types = update_connection_types
    update_dihedral_types = update_connection_types
    update_angles = update_bonds = update_angles = update_connections
