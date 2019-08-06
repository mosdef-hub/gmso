import warnings

import numpy as np
import unyt as u
from boltons.setutils import IndexedSet

from topology.core.bond import Bond
from topology.core.angle import Angle
from topology.core.potential import Potential
from topology.core.bond_type import BondType
from topology.core.angle_type import AngleType
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
        self._box = box
        self._sites = IndexedSet()
        self._connections = IndexedSet()
        self._bonds = IndexedSet()
        self._angles = IndexedSet()

        self._subtops = IndexedSet()

        self._typed = False
        self._atom_types = IndexedSet()
        self._connection_types = IndexedSet()
        self._bond_types = IndexedSet()
        self._angle_types = IndexedSet()

        self._combining_rule = 'lorentz'

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

    def is_typed(self):
        self.update_atom_types()
        self.update_connection_types()
        self.update_bond_types()
        self.update_angle_types()

        if len(self.atom_types) > 0 or len(self.connection_types) > 0:
            self._typed = True
        else:
            self._typed = False
        return self.typed

    @property
    def combining_rule(self):
        return self._combining_rule

    @combining_rule.setter
    def combining_rule(self, rule):
        if rule not in ['lorentz', 'geometric']:
            raise TopologyError('Combining rule must be `lorentz` or `geometric`')
        self._combining_rule = rule

    def positions(self):
        xyz = np.empty(shape=(self.n_sites, 3)) * u.nm
        for i, site in enumerate(self.sites):
            xyz[i, :] = site.position
        return xyz

    def add_site(self, site, update_types=True):
        # Might be a more elegant way of handling this, see PR #128
        if site in self.sites:
            warnings.warn("Redundantly adding Site {}".format(site))
        self._sites.add(site)
        if update_types:
            if not self.typed:
                self.typed = True
            self.update_atom_types()

    def add_connection(self, connection, update_types=True):
        # Might be a more elegant way of handling this, see PR #128
        if connection in self.connections:
            warnings.warn("Redundantly adding Connection {}".format(connection))

        for conn_member in connection.connection_members:
            if conn_member not in self.sites:
                self.add_site(conn_member)

        if update_types and not self.typed:
            self.typed = True

        self._connections.add(connection)
        if isinstance(connection, Bond):
            self.update_bonds()
        elif isinstance(connection, Angle):
            self.update_angles()
        self.update_connections()

        if update_types:
            if not self.typed:
                self.typed = True
            if isinstance(connection, Bond):
                self.update_bond_types()
            elif isinstance(connection, Angle):
                self.update_angle_types()
            self.update_connection_types()

    def add_subtopology(self, subtop):
        self._subtops.add(subtop)
        subtop.parent = self
        # Note: would remove duplicates but there should be none
        self._sites.union(subtop.sites)

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

    def enumerate_connectivity(self):
        """ Enumerate connectivity within a topology 
        
        Assumes bonds are provided, will infer angles, dihedrals, impropers"""
        from topology.core.connectivity import identify_connections
        identify_connections(self)

    def update_top(self, update_types=True):
        """ Update the entire topology's attributes

        Notes
        -----
        Will update: sites, connections, bonds, angles,
        atom_types, connectiontypes, bondtypes, angletypes
        """
        self.update_sites()
        self.update_connections()
        self.update_bonds()
        self.update_angles()

        if update_types:
            self.update_atom_types()
            self.update_connection_types()
            self.update_bond_types()
            self.update_angle_types()
            self.is_typed()

    def update_sites(self):
        """ (Is this necessary?)
        Update site list based on the connection members """
        for connection in self.connections:
            for con_member in connection.connection_members:
                if con_member not in self.sites:
                    self.add_site(con_member)

    def update_connections(self):
        """ Update connection list based on the site list """
        #self._connections = []
        for site in self.sites:
            for connection in site.connections:
                if connection not in self.connections:
                    self.add_connection(connection)

    def update_bonds(self):
        """ Rebuild the bond list by filtering through connection list """
        self._bonds = [b for b in self.connections if isinstance(b, Bond)]

    def update_angles(self):
        """ Rebuild the angle list by filtering through connection list """
        self._angles = [a for a in self.connections if isinstance(a, Angle)]

    def update_atom_types(self):
        """ Update the atom types based on the site list """
        #self._atom_types = []
        for site in self.sites:
            if site.atom_type is None:
                warnings.warn("Site {} detected with no AtomType".format(site))
            elif site.atom_type not in self.atom_types:
                self.atom_types.add(site.atom_type)

    def update_connection_types(self):
        """ Update the connection types based on the connection list """
        #self._connection_types = []
        for c in self.connections:
            if c.connection_type is None:
                warnings.warn("Non-parametrized Connection {} detected".format(c))
            elif not isinstance(c.connection_type, Potential):
                raise TopologyError("Non-Potential {} found "
                        "in Connection {}".format(c.connection_type, c))
            elif c.connection_type not in self.connection_types:
                self.connection_types.add(c.connection_type)

    def update_bond_types(self):
        """ Update the bond types based on the bond list """
        #self._bond_types = []
        for b in self.bonds:
            if b.connection_type is None:
                warnings.warn("Non-parametrized Bond {} detected".format(b))
            elif not isinstance(b.connection_type, BondType):
                raise TopologyError("Non-BondType {} found in Bond {}".format(
                    b.connection_type, b))
            elif b.connection_type not in self.bond_types:
                self.bond_types.add(b.connection_type)

    def update_angle_types(self):
        """ Update the angle types based on the angle list """
        #self._angle_types = []
        for a in self.angles:
            if a.connection_type is None:
                warnings.warn("Non-parametrized Angle {} detected".format(a))
            elif not isinstance(a.connection_type, AngleType):
                raise TopologyError("Non-AngleType {} found in Angle {}".format(
                    a.connection_type, a))
            elif a.connection_type not in self.angle_types:
                self.angle_types.add(a.connection_type)

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ' ')
        descr.append('{:d} sites, '.format(self.n_sites))
        descr.append('{:d} connections, '.format(self.n_connections))
        descr.append('id: {}>'.format(id(self)))

        return ''.join(descr)

    def __eq__(self, other):
        """Compare a topology for equivalence."""

        if self is other:
            return True

        if not isinstance(other, Topology):
            return False

        if self.name != other.name:
            return False

        if self.n_sites != other.n_sites:
            return False

        if self.combining_rule != other.combining_rule:
            return False

        for (con1, con2) in zip(self.connections, other.connections):
            if con1 != con2:
                return False

        for (site1, site2) in zip(self.sites, other.sites):
            if site1 != site2:
                return False

        if self.box != other.box:
            return False

        return True
