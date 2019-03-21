import numpy as np
import unyt as u
import warnings

from topology.core.connection import Connection
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
        self._site_list = list()
        self._connection_list = list()
        self._bond_list = list()
        self._angle_list = list()

        self._atom_types = list()
        self._connection_types = list()
        self._bond_types = list()
        self._angle_types = list()

        self._atom_type_expressions = list()
        self._connection_type_expressions = list()
        self._bond_type_expressions = list()
        self._angle_type_expressions = list()

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

    def positions(self):
        xyz = np.empty(shape=(self.n_sites, 3)) * u.nm
        for i, site in enumerate(self.site_list):
            xyz[i, :] = site.position
        return xyz

    def add_site(self, site):
        if site in self.site_list:
            warnings.warn("Redundantly adding Site {}".format(site))
        self._site_list.append(site)
        self.update_atom_types()

    def add_connection(self, connection):
        if connection in self.connection_list:
            warnings.warn("Redundantly adding Connection {}".format(connection))

        if connection.connection_members[0] not in self.site_list:
            self.add_site(connection.connection_members[0])
        if connection.connection_members[1] not in self.site_list:
            self.add_site(connection.connection_members[1])

        self._connection_list.append(connection)

        #self.update_connection_list() Do we need to call this? Code should work either way
        self.update_connection_types()
        if isinstance(connection, Bond):
            self.update_bond_list()
            self.update_bond_types()
        elif isinstance(connection, Angle):
            self.update_angle_list()
            self.update_angle_types()

    @property
    def n_sites(self):
        return len(self.site_list)

    @property
    def n_connections(self):
        return len(self.connection_list)

    @property
    def n_bonds(self):
        return len(self.bond_list)

    @property
    def n_angles(self):
        return len(self.angle_list)

    @property
    def site_list(self):
        return self._site_list

    @property
    def connection_list(self):
        return self._connection_list

    @property
    def bond_list(self):
        return self._bond_list

    @property
    def angle_list(self):
        return self._angle_list

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

    def update_top(self):
        """ Update the entire topology's attributes

        Notes
        -----
        Will update: sites, connections, bonds, angles, 
        atom_types, connectiontypes, bondtypes, angletypes
        """
        self.update_site_list()
        self.update_connection_list()
        self.update_bond_list()
        self.update_angle_list()

        self.update_atom_types()
        self.update_connection_types()
        self.update_bond_types()
        self.update_angle_types()

    def update_site_list(self):
        """ (Is this necessary?) 
        Update site list based on the connection members """
        for connection in self.connection_list:
            for con_member in connection.connection_members:
                if con_member not in self.site_list:
                    self.add_site(con_member)
            
    def update_connection_list(self):
        """ Update connection list based on the site list """
        #self._connection_list = []
        for site in self.site_list:
            for connection in site.connections:
                if connection not in self.connection_list:
                    self.add_connection(connection)

    def update_bond_list(self):
        """ Rebuild the bond list by filtering through connection list """
        self._bond_list = [b for b in self.connection_list if isinstance(b, Bond)]

    def update_angle_list(self):
        """ Rebuild the angle list by filtering through connection list """
        self._angle_list = [a for a in self.connection_list if isinstance(a, Angle)]

    def update_atom_types(self):
        """ Update the atom types based on the site list """
        #self._atom_types = []
        for site in self.site_list:
            if site.atom_type is None:
                warnings.warn("Site {} detected with no AtomType".format(site))
            elif site.atom_type not in self.atom_types:
                self.atom_types.append(site.atom_type)

    def update_connection_types(self):
        """ Update the connection types based on the connection list """
        #self._connection_types = []
        for c in self.connection_list:
            if c.connection_type is None:
                warnings.warn("Non-parametrized Connection {} detected".format(c))
            elif not isinstance(c.connection_type, Potential):
                raise TopologyError("Non-Potential {} found "
                        "in Connection {}".format(c.connection_type, c))
            elif c.connection_type not in self.connection_types:
                self.connection_types.append(c.connection_type)

    def update_bond_types(self):
        """ Update the bond types based on the bond list """
        #self._bond_types = []
        for b in self.bond_list:
            if b.connection_type is None:
                warnings.warn("Non-parametrized Bond {} detected".format(b))
            elif not isinstance(b.connection_type, BondType):
                raise TopologyError("Non-BondType {} found in Bond {}".format(
                    b.connection_type, b))
            elif b.connection_type not in self.bond_types:
                self.bond_types.append(b.connection_type)

    def update_angle_types(self):
        """ Update the angle types based on the angle list """
        #self._angle_types = []
        for a in self.angle_list:
            if a.connection_type is None:
                warnings.warn("Non-parametrized Angle {} detected".format(a))
            elif not isinstance(a.connection_type, AngleType):
                raise TopologyError("Non-AngleType {} found in Angle {}".format(
                    a.connection_type, a))
            elif a.connection_type not in self.angle_types:
                self.angle_types.append(a.connection_type)

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ' ')
        descr.append('{:d} sites, '.format(self.n_sites))
        descr.append('{:d} connections, '.format(self.n_connections))
        descr.append('id: {}>'.format(id(self)))

        return ''.join(descr)
