"""Base data structure for GMSO chemical systems."""
import warnings
from pathlib import Path

import numpy as np
import unyt as u
from boltons.setutils import IndexedSet

from gmso.abc.abstract_site import Site
from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.pairpotential_type import PairPotentialType
from gmso.core.parametric_potential import ParametricPotential
from gmso.exceptions import GMSOError
from gmso.utils._constants import (
    ANGLE_TYPE_DICT,
    ATOM_TYPE_DICT,
    BOND_TYPE_DICT,
    DIHEDRAL_TYPE_DICT,
    IMPROPER_TYPE_DICT,
    PAIRPOTENTIAL_TYPE_DICT,
)
from gmso.utils.connectivity import (
    identify_connections as _identify_connections,
)


class Topology(object):
    """A topology.

    A topology represents a chemical structure wherein lie the collection
    of sites which together form a chemical structure containing connections
    (gmso.Bond, gmso.Angle and gmso.Dihedral (along with their associated types).
    A topology is the fundamental data structure in GMSO, from which we can gather
    various information about the chemical structure and apply a forcefield
    before converting the structure into a format familiar to various simulation
    engines.

    Parameters
    ----------
    name : str, optional, default='Topology'
        A name for the Topology.
    box : gmso.Box, optional, default=None
        A gmso.Box object bounding the topology

    Attributes
    ----------
    typed : bool
        True if the topology is typed

    combining_rule : str, ['lorentz', 'geometric']
        The combining rule for the topology, can be either 'lorentz' or 'geometric'

    scaling_factors : dict
        A collection of scaling factors used in the forcefield

    n_sites : int
        Number of sites in the topology

    n_connections : int
        Number of connections in the topology (Bonds, Angles, Dihedrals, Impropers)

    n_bonds : int
        Number of bonds in the topology

    n_angles: int
        Number of angles in the topology

    n_dihedrals : int
        Number of dihedrals in the topology

    n_impropers : int
        Number of impropers in the topology

    n_subtops : int
        Number of subtopolgies in the topology

    connections : tuple of gmso.Connection objects
        A collection of bonds, angles, dihedrals, and impropers in the topology

    bonds : tuple of gmso.Bond objects
        A collection of bonds in the topology

    angles : tuple of gmso.Angle objects
        A collection of angles in the topology

    dihedrals : tuple of gmso.Dihedral objects
        A collection of dihedrals in the topology

    impropers : tuple of gmso.Improper objects
        A collection of impropers in the topology

    connection_types : tuple of gmso.Potential objects
        A collection of BondTypes, AngleTypes, DihedralTypes, and ImproperTypes in the topology

    atom_types : tuple of gmso.AtomType objects
        A collection of AtomTypes in the topology

    bond_types : tuple of gmso.BondType objects
        A collection of BondTypes in the topology

    angle_types : tuple of gmso.AngleType objects
        A collection go AngleTypes in the topology

    dihedral_types : tuple of gmso.DihedralType objects
        A collection of DihedralTypes in the topology

    improper_types : tuple of gmso.ImproperType objects
        A collection of ImproperTypes in the topology

    pairpotential_types : tuple of gmso.PairPotentialType objects
        A collection of PairPotentialTypes in the topology

    atom_type_expressions : list of gmso.AtomType.expression objects
        A collection of all the expressions for the AtomTypes in topology

    connection_type_expressions : list of gmso.Potential.expression objects
        A collection of all the expressions for the Potential objects in the topology that represent a connection type

    bond_type_expressions : list of gmso.BondType.expression objects
        A collection of all the expressions for the BondTypes in topology

    angle_type_expressions : list of gmso.AngleType.expression objects
        A collection of all the expressions for the AngleTypes in topology

    dihedral_type_expressions : list of gmso.DihedralType.expression objects
        A collection of all the expression for the DihedralTypes in the topology

    improper_type_expressions : list of gmso.ImproperType.expression objects
        A collection of all the expression for the ImproperTypes in the topology

    pairpotential_type_expressions : list of gmso.PairPotentialType.expression objects
        A collection of all the expression for the PairPotentialTypes in the topology

    See Also
    --------
    gmso.SubTopology :
        A topology within a topology
    """

    def __init__(self, name="Topology", box=None):

        self.name = name
        self._box = box
        self._sites = IndexedSet()
        self._typed = False
        self._connections = IndexedSet()
        self._bonds = IndexedSet()
        self._angles = IndexedSet()
        self._dihedrals = IndexedSet()
        self._impropers = IndexedSet()
        self._subtops = IndexedSet()
        self._atom_types = {}
        self._atom_types_idx = {}
        self._connection_types = {}
        self._bond_types = {}
        self._bond_types_idx = {}
        self._angle_types = {}
        self._angle_types_idx = {}
        self._dihedral_types = {}
        self._dihedral_types_idx = {}
        self._improper_types = {}
        self._improper_types_idx = {}
        self._combining_rule = "lorentz"
        self._pairpotential_types = {}
        self._pairpotential_types_idx = {}
        self._scaling_factors = {
            "nonBonded12Scale": 0.0,
            "nonBonded13Scale": 0.0,
            "nonBonded14Scale": 0.5,
            "electrostatics12Scale": 0.0,
            "electrostatics13Scale": 0.0,
            "electrostatics14Scale": 0.5,
        }
        self._set_refs = {
            ATOM_TYPE_DICT: self._atom_types,
            BOND_TYPE_DICT: self._bond_types,
            ANGLE_TYPE_DICT: self._angle_types,
            DIHEDRAL_TYPE_DICT: self._dihedral_types,
            IMPROPER_TYPE_DICT: self._improper_types,
            PAIRPOTENTIAL_TYPE_DICT: self._pairpotential_types,
        }

        self._index_refs = {
            ATOM_TYPE_DICT: self._atom_types_idx,
            BOND_TYPE_DICT: self._bond_types_idx,
            ANGLE_TYPE_DICT: self._angle_types_idx,
            DIHEDRAL_TYPE_DICT: self._dihedral_types_idx,
            IMPROPER_TYPE_DICT: self._improper_types_idx,
            PAIRPOTENTIAL_TYPE_DICT: self._pairpotential_types_idx,
        }

        self._unique_connections = {}

    @property
    def name(self):
        """Return the name of the topology."""
        return self._name

    @name.setter
    def name(self, name):
        """Set the name of the topology."""
        self._name = str(name) if name else "Topology"

    @property
    def box(self):
        """Return the Box of the topology."""
        return self._box

    @box.setter
    def box(self, box):
        """Set the box of the topology."""
        self._box = box

    @property
    def typed(self):
        """Check if the system is parametrized."""
        return self._typed

    @typed.setter
    def typed(self, typed):
        """Set the status of the topology if parametrized."""
        self._typed = typed

    @property
    def combining_rule(self):
        """Return the combining rule for the topology."""
        return self._combining_rule

    @combining_rule.setter
    def combining_rule(self, rule):
        """Set the combining rule for the topology."""
        if rule not in ["lorentz", "geometric"]:
            raise GMSOError("Combining rule must be `lorentz` or `geometric`")
        self._combining_rule = rule

    @property
    def scaling_factors(self):
        """Return the scaling factors for the topology."""
        return self._scaling_factors

    @scaling_factors.setter
    def scaling_factors(self, scaling_factors):
        """Set the scaling factors for the topology."""
        expected_items = [
            "nonBonded12Scale",
            "nonBonded13Scale",
            "nonBonded14Scale",
            "electrostatics12Scale",
            "electrostatics13Scale",
            "electrostatics14Scale",
        ]
        if not isinstance(scaling_factors, dict):
            raise GMSOError("Scaling factors should be a dictionary")
        for item in expected_items:
            if item not in scaling_factors.keys():
                raise GMSOError(
                    f"Expected {expected_items} as keys in the scaling factors"
                )
        for val in scaling_factors.values():
            if val < 0.0 or val > 1.0:
                raise GMSOError("Scaling factors should be between 0.0 and 1.0")

        self._scaling_factors = scaling_factors

    @property
    def positions(self):
        """Return the positions of the sites in the topology."""
        xyz = np.empty(shape=(self.n_sites, 3)) * u.nm
        for i, site in enumerate(self._sites):
            xyz[i, :] = site.position
        return xyz

    @property
    def n_sites(self):
        """Return the number of sites in the topology."""
        return len(self.sites)

    @property
    def n_connections(self):
        """Return the number of connections in the topology."""
        return len(self.connections)

    @property
    def n_bonds(self):
        """Return the number of bonds in the topology."""
        return len(self.bonds)

    @property
    def n_angles(self):
        """Return the amount of angles in the topology."""
        return len(self.angles)

    @property
    def n_dihedrals(self):
        """Return the amount of dihedrals in the topology."""
        return len(self.dihedrals)

    @property
    def n_impropers(self):
        """Return the number of impropers in the topology."""
        return len(self.impropers)

    @property
    def subtops(self):
        """Return the subtopologies in the topology."""
        return self._subtops

    @property
    def n_subtops(self):
        """Return number of subtopolgies."""
        return len(self._subtops)

    @property
    def sites(self):
        """Return all sites in the topology."""
        return tuple(self._sites)

    @property
    def connections(self):
        """Return all connections in topology."""
        return tuple(self._connections)

    @property
    def bonds(self):
        """Return all bonds in the topology."""
        return tuple(self._bonds)

    @property
    def angles(self):
        """Return all angles in the topology."""
        return tuple(self._angles)

    @property
    def dihedrals(self):
        """Return all dihedrals in the topology."""
        return tuple(self._dihedrals)

    @property
    def impropers(self):
        """Return all impropers in the topology."""
        return tuple(self._impropers)

    @property
    def atom_types(self):
        """Return all atom_types in the topology."""
        return tuple(self._atom_types.values())

    @property
    def connection_types(self):
        """Return all connection_types in the topology."""
        return tuple(self._connection_types.values())

    @property
    def bond_types(self):
        """Return all bond_types in the topology."""
        return tuple(self._bond_types.values())

    @property
    def angle_types(self):
        """Return all angle_types in the topology."""
        return tuple(self._angle_types.values())

    @property
    def dihedral_types(self):
        """Return all dihedral_types in the topology."""
        return tuple(self._dihedral_types.values())

    @property
    def improper_types(self):
        """Return all improper_types in the topology."""
        return tuple(self._improper_types.values())

    @property
    def pairpotential_types(self):
        return tuple(self._pairpotential_types.values())

    @property
    def atom_type_expressions(self):
        """Return all atom_type expressions in the topology."""
        return list(set([atype.expression for atype in self.atom_types]))

    @property
    def connection_type_expressions(self):
        """Return all connection_type expressions in the topology."""
        return list(
            set([contype.expression for contype in self.connection_types])
        )

    @property
    def bond_type_expressions(self):
        """Return all bond_type expressions in the topology."""
        return list(set([btype.expression for btype in self.bond_types]))

    @property
    def angle_type_expressions(self):
        """Return all angle_type expressions in the topology."""
        return list(set([atype.expression for atype in self.angle_types]))

    @property
    def dihedral_type_expressions(self):
        """Return all dihedral_type expressions in the topology."""
        return list(set([atype.expression for atype in self.dihedral_types]))

    @property
    def improper_type_expressions(self):
        """Return all improper_type expressions in the topology."""
        return list(set([atype.expression for atype in self.improper_types]))

    @property
    def pairpotential_type_expressions(self):
        return list(
            set([atype.expression for atype in self.pairpotential_types])
        )

    def add_site(self, site, update_types=True):
        """Add a site to the topology.

        This method will add a site to the existing topology, since
        sites are stored in an indexed set, adding redundant site
        will have no effect. If the update_types parameter is set to
        true (default behavior), this method will also check if there
        is an gmso.AtomType associated with the site and it to the
        topology's AtomTypes collection.

        Parameters
        ----------
        site : gmso.core.Site
            Site to be added to this topology
        update_types : (bool), default=True
            If true, add this site's atom type to the topology's set of AtomTypes
        """
        self._sites.add(site)
        if update_types and site.atom_type:
            site.atom_type.topology = self
            if site.atom_type in self._atom_types:
                site.atom_type = self._atom_types[site.atom_type]
            else:
                self._atom_types[site.atom_type] = site.atom_type
                self._atom_types_idx[site.atom_type] = len(self._atom_types) - 1
            self.is_typed(updated=False)

    def update_sites(self):
        """Update the sites of the topology.

        This method will update the sites in the topology
        based on the connection members, For example- if you
        add a bond to a topology, without adding the constituent
        sites, this method can be called to add the sites which are the
        connection members of the bond as shown below.

            >>> import gmso
            >>> site1 = gmso.Site(name='MySite1')
            >>> site2 = gmso.Site(name='MySite2')
            >>> bond1 = gmso.Bond(name='site1-site2', connection_members=[site1, site2])
            >>> this_topology = gmso.Topology('TwoSitesTopology')
            >>> this_topology.add_connection(bond1)
            >>> this_topology.update_sites()

        See Also
        --------
        gmso.Topology.add_site : Add a site to the topology.
        gmso.Topology.add_connection : Add a Bond, an Angle or a Dihedral to the topology.
        gmso.Topology.update_topology : Update the entire topology.
        """
        for connection in self.connections:
            for member in connection.connection_members:
                if member not in self._sites:
                    self.add_site(member)

    def add_connection(self, connection, update_types=True):
        """Add a gmso.Connection object to the topology.

        This method will add a gmso.Connection object to the
        topology, it can be used to generically include any
        Connection object i.e. Bond or Angle or Dihedral to
        the topology. According to the type of object added,
        the equivalent collection in the topology is updated.
        For example- If you add a Bond, this method will update
        topology.connections and topology.bonds object. Additionally,
        if update_types is True (default behavior), it will also
        update any Potential objects associated with the connection.

        Parameters
        ----------
        connection : one of gmso.Connection, gmso.Bond, gmso.Angle, gmso.Dihedral, or gmso.Improper object
        update_types : bool, default=True
            If True also add any Potential object associated with connection to the
            topology.

        Returns
        -------
        gmso.Connection
            The Connection object or equivalent Connection object that
            is in the topology
        """
        # Check if an equivalent connection is in the topology
        equivalent_members = connection._equivalent_members_hash()
        if equivalent_members in self._unique_connections:
            warnings.warn(
                "An equivalent connection already exists. "
                "Providing the existing equivalent Connection."
            )
            connection = self._unique_connections[equivalent_members]

        for conn_member in connection.connection_members:
            if conn_member not in self.sites:
                self.add_site(conn_member)
        self._connections.add(connection)
        self._unique_connections.update({equivalent_members: connection})
        if isinstance(connection, Bond):
            self._bonds.add(connection)
        if isinstance(connection, Angle):
            self._angles.add(connection)
        if isinstance(connection, Dihedral):
            self._dihedrals.add(connection)
        if isinstance(connection, Improper):
            self._impropers.add(connection)
        if update_types:
            self.update_connection_types()

        return connection

    def identify_connections(self):
        """Identify all connections in the topology."""
        _identify_connections(self)

    def update_connection_types(self):
        """Update the connection types based on the connection collection in the topology.

        This method looks into all the connection objects (Bonds, Angles, Dihedrals, Impropers) to
        check if any Potential object (BondType, AngleType, DihedralType, ImproperType) is not in the
        topology's respective collection and will add those objects there.

        See Also
        --------
        gmso.Topology.update_atom_types : Update atom types in the topology.
        """
        for c in self.connections:
            if c.connection_type is None:
                warnings.warn(
                    "Non-parametrized Connection {} detected".format(c)
                )
            elif not isinstance(c.connection_type, ParametricPotential):
                raise GMSOError(
                    "Non-Potential {} found"
                    "in Connection {}".format(c.connection_type, c)
                )
            elif c.connection_type not in self._connection_types:
                c.connection_type.topology = self
                self._connection_types[c.connection_type] = c.connection_type
                if isinstance(c.connection_type, BondType):
                    self._bond_types[c.connection_type] = c.connection_type
                    self._bond_types_idx[c.connection_type] = (
                        len(self._bond_types) - 1
                    )
                if isinstance(c.connection_type, AngleType):
                    self._angle_types[c.connection_type] = c.connection_type
                    self._angle_types_idx[c.connection_type] = (
                        len(self._angle_types) - 1
                    )
                if isinstance(c.connection_type, DihedralType):
                    self._dihedral_types[c.connection_type] = c.connection_type
                    self._dihedral_types_idx[c.connection_type] = (
                        len(self._dihedral_types) - 1
                    )
                if isinstance(c.connection_type, ImproperType):
                    self._improper_types[c.connection_type] = c.connection_type
                    self._improper_types_idx[c.connection_type] = (
                        len(self._improper_types) - 1
                    )
            elif c.connection_type in self.connection_types:
                if isinstance(c.connection_type, BondType):
                    c.connection_type = self._bond_types[c.connection_type]
                if isinstance(c.connection_type, AngleType):
                    c.connection_type = self._angle_types[c.connection_type]
                if isinstance(c.connection_type, DihedralType):
                    c.connection_type = self._dihedral_types[c.connection_type]
                if isinstance(c.connection_type, ImproperType):
                    c.connection_type = self._improper_types[c.connection_type]

    def add_pairpotentialtype(self, pairpotentialtype, update=True):
        """add a PairPotentialType to the topology

        This method checks whether the member_type of a PairPotentialType
        object is already stored in pairpotential_types. If so, update the
        pair potential between the member_type, and if not, add the
        PairPotentialType to the topology.

        Parameters
        ----------
        pairpotentialtype: gmso.core.PairPotentialType
            The PairPotentialType object to be added
        update: Boolean, default=True

        See Also:
        --------
        gmso.core.pairpotential_type: Pairwise potential that does not follow
        combination rules
        """
        if update:
            self.update_atom_types()
        if not isinstance(pairpotentialtype, PairPotentialType):
            raise GMSOError(
                "Non-PairPotentialType {} provided".format(pairpotentialtype)
            )
        for atype in pairpotentialtype.member_types:
            if atype not in [t.name for t in self.atom_types]:
                if atype not in [t.atomclass for t in self.atom_types]:
                    raise GMSOError(
                        "There is no name/atomclass of AtomType {} in current topology".format(
                            atype
                        )
                    )
        self._pairpotential_types[pairpotentialtype] = pairpotentialtype
        self._pairpotential_types_idx[pairpotentialtype] = (
            len(self._pairpotential_types) - 1
        )

    def remove_pairpotentialtype(self, pair_of_types):
        """Remove the custom pairwise potential between two AtomTypes/Atomclasses

        Parameters
        ----------
        pair_of_types: list-like of strs
            The pair (or set) of names or atomclasses of gmso.AtomTypes of which
            the custom pairwise potential should be removed
        """
        to_delete = []
        for t in self._pairpotential_types:
            if t.member_types == tuple(pair_of_types):
                to_delete.append(t)
        if len(to_delete) > 0:
            for t in to_delete:
                del self._pairpotential_types[t]
            self._reindex_connection_types(PAIRPOTENTIAL_TYPE_DICT)
        else:
            warnings.warn(
                "No pair potential specified for such pair of AtomTypes/atomclasses"
            )

    def update_atom_types(self):
        """Update atom types in the topology.

        This method checks all the sites in the topology which have an
        associated AtomType and if that AtomType is not in the topology's
        AtomTypes collection, it will add it there.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology
        """
        for site in self._sites:
            if site.atom_type is None:
                warnings.warn("Non-parametrized site detected {}".format(site))
            elif not isinstance(site.atom_type, AtomType):
                raise GMSOError(
                    "Non AtomType instance found in site {}".format(site)
                )
            elif site.atom_type not in self._atom_types:
                site.atom_type.topology = self
                self._atom_types[site.atom_type] = site.atom_type
                self._atom_types_idx[site.atom_type] = len(self._atom_types) - 1
            elif site.atom_type in self._atom_types:
                site.atom_type = self._atom_types[site.atom_type]
        self.is_typed(updated=True)

    def add_subtopology(self, subtop, update=True):
        """Add a sub-topology to this topology.

        This methods adds a gmso.Core.SubTopology object to the topology
        All the sites in this sub-topology are added to the collection of current
        sites in this topology.

        Parameters
        ----------
        subtop : gmso.SubTopology
            The sub-topology object to be added.
        update : bool, default=True

        See Also
        --------
        gmso.SubTopology : A topology within a topology
        """
        self._subtops.add(subtop)
        subtop.parent = self
        self._sites.union(subtop.sites)
        if update:
            self.update_topology()

    def is_typed(self, updated=False):
        """Verify if the topology is parametrized."""
        if not updated:
            self.update_connection_types()
            self.update_atom_types()

        if len(self.atom_types) > 0 or len(self.connection_types) > 0:
            self._typed = True
        else:
            self._typed = False
        return self._typed

    def is_fully_typed(self, updated=False, group="topology"):
        """Check if the topology or a specifc group of objects that make up the topology are fully typed

        Parameters
        ----------
        updated : bool, optional, default=False
            If False, will update the atom type and connection type list of the
            Topology before the checking step.
        group : str, optional, default='topology'
            Specific objects to be checked. Options include:
            'topology'  : check for status of all the topology constituents
            'sites'     : check for status of all topology.sites
            'bonds'     : check for status of all topology.bonds
            'angles'    : check for status of all topology.angles
            'dihedrals' : check for status of all topology.dihedrals
            'impropers' : check for status of all topology.impropers

        Returns
        -------
        bool
            Status of the check
        Notes
        -----
        `self._type` is set to True as long as the Topology is at least
        partially typed.
        """
        if not updated:
            self.update_connection_types()
            self.update_atom_types()

        typed_status = {
            "sites": lambda top: all(site.atom_type for site in top._sites),
            "bonds": lambda top: all(bond.bond_type for bond in top._bonds),
            "angles": lambda top: all(
                angle.angle_type for angle in top._angles
            ),
            "dihedrals": lambda top: all(
                dihedral.dihedral_type for dihedral in top._dihedrals
            ),
            "impropers": lambda top: all(
                improper.improper_type for improper in top._impropers
            ),
        }

        if group == "topology":
            result = list()
            for subgroup in typed_status:
                result.append(typed_status[subgroup](self))
            return all(result)
        elif group in typed_status:
            return typed_status[group](self)
        else:
            raise ValueError(
                f"Could not check typing status of {group}. "
                "Available options: 'topology', 'sites', 'bonds', "
                "'angles', 'dihedrals', 'impropers'."
            )

    def get_untyped(self, group):
        """Get the untyped (non-parametrized) objects of the Topology.

        Parameters
        ----------
        group : {'sites', 'bonds', 'angles', 'dihedrals', 'impropers', 'topology'}
            The group of objects to be checked. The 'topology' option will return
            all untyped object of the topology.

        Returns
        -------
        untyped : dict
            Dictionary of all untyped object, key of the dictionary corresponds to
            object group names define above.
        """
        untyped = dict()
        untyped_extractors = {
            "sites": self._get_untyped_sites,
            "bonds": self._get_untyped_bonds,
            "angles": self._get_untyped_angles,
            "dihedrals": self._get_untyped_dihedrals,
            "impropers": self._get_untyped_impropers,
        }
        if group == "topology":
            for subgroup in untyped_extractors:
                untyped.update(untyped_extractors[subgroup]())
        elif isinstance(group, (list, tuple, set)):
            for subgroup in group:
                untyped.update(untyped_extractors[subgroup]())
        elif isinstance(group, str) and group in untyped_extractors:
            untyped = untyped_extractors[group]()
        else:
            raise ValueError(
                f"Cannot get untyped {group}. "
                f"Available options: {[untyped_extractors.keys()]}."
            )
        return untyped

    def _get_untyped_sites(self):
        "Return a list of untyped sites"
        untyped = {"sites": list()}
        for site in self._sites:
            if not site.atom_type:
                untyped["sites"].append(site)
        return untyped

    def _get_untyped_bonds(self):
        "Return a list of untyped bonds"
        untyped = {"bonds": list()}
        for bond in self._bonds:
            if not bond.bond_type:
                untyped["bonds"].append(bond)
        return untyped

    def _get_untyped_angles(self):
        "Return a list of untyped angles"
        untyped = {"angles": list()}
        for angle in self._angles:
            if not angle.angle_type:
                untyped["angles"].append(angle)
        return untyped

    def _get_untyped_dihedrals(self):
        "Return a list of untyped dihedrals"
        untyped = {"dihedrals": list()}
        for dihedral in self._dihedrals:
            if not dihedral.dihedral_type:
                untyped["dihedrals"].append(dihedral)
        return untyped

    def _get_untyped_impropers(self):
        "Return a list of untyped impropers"
        untyped = {"impropers": list()}
        for improper in self._impropers:
            if not improper.improper_type:
                untyped["impropers"].append(improper)
        return untyped

    def update_angle_types(self):
        """Use gmso.Topology.update_connection_types to update AngleTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def update_bond_types(self):
        """Use gmso.Topology.update_connection_types to update BondTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def update_dihedral_types(self):
        """Use gmso.Topology.update_connection_types to update DihedralTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def update_improper_types(self):
        """Use gmso.Topology.update_connection_types to update ImproperTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def update_topology(self):
        """Update the entire topology."""
        self.update_sites()
        self.update_atom_types()
        self.update_connection_types()
        self.is_typed(updated=True)

    def _get_bonds_for(self, site):
        """Return a list of bonds in this Topology that the site is a part of."""
        bonds = []
        for bond in self.bonds:
            if site in bond.connection_members:
                bonds.append(bond)
        return bonds

    def _get_angles_for(self, site):
        """Return a list of angles in this Topology that the site is a part of."""
        angles = []
        for angle in self.angles:
            if site in angle.connection_members:
                angles.append(angle)
        return angles

    def _get_dihedrals_for(self, site):
        """Return a list of dihedrals in this Topology that the site is a part of."""
        dihedrals = []
        for dihedral in self.dihedrals:
            if site in dihedral.connection_members:
                dihedrals.append(dihedral)
        return dihedrals

    def get_index(self, member):
        """Get index of a member in the topology.

        Parameters
        ----------
        member : gmso Topology objects
            The member to for which to return index for.
            `member` can be of type gmso.Site, gmso.Bond, gmso.Angle, gmso.Dihedral, gmso.Improper,
            gmso.AtomType, gmso.BondType, gmso.AngleType, gmso.DihedralType or gmso.ImproperType.

        Returns
        -------
        int
            The index of the member in the topology's collection objects
        """
        refs = {
            Atom: self._sites,
            Bond: self._bonds,
            Angle: self._angles,
            Dihedral: self._dihedrals,
            Improper: self._impropers,
            AtomType: self._atom_types_idx,
            BondType: self._bond_types_idx,
            AngleType: self._angle_types_idx,
            DihedralType: self._dihedral_types_idx,
            ImproperType: self._improper_types_idx,
            PairPotentialType: self._pairpotential_types_idx,
        }

        member_type = type(member)

        if member_type not in refs.keys():
            raise TypeError(
                f"Cannot index member of type {member_type.__name__}"
            )

        try:
            index = refs[member_type].index(member)
        except AttributeError:
            index = refs[member_type][member]

        return index

    def _reindex_connection_types(self, ref):
        """Re-generate the indices of the connection types in the topology."""
        if ref not in self._index_refs:
            raise GMSOError(
                f"cannot reindex {ref}. It should be one of "
                f"{ANGLE_TYPE_DICT}, {BOND_TYPE_DICT}, "
                f"{ANGLE_TYPE_DICT}, {DIHEDRAL_TYPE_DICT}, {IMPROPER_TYPE_DICT},"
                f"{PAIRPOTENTIAL_TYPE_DICT}"
            )
        for i, ref_member in enumerate(self._set_refs[ref].keys()):
            self._index_refs[ref][ref_member] = i

    def iter_sites(self, key, value):
        """Iterate through this topology's sites based on certain attribute and their values.

        Parameters
        ----------
        key: str
            The attribute of the site to look for
        value:
            The value that the given attribute should be equal to

        Yields
        ------
        gmso.abc.abstract_site.Site
            The site where getattr(site, key) == value
        """
        if key not in Site.__iterable_attributes__:

            raise ValueError(
                f"`{key}` is not an iterable attribute for Site. "
                f"To check what the iterable attributes are see gmso.abc.abstract_site module."
            )

        if value is None:
            raise ValueError(
                "Expected `value` to be something other than None. Provided None."
            )

        for site in self.sites:
            if getattr(site, key) == value:
                yield site

    def iter_sites_by_residue_name(self, name):
        """Iterate through this topology's sites which contain this specific residue `name`.

        See Also
        --------
        gmso.core.topology.Topology.iter_sites
            The method to iterate over Topology's sites
        """
        return self.iter_sites("residue_name", name)

    def iter_sites_by_residue_number(self, number):
        """Iterate through this topology's sites which contain this specific residue `number`.

        See Also
        --------
        gmso.core.topology.Topology.iter_sites
            The method to iterate over Topology's sites
        """
        return self.iter_sites("residue_number", number)

    def save(self, filename, overwrite=False, **kwargs):
        """Save the topology to a file.

        Parameters
        ----------
        filename: str, pathlib.Path
            The file to save the topology as
        overwrite: bool, default=True
            If True, overwrite the existing file if it exists
        **kwargs:
            The arguments to specific file savers listed below(as extensions):
            * json: types, update, indent
        """
        if not isinstance(filename, Path):
            filename = Path(filename).resolve()

        if filename.exists() and not overwrite:
            raise FileExistsError(
                f"The file {filename} exists. Please set "
                f"overwrite=True if you wish to overwrite the existing file"
            )

        from gmso.formats import SaversRegistry

        saver = SaversRegistry.get_callable(filename.suffix)
        saver(self, filename, **kwargs)

    def __repr__(self):
        """Return custom format to represent topology."""
        return (
            f"<Topology {self.name}, {self.n_sites} sites,\n "
            f"{self.n_connections} connections,\n "
            f"{len(self.connection_types)} potentials,\n "
            f"id: {id(self)}>"
        )

    def __str__(self):
        """Return custom format to represent topology as a string."""
        return f"<Topology {self.name}, {self.n_sites} sites, id: {id(self)}>"

    @classmethod
    def load(cls, filename, **kwargs):
        """Load a file to a topology"""
        filename = Path(filename).resolve()
        from gmso.formats import LoadersRegistry

        loader = LoadersRegistry.get_callable(filename.suffix)
        return loader(filename, **kwargs)
